module _rMSEn
export rMSEn
using Statistics: std, mean, median, var
using DSP.Filters: filtfilt, Butterworth, Lowpass, digitalfilter
using Plots
    
    """
        MSx, CI = rMSEn(Sig, Mobj) 

    Returns a vector of refined multiscale entropy values (`MSx`) and the complexity 
    index (`CI`) of the data sequence (`Sig`) using the parameters specified by
    the multiscale object (`Mobj`) and the following default parameters:
    Scales = 3, Butterworth LPF Order = 6, Butterworth LPF cutoff frequency
    at scale (T): Fc = 0.5/T. 
    If the entropy function specified by `Mobj` is SampEn or ApEn, rMSEn 
    updates the threshold radius of the data sequence (Xt) at each scale
    to 0.2std(Xt) if no `r` value is provided by Mobj, or r.std(Xt) if `r`
    is specified.

        MSx, CI = rMSEn(Sig::AbstractArray{T,1} where T<:Real, Mobj::NamedTuple; Scales::Int=3, 
                            F_Order::Int=6, F_Num::Float64=0.5, RadNew::Int=0, Plotx::Bool=false)

    Returns a vector of refined multiscale entropy values (`MSx`) and the complexity 
    index (`CI`) of the data sequence (`Sig`) using the parameters specified by
    the multiscale object (`Mobj`) and the following 'keyword' arguments:

    # Arguments:
    `Scales`   - Number of temporal scales, an integer > 1 (default = 3) \n
    `F_Order`  - Butterworth low-pass filter order, a positive integer
                 (default: 6) \n
    `F_Num`    - Numerator of Butterworth low-pass filter cutoff frequency,
                 a scalar value in range [0 < `F_Num` < 1]. The cutoff frequency
                 at each scale (T) becomes: Fc = `F_Num/T.  (default: 0.5) \n
    `RadNew`   - Radius rescaling method, an integer in the range [1 4].
                 When the entropy specified by `Mobj` is `SampEn` or `ApEn`, 
                 `RadNew` allows the radius threshold to be updated at each 
                 time scale (Xt). If a radius value is specified by `Mobj` (`r`),
                 this becomes the rescaling coefficient, otherwise it is set
                 to 0.2 (default). The value of `RadNew` specifies one of the 
                 following methods:\n
                 [1]    Standard Deviation          - r*std(Xt)\n
                 [2]    Variance                    - r*var(Xt) \n
                 [3]    Mean Absolute Deviation     - r*mean_ad(Xt) \n
                 [4]    Median Absolute Deviation   - r*med_ad(Xt)\n
    `Plotx`    - When `Plotx` == true, returns a plot of the entropy value at each
                 time scale (i.e. the multiscale entropy curve)  [default: false] \n

    # See also `MSobject`, `MSEn`, `cMSEn`, `hMSEn`, `SampEn`, `ApEn`, `XMSEn`

    # References:
        [1] Madalena Costa, Ary Goldberger, and C-K. Peng,
            "Multiscale entropy analysis of complex physiologic time series."
            Physical review letters
            89.6 (2002): 068102.

        [2] Vadim V. Nikulin, and Tom Brismar,
            "Comment on “Multiscale entropy analysis of complex physiologic
            time series”." 
            Physical review letters 
            92.8 (2004): 089803.

        [3] Madalena Costa, Ary L. Goldberger, and C-K. Peng. 
            "Costa, Goldberger, and Peng reply." 
            Physical Review Letters
            92.8 (2004): 089804.

        [4] José Fernando Valencia, et al.,
            "Refined multiscale entropy: Application to 24-h holter 
            recordings of heart period variability in healthy and aortic 
            stenosis subjects." 
            IEEE Transactions on Biomedical Engineering 
            56.9 (2009): 2202-2213.

        [5] Puneeta Marwaha and Ramesh Kumar Sunkaria,
            "Optimal selection of threshold value ‘r’for refined multiscale
            entropy." 
            Cardiovascular engineering and technology 
            6.4 (2015): 557-576.

    """
    function rMSEn(Sig::AbstractArray{T,1} where T<:Real, Mobj::NamedTuple; Scales::Int=3, 
        F_Order::Int=6, F_Num::Float64=0.5, RadNew::Int=0, Plotx::Bool=false)

    (size(Sig,1)>10) ? nothing : error("Sig:   must be a numeric vector" )
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
         with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (F_Order>1) ? nothing :  error("F_Order:     must be an integer > 1")
    (0 < F_Num < 1) ? nothing : error("F_Num:     must be a scalar in range 0 < F_Num < 1")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("SampEn","ApEn"))) ? nothing :
    error("RadNew:  must be 0, or an integer in range [1 4] with entropy function `SampEn` or `ApEn`")
        
    MSx = zeros(Scales)
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
    (RadNew==0 && String(Symbol(Mobj.Func)) in ("SampEn","ApEn")) ? RadNew=1 : nothing
    if RadNew > 0
        if RadNew == 1
            Rnew = x -> std(x, corrected=false)
        elseif RadNew == 2
            Rnew = x -> var(x, corrected=false)
        elseif RadNew == 3
            Rnew = x -> mean(abs.(x .- mean(x)))
        elseif RadNew == 4
            Rnew = x -> median(abs.(x .- median(x)))    
        end

        if haskey(Mobj,:r)
            Cx = Mobj.r
        else
            Cy = ("Standard Deviation","Variance","Mean Abs Deviation", "Median Abs Deviation")
            @warn("No radius value provided in Mobj.
                Default set to 0.2*$(Cy[RadNew]) of each new time-series.")            
            Cx = .2
        end
    end

    for T = 1:Scales
        print(". ")
        Temp = refined(Sig,T,F_Order,F_Num)
        RadNew > 0 ? Args = (Args..., r=Cx*Rnew(Temp[:])) : nothing      
        Tempx = Mobj.Func(Temp[:]; Args...)
        #MSx[T] = Tempx[1][end]
        typeof(Tempx)<:Tuple ? MSx[T] = Tempx[1][end] : MSx[T] = Tempx[end]
    end
    CI = sum(MSx)
    print("\n")
    if any(isnan.(MSx))
        println("Some entropy values may be undefined.")
    end

    if Plotx
        p1 = plot(1:Scales, MSx, c=RGB(8/255, 63/255, 77/255), lw=3)
        scatter!(1:Scales, MSx, markersize=6, c=RGB(1, 0, 1),
        xlabel = "Scale Factor", ylabel = "Entropy Value", 
        guidefont = font(12, "arial", RGB(7/255, 54/255, 66/255)),
        tickfontsize = 10, tickfontfamily="arial", legend=false,
        title = "Refined Multiscale $(Mobj.Func)",
        plot_titlefontsize=16, plot_titlefontcolor=RGB(7/255, 54/255, 66/255))
        display(p1)
    end

    return MSx, CI
    end

    function refined(Z,sx,P1,P2)
        Yt = filtfilt(digitalfilter(Lowpass(P2/sx), Butterworth(P1)), Z)
        return Yt[1:sx:end]
    end

end

"""


Copyright 2021 Matthew W. Flood, EntropyHub

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

For Terms of Use see https://github.com/MattWillFlood/EntropyHub

"""