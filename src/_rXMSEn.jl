module _rXMSEn
export rXMSEn
using DSP.Filters: filtfilt, Butterworth, Lowpass, digitalfilter
using Statistics: std, mean, median, var
using Plots
    """
        MSx, CI = rXMSEn(Sig1, Sig2, Mobj)

    Returns a vector of refined multiscale cross-entropy values (`MSx`) and
    the complexity index (`CI`) between the data sequences contained in `Sig1` and `Sig2`
    using the parameters specified by the multiscale object (`Mobj`) and the
    following default parameters:   Scales = 3, Butterworth LPF Order = 6,
    Butterworth LPF cutoff frequency at scale (T): Fc = 0.5/T. 
    If the entropy function specified by `Mobj` is `XSampEn` or `XApEn`, 
    `rMSEn` updates the threshold radius of the data sequences (Xt) at each scale
    to 0.2*SDpooled(Xa, Xb) when no `r` value is provided by `Mobj`, or 
    `r`*SDpooled(Xa, Xb) if `r` is specified.
     
        MSx, CI = rXMSEn(Sig1::AbstractVector{T} where T<:Real, Sig2::AbstractVector{T} where T<:Real, Mobj::NamedTuple;
                         Scales::Int=3, F_Order::Int=6, F_Num::Float64=0.5, RadNew::Int=0, Plotx::Bool=false)

    Returns a vector of refined multiscale cross-entropy values (`MSx`) and 
    the complexity index (`CI`) between the data sequences contained in `Sig1` and `Sig2`
    using the parameters specified by the multiscale object (`Mobj`) and the
    following keyword arguments:

    # Arguments:
    `Scales`   - Number of temporal scales, an integer > 1 (default: 3)  \n
    `F_Order`  - Butterworth low-pass filter order, a positive integer (default: 6)  \n
    `F_Num`    - Numerator of Butterworth low-pass filter cutoff frequency,
                 a scalar value in range [0 < `F_Num` < 1]. The cutoff frequency
                 at each scale (T) becomes: Fc = `F_Num`/T.  (default: 0.5) \n
    `RadNew`   - Radius rescaling method, an integer in the range [1 4].
                 When the entropy specified by `Mobj` is `XSampEn` or `XApEn`, 
                 `RadNew` allows the radius threshold to be updated at each 
                 time scale (Xt). If a radius value is specified by `Mobj` (`r`),
                 this becomes the rescaling coefficient, otherwise it is set
                 to 0.2 (default). The value of `RadNew` specifies one of the 
                 following methods: \n
                 [1]    Pooled Standard Deviation          - r*std(Xt) \n
                 [2]    Pooled Variance                    - r*var(Xt) \n
                 [3]    Total Mean Absolute Deviation      - r*mean_ad(Xt) \n
                 [4]    Total Median Absolute Deviation    - r*med_ad(Xt) \n
    `Plotx`    - When `Plotx` == true, returns a plot of the entropy value at 
                 each time scale (i.e. the multiscale entropy curve)
                 [default = false] \n
  
    # See also `MSobject`, `XMSEn`, `cXMSEn`, `hXMSEn`, `XSampEn`, `XApEn`, `MSEn`
  
    # References:
        [1] Matthew W. Flood (2021), 
            "EntropyHub - An open source toolkit for entropic time series analysis"
            PLoS ONE 16(11):e0295448, 
            DOI:  10.1371/journal.pone.0259448
            https://www.EntropyHub.xyz
  
        [2]   Rui Yan, Zhuo Yang, and Tao Zhang,
            "Multiscale cross entropy: a novel algorithm for analyzing two
            time series." 
            5th International Conference on Natural Computation. 
            Vol. 1, pp: 411-413 IEEE, 2009.
  
        [3] José Fernando Valencia, et al.,
            "Refined multiscale entropy: Application to 24-h holter 
            recordings of heart period variability in healthy and aortic 
            stenosis subjects." 
            IEEE Transactions on Biomedical Engineering 
            56.9 (2009): 2202-2213.
  
        [4] Puneeta Marwaha and Ramesh Kumar Sunkaria,
            "Optimal selection of threshold value ‘r’for refined multiscale
            entropy." 
            Cardiovascular engineering and technology 
            6.4 (2015): 557-576.
  
        [5] Yi Yin, Pengjian Shang, and Guochen Feng, 
            "Modified multiscale cross-sample entropy for complex time 
            series."
            Applied Mathematics and Computation 
            289 (2016): 98-110.
  
        [6] Antoine Jamin, et al,
            "A novel multiscale cross-entropy method applied to navigation 
            data acquired with a bike simulator." 
            41st annual international conference of the IEEE EMBC
            IEEE, 2019.
  
        [7] Antoine Jamin and Anne Humeau-Heurtier. 
            "(Multiscale) Cross-Entropy Methods: A Review." 
            Entropy 
            22.1 (2020): 45.
   
    """
    function rXMSEn(Sig1::AbstractVector{T} where T<:Real, Sig2::AbstractVector{T} where T<:Real, Mobj::NamedTuple;
        Scales::Int=3, F_Order::Int=6, F_Num::Float64=0.5, RadNew::Int=0, Plotx::Bool=false)
        
    (size(Sig1,1)>=10)  && (size(Sig2,1)>=10) ? nothing : error("Sig1/Sig2:   sequences must have >= 10 values")
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
         with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (F_Order>1) ? nothing :  error("F_Order:     must be an integer > 1")
    (0 < F_Num < 1) ? nothing : error("F_Num:     must be a scalar in range 0 < F_Num < 1")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("XSampEn","XApEn"))) ? nothing :
    error("RadNew:  must be 0, or an integer in range [1 4] with entropy function `XSampEn` or `XApEn`")
        
    String(Symbol(Mobj.Func))=="XSampEn" ? Mobj = merge(Mobj,(Vcp=false,)) : nothing

    MSx = zeros(Scales)
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
    (RadNew==0 && String(Symbol(Mobj.Func)) in ("XSampEn","XApEn")) ? RadNew=1 : nothing
    if RadNew > 0
        if RadNew == 1
            Rnew = (x,y) -> sqrt((var(x)*(size(x,1)-1) + var(y)*(size(y,1)-1))/(size(x,1)+size(y,1)-1))
        elseif RadNew == 2
            Rnew = (x,y) -> ((var(x)*(size(x,1)-1) + var(y)*(size(y,1)-1))/(size(x,1)+size(y,1)-1))
        elseif RadNew == 3
            Rnew = (x,y) -> mean(abs.(vcat(x,y) .- mean(vcat(x,y))))
        elseif RadNew == 4
            Rnew = (x,y) -> median(abs.(vcat(x,y) .- median(vcat(x,y))))    
        end

        if haskey(Mobj,:r)
            Cx = Mobj.r
        else
            Cy = ("Pooled Standard Deviation","Pooled Variance","Total Mean Abs Deviation", "Total Median Abs Deviation")
            @warn("No radius value provided in Mobj.
                Default set to 0.2*$(Cy[RadNew]) of each new time-series.")            
            Cx = .2
        end
    end

    for T = 1:Scales
        print(" .")
        TempA, TempB = refined(Sig1, Sig2,T,F_Order,F_Num)
        RadNew > 0 ? Args = (Args..., r=Cx*Rnew(TempA, TempB)) : nothing      
        Tempx = Mobj.Func(TempA, TempB; Args...)
        typeof(Tempx)<:Tuple ? MSx[T] = Tempx[1][end] : MSx[T] = Tempx[end]
        #MSx[T] = Tempx[1][end]
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

    function refined(Za, Zb,sx,P1,P2)
        Y1 = filtfilt(digitalfilter(Lowpass(P2/sx), Butterworth(P1)), Za)[:]
        Y2 = filtfilt(digitalfilter(Lowpass(P2/sx), Butterworth(P1)), Zb)[:]        
        return Y1, Y2
    end
end

"""
Copyright 2024 Matthew W. Flood, EntropyHub
  
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