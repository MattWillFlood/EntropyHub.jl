module _rXMSEn
export rXMSEn
using DSP.Filters: filtfilt, Butterworth, Lowpass, digitalfilter
using Statistics: std, mean, median, var
using Plots
    """
    # MSx, CI = rXMSEn(`Sig`, `Mobj`)

    Returns a vector of refined multiscale cross-entropy values (`MSx`) and
    the complexity index (`CI`) between the data sequences contained in `Sig`
    using the parameters specified by the multiscale object (`Mobj`) and the
    following default parameters:   Scales = 3, Butterworth LPF Order = 6,
    Butterworth LPF cutoff frequency at scale (T): Fc = 0.5/T. 
    If the entropy function specified by `Mobj` is `XSampEn` or `XApEn`, 
    `rMSEn` updates the threshold radius of the data sequences (Xt) at each scale
    to 0.2*std(Xt) when no `r` value is provided by `Mobj`, or r*std(Xt) if 
    `r` is specified.
     
    # MSx, CI = rXMSEn(`Sig`, `Mobj`, `keyword` = value, ...)
    Returns a vector of refined multiscale cross-entropy values (`MSx`) and 
    the complexity index (`CI`) between the data sequences contained in `Sig`
    using the parameters specified by the multiscale object (`Mobj`) and the
    following keyword arguments:
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
                 following methods:  \n
                 [1]    Standard Deviation          - r*std(Xt)  (default)  \n
                 [2]    Variance                    - r*var(Xt)  \n
                 [3]    Mean Absolute Deviation     - r*mean_ad(Xt)  \n
                 [4]    Median Absolute Deviation   - r*med_ad(Xt,1) \n
    `Plotx`    - When `Plotx` == true, returns a plot of the entropy value at 
                 each time scale (i.e. the multiscale entropy curve)
                 [default = false] \n
  
    # See also `MSobject`, `XMSEn`, `cXMSEn`, `hXMSEn`, `XSampEn`, `XApEn`, `MSEn`
  
    # References:
      [1]   Matthew W. Flood,
            "rXMSEn - EntropyHub Project"
            2021, https://github.com/MattWillFlood/EntropyHub
  
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
    function rXMSEn(Sig::AbstractArray{T,2} where T<:Real, Mobj::NamedTuple; Scales::Int=3, 
        F_Order::Int=6, F_Num::Float64=0.5, RadNew::Int=0, Plotx::Bool=false)
        
    size(Sig,1) == 2  ?  Sig = Sig' : nothing

    (size(Sig,1)>10) ? nothing : error("Sig:   must be a numeric vector" )
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
         with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (F_Order>1) ? nothing :  error("F_Order:     must be an integer > 1")
    (0 < F_Num < 1) ? nothing : error("F_Num:     must be a scalar in range 0 < F_Num < 1")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("XSampEn","XApEn"))) ? nothing :
    error("RadNew:  must be 0, or an integer in range [1 4] with entropy function `XSampEn` or `XApEn`")
        
    MSx = zeros(Scales)
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
    (RadNew==0 && String(Symbol(Mobj.Func)) in ("XSampEn","XApEn")) ? RadNew=1 : nothing
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
            Cy = ("Standard Deviation","Variance","Mean Abs Deviation",
                    "Median Abs Deviation")
            @warn("No radius value provided in Mobj.
                Default set to 0.2*$(Cy[RadNew]) of each new time-series.")            
            Cx = .2
        end
    end

    for T = 1:Scales
        print(" .")
        Temp = refined(Sig,T,F_Order,F_Num)
        RadNew > 0 ? Args = (Args..., r=Cx*Rnew(Temp[:])) : nothing      
        Tempx = Mobj.Func(Temp; Args...)
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

    function refined(Z,sx,P1,P2)
        Yt = filtfilt(digitalfilter(Lowpass(P2/sx), Butterworth(P1)), Z)
        return Yt[1:sx:end,:]
    end
end