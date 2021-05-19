module _cXMSEn
export cXMSEn
using Statistics: std, mean, median, var
using Plots
    """
    # MSx, CI = cXMSEn(`Sig`, `Mobj`) 
    Returns a vector of composite multiscale cross-entropy values (`MSx`) 
    between two univariate data sequences contained in `Sig` using the 
    parameters specified by the multiscale object (`Mobj`) using the composite 
    multiscale method (cMSE) over 3 temporal scales.

    # MSx, CI = cXMSEn(`Sig`, `Mobj`, 'keyword' = value, ...)
    Returns a vector of composite multiscale cross-entropy values (`MSx`) 
    between the data sequences contained in `Sig` using the parameters
    specified by the multiscale object (`Mobj`) and the following keyword arguments:
    `Scales`   - Number of temporal scales, an integer > 1   (default: 3)
    `RadNew`   - Radius rescaling method, an integer in the range [1 4].
                 When the entropy specified by Mobj is `XSampEn` or `XApEn`, 
                 RadNew rescales the radius threshold of the sub-sequences
                 at each time scale (Ykj). If a radius value is specified by
                 `Mobj` (`r`), this becomes the rescaling coefficient, otherwise
                 it is set to 0.2 (default). The value of RadNew specifies
                 one of the following methods:
                 [1]    Standard Deviation          - r*std(Ykj)
                 [2]    Variance                    - r*var(Ykj)
                 [3]    Mean Absolute Deviation     - r*mean_ad(Ykj)
                 [4]    Median Absolute Deviation   - r*med_ad(Ykj,1)
    `Refined`  - Refined-composite XMSEn method. When `Refined` == true and the 
                 entropy function specified by `Mobj` is 'XSampEn', cXMSEn
                 returns the refined-composite multiscale entropy (rcXMSEn).
                 (default: false)
    `Plotx`    - When `Plotx` == true, returns a plot of the entropy value at each 
                 time scale (i.e. the multiscale entropy curve) [default: false]

    # See also `MSobject`, `XMSEn`, `rXMSEn`, `hXMSEn`, `XSampEn`, `XApEn`, `cMSEn`

    # References:
       [1] Rui Yan, Zhuo Yang, and Tao Zhang,
            "Multiscale cross entropy: a novel algorithm for analyzing two
            time series." 
            5th International Conference on Natural Computation. 
            Vol. 1, pp: 411-413 IEEE, 2009.

       [2] Yi Yin, Pengjian Shang, and Guochen Feng, 
            "Modified multiscale cross-sample entropy for complex time 
            series."
            Applied Mathematics and Computation 
            289 (2016): 98-110.

       [3] Madalena Costa, Ary Goldberger, and C-K. Peng,
            "Multiscale entropy analysis of complex physiologic time series."
            Physical review letters
            89.6 (2002): 068102.

       [4] Antoine Jamin, et al,
            "A novel multiscale cross-entropy method applied to navigation 
            data acquired with a bike simulator." 
            41st annual international conference of the IEEE EMBC
            IEEE, 2019.

       [5] Antoine Jamin and Anne Humeau-Heurtier. 
            "(Multiscale) Cross-Entropy Methods: A Review." 
            Entropy 
            22.1 (2020): 45.

       [6] Shuen-De Wu, et al.,
            "Time series analysis using composite multiscale entropy." 
            Entropy 
            15.3 (2013): 1069-1084.

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
    function cXMSEn(Sig::AbstractArray{T,2} where T<:Real, Mobj::NamedTuple; 
        Scales::Int=3, RadNew::Int=0, Refined::Bool=false, Plotx::Bool=false)

    size(Sig,1) == 2 ? Sig = Sig' : nothing

    (size(Sig,1)>10) ? nothing : error("Sig:   must be a numeric vector" )
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy 
            object created with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("XSampEn","XApEn"))) ? 
            nothing : error("RadNew:     must be 0, or an integer in range [1 4] with 
                entropy function 'XSampEn' or 'XApEn'")
   
    MSx = zeros(Scales)
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
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
        Temp = modified(Sig,T) 
        N = Int(T*floor(size(Temp,1)/T))
        Temp3 = zeros(T)
        Temp2 = zeros(T)

        for k = 1:T
            print(" .")
            RadNew > 0 ? Args = (Args..., r=Cx*Rnew(Temp[k:T:N,:])) : nothing

            if Refined == 1  
                _, Ma, Mb = Mobj.Func(Temp[k:T:N,:]; Args...)
                Temp2[k] = Ma[end]
                Temp3[k] = Mb[end]
            else
                Temp2 = Mobj.Func(Temp[k:T:N,:]; Args...)
                typeof(Temp2)<:Tuple ? Temp3[k] = Temp2[1][end] : Temp3[k] = Temp2[end]
                #Temp3[k] = Temp2[1][end]
            end
        end    
        
        Refined == 1 ? MSx[T] = -log(sum(Temp2)/sum(Temp3)) : MSx[T] = mean(Temp3)
    end
    CI = sum(MSx)
    print("\n")
    if any(isnan.(MSx))
        println("Some entropy values may be undefined.")
    end

    if Plotx
        Refined ? strx = "Refined-Composite" : strx = "Composite"

        p1 = plot(1:Scales, MSx, c=RGB(8/255, 63/255, 77/255), lw=3)
        scatter!(1:Scales, MSx, markersize=6, c=RGB(1, 0, 1),
        xlabel = "Scale Factor", ylabel = "Entropy Value", 
        guidefont = font(12, "arial", RGB(7/255, 54/255, 66/255)),
        tickfontsize = 10, tickfontfamily="arial", legend=false,
        title = strx*" Multiscale $(Mobj.Func)",
        plot_titlefontsize=16, plot_titlefontcolor=RGB(7/255, 54/255, 66/255)) 
        display(p1)
    end

    return MSx, CI
    end

    function modified(Z,sx)
        Ns = size(Z,1) - sx + 1
        Y = zeros(Ns,2)
        for k = 1:Ns
            Y[k,1] = mean(Z[k:k+sx-1,1])
            Y[k,2] = mean(Z[k:k+sx-1,2])
        end
        return Y 
    end

end