module _cMSEn
export cMSEn
using Statistics: std, mean, median, var
using Plots
    """
    # MSx, CI = cMSEn(`Sig`, `Mobj`) 

    Returns a vector of composite multiscale entropy values (`MSx`) for the data 
    sequence (`Sig`) using the parameters specified by the multiscale object 
    (`Mobj`) using the composite multiscale entropy method over 3 temporal scales.

    # MSx, CI = cMSEn(`Sig`, `Mobj`, 'keyword' = value, ...)

    Returns a vector of composite multiscale entropy values (`MSx`) of the 
    data sequence (`Sig`) using the parameters specified by the multiscale 
    object (`Mobj`) and the following 'keyword' arguments:

    # Arguments:
    `Scales`   - Number of temporal scales, an integer > 1   (default: 3) \n
    `RadNew`   - Radius rescaling method, an integer in the range [1 4].
                 When the entropy specified by `Mobj` is `SampEn` or `ApEn`, 
                 RadNew allows the radius threshold to be updated at each 
                 time scale (Xt). If a radius value is specified by `Mobj` (`r`),
                 this becomes the rescaling coefficient, otherwise it is set
                 to 0.2 (default). The value of RadNew specifies one of the 
                 following methods:\n
                 [1]    Standard Deviation          - r*std(Xt)\n
                 [2]    Variance                    - r*var(Xt) \n
                 [3]    Mean Absolute Deviation     - r*mean_ad(Xt) \n
                 [4]    Median Absolute Deviation   - r*med_ad(Xt)\n
    `Refined`  - Refined-composite MSEn method. When `Refined` == true and the 
                 entropy function specified by Mobj is `SampEn`, cMSEn returns 
                 the refined-composite multiscale entropy (rcMSEn) [default: false]\n
    `Plotx`    - When `Plotx` == true, returns a plot of the entropy value at each
                time scale (i.e. the multiscale entropy curve) [default: false]\n

    # See also `MSobject`, `rMSEn`, `MSEn`, `hMSEn`, `SampEn`, `ApEn`, `XMSEn`

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

     [4] Shuen-De Wu, et al.,
          "Time series analysis using composite multiscale entropy." 
          Entropy 
          15.3 (2013): 1069-1084.

     [5] Shuen-De Wu, et al.,
          "Analysis of complex time series using refined composite 
          multiscale entropy." 
          Physics Letters A 
          378.20 (2014): 1369-1374.

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
    function cMSEn(Sig::AbstractArray{T,1} where T<:Real, Mobj::NamedTuple;  
            Scales::Int=3, RadNew::Int=0, Refined::Bool=false, Plotx::Bool=false)

    (size(Sig,1)>10) ? nothing : error("Sig:   must be a numeric vector" )
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
                with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("SampEn","ApEn"))) ? nothing :
        error("RadNew:     must be 0, or an integer in range [1 4] with 
                entropy function 'SampEn' or 'ApEn'")
                       
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
        N = Int(T*floor(length(Temp)/T))
        Temp2 = zeros(T)
        Temp3 = zeros(T)

        for k = 1:T
            print(". ")
            RadNew > 0 ? Args = (Args..., r=Cx*Rnew(Temp[k:T:N])) : nothing

            if Refined
                _, Ma, Mb = Mobj.Func(Temp[k:T:N]; Args...)
                Temp2[k] = Ma[end]
                Temp3[k] = Mb[end]
            else
                Temp2 = Mobj.Func(Temp[k:T:N]; Args...)
                typeof(Temp2)<:Tuple ? Temp3[k] = Temp2[1][end] : Temp3[k] = Temp2[end]
                #Temp3[k] = Temp2[1][end]
            end
        end    
        
        Refined ? MSx[T] = -log(sum(Temp2)/sum(Temp3)) : MSx[T] = mean(Temp3)
        
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
        plot_titlefontsize=16, plot_titlefontcolor=RGB(7/255, 54/255, 66/255)) #ylim=(0,maximum(MSx)+.2),
        display(p1)
    end

    return MSx, CI
    end

    function modified(Z,sx)
        Ns = length(Z) - sx + 1
        Y = zeros(Ns)
        for k = 1:Ns
            Y[k] = mean(Z[k:k+sx-1])
        end
        return Y 
    end

end