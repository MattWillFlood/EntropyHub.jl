module _cXMSEn
export cXMSEn
using Statistics: std, mean, median, var
using Plots
using DSP: conv
    """
        MSx, CI = cXMSEn(Sig1, Sig2, Mobj) 

    Returns a vector of composite multiscale cross-entropy values (`MSx`) 
    between two univariate data sequences contained in `Sig1` and `Sig2` using the 
    parameters specified by the multiscale object (`Mobj`) using the composite 
    multiscale method (cMSE) over 3 temporal scales.

        MSx, CI = cXMSEn(Sig1::AbstractVector{T} where T<:Real, Sig2::AbstractVector{T} where T<:Real, Mobj::NamedTuple; 
                              Scales::Int=3, RadNew::Int=0, Refined::Bool=false, Plotx::Bool=false)

    Returns a vector of composite multiscale cross-entropy values (`MSx`) 
    between the data sequences contained in `Sig1` and `Sig2` using the parameters
    specified by the multiscale object (`Mobj`) and the following keyword arguments:

    # Arguments:
    `Scales`   - Number of temporal scales, an integer > 1   (default: 3)\n
    `RadNew`   - Radius rescaling method, an integer in the range [1 4].
                 When the entropy specified by Mobj is `XSampEn` or `XApEn`, 
                 RadNew rescales the radius threshold of the sub-sequences
                 at each time scale (Ykj). If a radius value is specified by
                 `Mobj` (`r`), this becomes the rescaling coefficient, otherwise
                 it is set to 0.2 (default). The value of RadNew specifies
                 one of the following methods:\n
                 [1]    Pooled Standard Deviation          - r*std(Ykj)\n
                 [2]    Pooled Variance                    - r*var(Ykj)\n
                 [3]    Total Mean Absolute Deviation      - r*mean_ad(Ykj)\n
                 [4]    Total Median Absolute Deviation    - r*med_ad(Ykj,1)\n
    `Refined`  - Refined-composite XMSEn method. When `Refined` == true and the 
                 entropy function specified by `Mobj` is `XSampEn` or `XFuzzEn`, cXMSEn
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


    """
    function cXMSEn(Sig1::AbstractVector{T} where T<:Real, Sig2::AbstractVector{T} where T<:Real, Mobj::NamedTuple; 
                        Scales::Int=3, RadNew::Int=0, Refined::Bool=false, Plotx::Bool=false)

    (size(Sig1,1)>=10)  && (size(Sig2,1)>=10) ? nothing : error("Sig1/Sig2:   sequences must have >= 10 values")
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("XSampEn","XApEn"))) ? 
            nothing : error("RadNew:     must be 0, or an integer in range [1 4] with entropy function 'XSampEn' or 'XApEn'")
    (Refined && String(Symbol(Mobj.Func)) in ("XSampEn","XFuzzEn")) || Refined==false ? 
            nothing : error("Refined:     If Refined == true, the base entropy function must be 'XSampEn' or 'XFuzzEn'")
   
    String(Symbol(Mobj.Func))=="XSampEn" ? Mobj = merge(Mobj,(Vcp=false,)) : nothing

    if Refined && String(Symbol(Mobj.Func))=="XFuzzEn"
        Tx = 1;
        "Logx" in String.(Symbol.(keys(Mobj))) ? Logx = Mobj.Logx : Logx = exp(1);

    elseif Refined && String(Symbol(Mobj.Func))=="XSampEn"
        Tx = 0;
        "Logx" in String.(Symbol.(keys(Mobj))) ? Logx = Mobj.Logx : Logx = exp(1);
    else
        Tx = 0;
    end
    
    MSx = zeros(Scales)
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
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
        TempA, TempB = modified(Sig1, Sig2, T, Tx) 
        N1 = Int(T*floor(size(TempA,1)/T))
        N2 = Int(T*floor(size(TempB,1)/T)) 
        Temp3 = zeros(T)
        Temp2 = zeros(T)

        for k = 1:T
            print(" .")
            RadNew > 0 ? Args = (Args..., r=Cx*Rnew(TempA[k:T:N1], TempB[k:T:N2])) : nothing

            if Refined == 1  
                _, Ma, Mb = Mobj.Func(TempA[k:T:N1], TempB[k:T:N2]; Args...)
                Temp2[k] = Ma[end]
                Temp3[k] = Mb[end]
            else
                Temp2 = Mobj.Func(TempA[k:T:N1], TempB[k:T:N2]; Args...)
                typeof(Temp2)<:Tuple ? Temp3[k] = Temp2[1][end] : Temp3[k] = Temp2[end]
                #Temp3[k] = Temp2[1][end]
            end
        end    
        
        # Refined == 1 ? MSx[T] = -log(sum(Temp2)/sum(Temp3)) : MSx[T] = mean(Temp3)
   
        if Refined && Tx==0
            MSx[T] = -log(sum(Temp2)/sum(Temp3))/log(Logx)
        elseif Refined && Tx==1
            MSx[T] = -log(sum(Temp3)/sum(Temp2))/log(Logx)
        else
            MSx[T] = mean(Temp3)
        end   
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

    function modified(Za, Zb, sx, Tx)
        if Tx==1
            Y1 = [std(Za[x:x+sx-1], corrected=false) for x in 1:(length(Za)-sx+1)]
            Y2 = [std(Zb[x:x+sx-1], corrected=false) for x in 1:(length(Zb)-sx+1)]
        else
            Y1 = (conv(Za,ones(Int, sx))/sx)[sx:end-sx+1][:]
            Y2 = (conv(Zb,ones(Int, sx))/sx)[sx:end-sx+1][:]
        end
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