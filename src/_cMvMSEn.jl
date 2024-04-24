module _cMvMSEn
export cMvMSEn
using Statistics: std, mean, median, var
using Plots
using DSP: conv
    """
        MSx, CI = cMvMSEn(Data, Mobj) 
    
     Returns a vector of composite multivariate multiscale entropy values (`MSx`) and the complexity 
     index (`CI`) of the data sequences in `Data` using the parameters specified 
     by the multiscale object (`Mobj`) over 3 temporal scales with coarse-graining (default). 
    
     ________________________________________________________________________
     
    !!! note

        By default, the `MvSampEn` and `MvFuzzEn` multivariate entropy algorithms
        estimate entropy values using the "full"  method by comparing delay vectors 
        across all possible `m+1` expansions of the embedding space as applied in [1].
        These methods are not lower-bounded to 0, like most entropy algorithms,
        so `MvMSEn` may return negative entropy values if the base multivariate 
        entropy function is `MvSampEn` and `MvFuzzEn`, even for stochastic processes...

     ________________________________________________________________________
           
        MSx, CI = cMvMSEn(Data, Mobj, Refined = True) 
        
     Returns a vector of refined-composite multiscale entropy values (`MSx`) for the data 
     sequences in (`Data`) using the parameters specified by the multiscale object 
     (`Mobj`) using the refined-composite multivariate multiscale entropy method (rcMSE) over 3 temporal
     scales. When `Refined == true`, the base entropy method must be `MvSampEn` or `MvFuzzEn`.
     If the entropy method is `MvSampEn`, cMvMSEn employs the method described in [1]. 
     If the entropy method is `MvFuzzEn`, cMvMSEn employs the method described in [5]. 

        MSx, CI = cMVMSEn(Data::AbstractArray{T,2} where T<:Real, Mobj::NamedTuple; Scales::Int=3, Refined::Bool="false", Plotx::Bool=false)

     Returns a vector of multivariate multiscale entropy values (`MSx`) and the complexity 
     index (`CI`) of the data sequences in `Data` using the parameters specified by
     the multiscale object (`Mobj`) and the following keyword arguments:

     # Arguments:
     `Scales`   - Number of temporal scales, an integer > 1   (default: 3) \n
     `Refined`  - Refined-composite MvMSEn method. When `Refined == True` and the 
                  entropy function specified by `Mobj` is `MvSampEn` or `MvFuzzEn`, 
                  `cMvMSEn` returns the refined-composite multivariate multiscale entropy (rcMSEn) [default: False]
     `Plotx`    - When Plotx == true, returns a plot of the entropy value at each
                  time scale (i.e. the multiscale entropy curve) [default: false]\n

     # See also  `MvMSEn`, `MSobject`, `MvFuzzEn`, `MvSampEn`, `MvPermEn`, `MvCoSiEn`,  `MvDispEn`
    
     # References:
        [1] Shuen-De Wu, et al.,
            "Time series analysis using composite multiscale entropy."
            Entropy
            15.3 (2013): 1069-1084.

        [2] Shuen-De Wu, et al.,
            "Analysis of complex time series using refined composite
            multiscale entropy."
            Physics Letters A
            378.20 (2014): 1369-1374.

        [3] Ahmed Mosabber Uddin, Danilo P. Mandic
            "Multivariate multiscale entropy: A tool for complexity
            analysis of multichannel data."
            Physical Review E 84.6 (2011): 061918.

        [4] Ahmed Mosabber Uddin, Danilo P. Mandic
            "Multivariate multiscale entropy analysis."
            IEEE signal processing letters 19.2 (2011): 91-94.

        [5] Azami, Alberto Fernández, Javier Escudero.
            "Refined multiscale fuzzy entropy based on standard deviation for
            biomedical signal analysis."
            Medical & biological engineering & computing 55 (2017): 2037-2052.

    """
    function cMvMSEn(Data::AbstractArray{T,2} where T<:Real, Mobj::NamedTuple; Scales::Int=3, Refined::Bool=false, Plotx::Bool=false)

    N, Dn = size(Data)   
    (N>10) && (Dn>1) ? nothing : error("Data:   must be an NxM matrix where N>10 and M>1")
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
                with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (Refined && String(Symbol(Mobj.Func)) in ("MvSampEn","MvFuzzEn")) || Refined==false ? 
    nothing : error("Refined:     If Refined == true, the base entropy function must be 'MvSampEn' or 'MvFuzzEn'")
    String(Symbol(Mobj.Func))[1:2] != "Mv" ? error("Base entropy estimator must be a multivariate entropy method. ",
    "To perform univariate multiscale entropy estimation, use MSEn().") : nothing    

    if Refined
        if  string(Mobj.Func)== "MvFuzzEn"
            Tx = 1
        elseif string(Mobj.Func) == "MvSampEn"
            Tx = 0
        end        
        "Logx" in String.(Symbol.(keys(Mobj))) ? Logx = Mobj.Logx : Logx = exp(1)        
    else
        Tx = 0
    end
  
    MSx = zeros(Scales)
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
    for T = 1:Scales
        print(". ")
        Temp = modified(Data,T,Tx,Dn)       
        N = Int(T*floor(size(Temp,1)/T))
        Ma = zeros(T)
        Mb = zeros(T)
        for k = 1:T
            print(". ")

            if Refined
                _, Ma[k], Mb[k], _ = Mobj.Func(Temp[k:T:N,:]; Args...)
            else
                Ma[k], _, _, _ = Mobj.Func(Temp[k:T:N,:]; Args...)
            end
        end
        Refined ? MSx[T] = -log(Logx, sum(Mb)/sum(Ma)) : MSx[T] = mean(Ma)  
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
        title = "$(strx) Multivariate Multiscale $(string(Mobj.Func)[3:end])",
        plot_titlefontsize=16, plot_titlefontcolor=RGB(7/255, 54/255, 66/255)) #ylim=(0,maximum(MSx)+.2),
        display(p1)
    end

    return MSx, CI
    end


    function modified(Z, sx, Tx, Dn)
        if Tx==0
            Y = (conv(Z,ones(Int,sx))/sx)[sx:end-sx+1,:]
        else
            Ns = size(Z,1)-sx+1
            Y = zeros(Ns,Dn)
            for k in 1:Dn
                Y[:,k] = map(x -> std(Z[x:x+sx-1,k], corrected=false), 1:Ns)
                #Y[:,k] = std(reshape(Z[1:sx*Ns, k],Ns,sx),corrected=false,dims=2)
            end
        end

        return Y
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