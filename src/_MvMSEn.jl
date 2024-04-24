module _MvMSEn
export MvMSEn
using Statistics: std, mean, median, var
using Plots
using DSP: conv
    """
        MSx, CI = MvMSEn(Data, Mobj) 

     Returns a vector of multivariate multiscale entropy values (`MSx`) and the complexity 
     index (`CI`) of the data sequences in `Data` using the parameters specified 
     by the multiscale object (`Mobj`) over 3 temporal scales with coarse-
     graining (default). 
         
    !!! note

        By default, the `MvSampEn` and `MvFuzzEn` multivariate entropy algorithms
        estimate entropy values using the "full"  method by comparing delay vectors 
        across all possible `m+1` expansions of the embedding space as applied in [1].
        These methods are not lower-bounded to 0, like most entropy algorithms,
        so `MvMSEn` may return negative entropy values if the base multivariate 
        entropy function is `MvSampEn` and `MvFuzzEn`, even for stochastic processes...

    -------------------------------------------------------------
    
        MSx, CI = MSEn(Data::AbstractArray{T,2} where T<:Real, Mobj::NamedTuple; Scales::Int=3, Methodx::String="coarse", Plotx::Bool=false)

     Returns a vector of multivariate multiscale entropy values (`MSx`) and the complexity 
     index (`CI`) of the data sequences in `Data` using the parameters specified by
     the multiscale object (`Mobj`) and the following keyword arguments:

     # Arguments:
     `Scales`   - Number of temporal scales, an integer > 1   (default: 3) \n
     `Method`   - Graining method, one of the following: \n
                  {`coarse`,`modified`,`generalized`} [default = `coarse`]  
                  For further info on these graining procedures, see the EntropyHub guide.  \n
     `Plotx`    - When Plotx == true, returns a plot of the entropy value at each
                  time scale (i.e. the multiscale entropy curve) [default: false]

    !!! tip 
    
        For further info on these graining procedures see the EntropyHub guide.

     # See also  `MSobject`, `cMvMSEn`, `MvFuzzEn`, `MvSampEn`, `MvPermEn`, `MvCoSiEn`,  `MvDispEn`
    
     # References:
        [1] Ahmed Mosabber Uddin, Danilo P. Mandic
             "Multivariate multiscale entropy analysis."
             IEEE signal processing letters 19.2 (2011): 91-94.
        
        [2] Madalena Costa, Ary Goldberger, and C-K. Peng,
             "Multiscale entropy analysis of complex physiologic time series."
             Physical review letters
             89.6 (2002): 068102.
        
        [3] Vadim V. Nikulin, and Tom Brismar,
             "Comment on “Multiscale entropy analysis of complex physiologic
             time series”." 
             Physical Review Letters 
             92.8 (2004): 089803.
        
        [4] Madalena Costa, Ary L. Goldberger, and C-K. Peng. 
             "Costa, Goldberger, and Peng reply." 
             Physical Review Letters
             92.8 (2004): 089804.
        
        [5] Madalena Costa, Ary L. Goldberger and C-K. Peng,
             "Multiscale entropy analysis of biological signals." 
             Physical review E 
             71.2 (2005): 021906.
        
        [6] Ranjit A. Thuraisingham and Georg A. Gottwald,
             "On multiscale entropy analysis for physiological data."
             Physica A: Statistical Mechanics and its Applications
             366 (2006): 323-332.
        
        [7] Ahmed Mosabber Uddin, Danilo P. Mandic
             "Multivariate multiscale entropy: A tool for complexity
             analysis of multichannel data."
             Physical Review E 84.6 (2011): 061918.


    """
    function MvMSEn(Data::AbstractArray{T,2} where T<:Real, Mobj::NamedTuple; Scales::Int=3, Methodx::String="coarse", Plotx::Bool=false)

    N, Dn = size(Data)   
    (N>10) && (Dn>1) ? nothing : error("Data:   must be an NxM matrix where N>10 and M>1")
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
                with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (lowercase(Methodx) in ["coarse","modified","generalized"]) ? nothing :
        error("Method:  must be one of the following string names - 'coarse','modified','generalized'")                
    String(Symbol(Mobj.Func))[1:2] != "Mv" ? error("Base entropy estimator must be a multivariate entropy method. ",
    "To perform univariate multiscale entropy estimation, use MSEn().") : nothing    

    MSx = zeros(Scales)
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
    Func2 = getfield(_MvMSEn,Symbol(lowercase(Methodx)))
    for T = 1:Scales
        print(". ")
        Temp = Func2(Data,T,Dn)         
        Tempx = Mobj.Func(Temp; Args...)
        typeof(Tempx)<:Tuple ? Temp2 = mean(Tempx[1]) : Temp2 = mean(Tempx)
        MSx[T] = Temp2
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
        title = "Multivariate Multiscale $(string(Mobj.Func)[3:end]) ($(titlecase(Methodx))-graining method)",
        plot_titlefontsize=16, plot_titlefontcolor=RGB(7/255, 54/255, 66/255)) #ylim=(0,maximum(MSx)+.2),
        display(p1)
    end

    return MSx, CI
    end


    function coarse(Z, sx, Dn)
        Ns = Int(floor(size(Z,1)/sx))
        Y = zeros(Ns,Dn)
        for k in 1:Dn
            Y[:,k] = mean(reshape(Z[1:sx*Ns,k],sx,Ns),dims=1)
        end
        return Y
    end

    function modified(Z, sx, Dn)
        Y = (conv(Z,ones(Int,sx))/sx)[sx:end-sx+1,:]
        return Y 
    end

    function generalized(Z, sx, Dn)
        Ns = Int(floor(size(Z,1)/sx))
        Y = zeros(Ns,Dn)
        for k in 1:Dn
            Y[:,k] = var(reshape(Z[1:sx*Ns, k],sx,Ns)',corrected=false,dims=2)
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