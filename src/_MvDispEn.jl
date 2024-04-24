module _MvDispEn
export MvDispEn
using Clustering: kmeans, assignments
using Statistics: std, mean
using StatsFuns: normcdf
using Combinatorics: combinations
"""
    MDisp, RDE = MvDispEn(Data) 
    
 Returns the multivariate dispersion entropy estimate (`MDisp`) and
 the reverse dispersion entropy (`RDE`) for the M multivariate sequences 
 in `Data` using the default parameters:
 embedding dimension = 2*ones(M,1), time delay = ones(M,1), # symbols = 3, 
 algorithm method = "v1" (see below), data transform = normalised cumulative density function (ncdf)
 logarithm = natural, data normalization = true,

 ________________________________________________________________________ 
 
!!! note 

    By default, `MvDispEn` uses the method termed `mvDEii` in [1],
    which follows the original multivariate embedding approach of Ahmed & Mandic [2].
    The `v1` method therefore returns a singular entropy estimate.

    If the `v2` method is selected (`Methodx=="v2"`), the main method
    outlined in [1] termed `mvDE` is applied. In this case, entropy is estimated
    using each combination of multivariate delay vectors with lengths 1:max(m),
    with each entropy value returned accordingly. See [1] for more info.

 ________________________________________________________________________ 
      
    MDisp, RDE = MvDispEn(Data::AbstractArray{T,2} where T<:Real; m::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, tau::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, c::Int=3, Methodx::String="v1", Typex::String="NCDF", Norm::Bool=false, Logx::Real=exp(1))

 Returns the multivariate dispersion entropy estimate (`MDisp`) for the M
 multivariate data sequences in `Data` using the specified keyword arguments:

 # Arguments:
 `Data`    - Multivariate dataset, NxM matrix of N (>10) observations (rows) and M (cols) univariate data sequences\n
 `m`       - Embedding Dimension, a vector of M positive integers\n
 `tau`     - Time Delay, a vector of M positive integers\n
 `c`       - Number of symbols in transform, an integer > 1 \n
 `Methodx` - The method of multivariate dispersion entropy estimation as outlined in [1], either:\n
        * `"v1"` - employs the method consistent with the original multivariate embedding approach of Ahmed &
                    Mandic [2], termed `mvDEii` in [1]. (default)
        * `"v2"` - employs the main method derived in [1],  termed `mvDE`.
 `Typex`   - Type of data-to-symbolic sequence transform, one of the following:\n
                {`'linear'`, `'kmeans'`, `'ncdf'`, `'equal'`}
            See the `EntropyHub Guide` for more info on these transforms.
 `Norm`    - Normalisation of `MDisp` and `RDE` values, a boolean:\n
                * [false]   no normalisation (default)
                * [true]    normalises w.r.t number of possible dispersion patterns (`c^m`).
 `Logx`   - Logarithm base, a positive scalar \n 
        
 # See also   `DispEn`, `DispEn2D`, `MvSampEn`, `MvFuzzEn`, `MvPermEn`, `MSEn`
 
 # References:
    [1] H Azami, A Fernández, J Escudero
          "Multivariate Multiscale Dispersion Entropy of Biomedical Times Series"
          Entropy 2019, 21, 913.

    [2] Ahmed Mosabber Uddin, Danilo P. Mandic
          "Multivariate multiscale entropy: A tool for complexity
          analysis of multichannel data."
          Physical Review E 84.6 (2011): 061918.

    [3] Mostafa Rostaghi and Hamed Azami,
           "Dispersion entropy: A measure for time-series analysis." 
           IEEE Signal Processing Letters 
           23.5 (2016): 610-614.

    [4] Hamed Azami and Javier Escudero,
           "Amplitude-and fluctuation-based dispersion entropy." 
           Entropy 
           20.3 (2018): 210.

    [5] Li Yuxing, Xiang Gao and Long Wang,
           "Reverse dispersion entropy: A new complexity measure for sensor signal." 
           Sensors 
           19.23 (2019): 5203.
"""
    function MvDispEn(Data::AbstractArray{T,2} where T<:Real; m::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, tau::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, 
                c::Int=3, Methodx::String="v1", Typex::String="NCDF", Norm::Bool=false, Logx::Real=exp(1))

    N, Dn = size(Data)
    isnothing(m) ? m = 2*ones(Int, Dn) : nothing
    isnothing(tau) ? tau = ones(Int, Dn) : nothing
    Logx==0 ? Logx = exp(1) : nothing

    (N>10) && (Dn>1) ? nothing : error("Data:   must be an NxM matrix where N>10 and M>1")
    ndims(m)==1  && length(m)==Dn  && all(m.>0) && eltype(m)<:Int ? nothing : error("m:  vector of M positive integers")
    ndims(tau)==1  && length(tau)==Dn  && all(tau.>0) && eltype(tau)<:Int ? nothing : error("tau:  vector of M positive integers")
    (c > 1) ? nothing :   error("c:     must be an integer > 1")
    (lowercase(Typex) in ["linear", "kmeans", "ncdf","equal"]) ? nothing :
                 error("Typex:    must be one of the following strings - 'linear','kmeans','ncdf','equal'") 
    (lowercase(Methodx) in ["v1", "v2"]) ? nothing :  error("Methodx:    must be either 'v1' or 'v2'") 
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")

    Sx = zeros(Int8, N, Dn)               
    for q in eachindex(1:Dn)
        Sig = Data[:,q]  

        if lowercase(Typex) == "linear"    
            Edges = range(minimum(Sig),maximum(Sig),length=c+1)
            Zi = map(x -> sum(Edges[1:c].<=x), Sig)
    
        elseif lowercase(Typex) == "kmeans"
            Temp = kmeans(transpose(Sig), c; maxiter=200)
            Zx = assignments(Temp)
            Clux = Temp.centers
            xx = sortperm(Clux[:]);   Zi = zeros(N)
            for k = 1:c
                Zi[Zx.==xx[k]] .= k;
            end
                
        elseif lowercase(Typex) == "ncdf"
            Zx = normcdf.(mean(Sig),std(Sig,corrected=false),Sig);
            Zi = map(x -> sum(range(0,1,length=c+1)[1:c].<=x), Zx)
        
        elseif lowercase(Typex) == "equal"
            ix = sortperm(Sig,alg=MergeSort);
            xx = Int.(round.(range(0,N,length=c+1)))
            Zi = zeros(N)
            for k = 1:c
                Zi[ix[xx[k]+1:xx[k+1]]] .= k
            end
        end

        Sx[:,q] = Zi[:]
    end

    Nx = N-maximum((m.-1).*tau)
    Vex = zeros(Int8, Nx, sum(m))
    q=1
    for k = 1:Dn
        for p = 1:m[k]
            Vex[:,q] = Sx[1+(p-1)*tau[k]:Nx+(p-1)*tau[k],  k]
            q += 1
        end    
    end

    if  lowercase(Methodx) == "v1"
        Px = unique(Vex,dims=1)      
        Counter = map(n -> sum(all((Vex .- Px[n:n,:]).==0, dims=2)), 1:size(Px,1))
        Counter = Counter[Counter.!=0]
        Ppi = Counter/sum(Counter)        
        if round(sum(Ppi),digits=5) != 1
            @warn ("Potential error with probability calculation")
        end

        MDisp = -sum(Ppi.*log.(Logx,Ppi))
        RDE = sum((Ppi .- (1/(c^sum(m)))).^2)            
        if Norm
            MDisp = MDisp/log(Logx, c^sum(m))
            RDE = RDE/(1-(1/(c^sum(m))))
        end
        #return MDisp, RDE
                
    elseif lowercase(Methodx) == "v2"
        P = sum(m)
        MDisp = zeros(maximum(m))
        RDE = zeros(maximum(m))
        for k in eachindex(1:maximum(m))
            print(" . ")
            Temp = collect(combinations(1:P,k))
            Vez = zeros(Int8, Nx*binomial(P,k),k)
            for q in eachindex(1:size(Temp,1))
                # Vez[q*Nx:(q+1)*Nx-1,:] = Vex[:,Temp[q][:]] 
                Vez[1+(q-1)*Nx:q*Nx,:] = Vex[:,Temp[q]] 
            end
            Px = unique(Vez, dims=1) 
            Counter = map(n -> sum(all((Vez .- Px[n:n,:]).==0,dims=2)), 1:size(Px,1))
            Counter = Counter[Counter.!=0]
            Ppi = Counter/sum(Counter)
        
            if round(sum(Ppi),digits=5) != 1
                @warn ("Potential error with probability calculation")
            end
            MDisp[k] = -sum(Ppi.*log.(Logx,Ppi))
            RDE[k] = sum((Ppi .- (1/(c^k))).^2)
         
            if Norm
                MDisp[k] = MDisp[k]/log(Logx, c^k)
                RDE[k] = RDE[k]/(1-(1/(c^k)))
            end
        end    
    end
    return MDisp, RDE                

end        
end

"""Copyright 2024 Matthew W. Flood, EntropyHub
  
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

For Terms of Use see https://github.com/MattWillFlood/EntropyHub"""