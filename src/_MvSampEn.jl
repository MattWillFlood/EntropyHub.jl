module _MvSampEn
export MvSampEn
using Statistics: std
# using LinearAlgebra: UpperTriangular, I
"""
    MSamp, B0, Bt, B1 = MvSampEn(Data) 

 Returns the multivariate sample entropy estimate (`MSamp`) and the
 average number of matched delay vectors (`m`: `B0`; joint total 
 `m+1` subspace: `Bt`; all possible `m+1` subspaces: `B1`),
 from the M multivariate sequences in `Data` using the default parameters: 
 embedding dimension = 2*ones(M), time delay = ones(M), radius threshold = 0.2,
 logarithm = natural, data normalization = false
 
 
!!! note 

    The entropy value returned as `MSamp` is estimated using the "full" 
    method [i.e.  -log(Bt/B0)] which compares delay vectors across all possible `m+1` 
    expansions of the embedding space as applied in [1][2]. Contrary to
    conventional definitions of sample entropy, this method does not provide a
    lower bound of 0!!
    Thus, it is possible to obtain negative entropy values for multivariate 
    sample entropy, even for stochastic processes...

    Alternatively, one can calculate `MSamp` via the "naive" method, 
    which ensures a lower bound of 0, by using the average number of matched
    vectors for an individual `m+1` subspace (B1) [e.g. -log(B1(1)/B0)],
    or the average for all `m+1` subspaces [i.e. -log(mean(B1)/B0)].

    To maximize the number of points in the embedding process, this algorithm 
    uses N - max(m * tau) delay vectors and _**not**_ N-max(m) * max(tau) as employed 
    in [1][2].
 
-------------------------------------------------------------
      
    MSamp, B0, Bt, B1 = MvSampEn(Data::AbstractArray{T} where T<:Real; m::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, tau::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, r::Real=0.2, Logx::Real=exp(1), Norm::Bool=false)
        
 Returns the multivariate sample entropy estimates (`MSamp`) estimated
 from the M multivariate data sequences in `Data` using the specified 
 keyword arguments:

 # Arguments:
 `Data`  - Multivariate dataset, NxM matrix of N (>10) observations (rows) and M (cols) univariate data sequences\n
 `m`     - Embedding Dimension, a vector of M positive integers\n
 `tau`   - Time Delay, a vector of M positive integers\n
 `r`     - Radius Distance Threshold, a positive scalar  \n
 `Logx`  - Logarithm base, a positive scalar \n 
 `Norm`  - Normalisation of all M sequences to unit variance, a boolean\n
        
 # See also `SampEn`, `XSampEn`, `SampEn2D`, `MSEn`, `MvFuzzEn`, `MvPermEn`

 # References:
    [1] Ahmed Mosabber Uddin, Danilo P. Mandic
        "Multivariate multiscale entropy: A tool for complexity
        analysis of multichannel data."
        Physical Review E 84.6 (2011): 061918.

    [2] Ahmed Mosabber Uddin, Danilo P. Mandic
        "Multivariate multiscale entropy analysis."
        IEEE signal processing letters 19.2 (2011): 91-94.
"""
    function MvSampEn(Data::AbstractArray{T,2} where T<:Real; m::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, tau::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, 
                        r::Real=0.2, Logx::Real=exp(1), Norm::Bool=false)

    N, Dn = size(Data)
    isnothing(m) ? m = 2*ones(Int, Dn) : nothing
    isnothing(tau) ? tau = ones(Int, Dn) : nothing
    Norm ? Data = Data./std(Data, dims=1, corrected=false) : nothing

    (N>10) && (Dn>1) ? nothing : error("Data:   must be an NxM matrix where N>10 and M>1")
    ndims(m)==1  && length(m)==Dn  && all(m.>0) && eltype(m)<:Int ? nothing : error("m:  vector of M positive integers")
    ndims(tau)==1  && length(tau)==Dn  && all(tau.>0) && eltype(tau)<:Int ? nothing : error("tau:  vector of M positive integers")
    (r>=0) ? nothing : error("r:     must be a positive scalar value")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")

    Nx = N - maximum((m.-1).*tau)
    Ny = N - maximum(m.*tau)
    Vex = zeros(Nx,sum(m))
    q = 1;
    for k = 1:Dn
        for p = 1:m[k]
            Vex[:,q] = Data[1+(p-1)*tau[k]:Nx+(p-1)*tau[k],  k]
            q += 1
        end    
    end
    Count0 = Distx(Vex,r)
    B0 = sum(Count0)/(Nx*(Nx-1)/2)

    B1 = zeros(Dn)
    Temp = cumsum(m)
    Vez = Inf.*ones(1,sum(m)+1)
    for k = 1:Dn
        Sig = Data[1+m[k]*tau[k]:Ny+m[k]*tau[k], k]
        Vey = hcat(Vex[1:Ny, 1:Temp[k]], Sig,  Vex[1:Ny, Temp[k]+1:end])
        Vez = vcat(Vez, Vey)
        Count1 = Distx(Vey, r);
        B1[k] = sum(Count1)/(Ny*(Ny-1)/2)
    end
    Vez = Vez[2:end,:]
    Count1 = Distx(Vez, r);
    Bt = sum(Count1)/(Dn*Ny*((Dn*Ny)-1)/2)
    
    MSamp = -log(Bt/B0)/log(Logx)
    return MSamp, B0, Bt, B1
    end


    function Distx(Vex, r)
        Nt = size(Vex)[1]
        Counter = zeros(Bool, Nt-1,Nt-1)
        for x=1:Nt-1
            Counter[x,x:end] = all(abs.(Vex[x+1:end,:] .- Vex[x:x,:]) .<= r, dims=2)
        end
        return Counter
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