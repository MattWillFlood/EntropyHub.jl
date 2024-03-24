module _SlopEn
export SlopEn
using GroupSlices
    """
        Slop = SlopEn(Sig) 

    Returns the slope entropy (`Slop`) estimates for embedding dimensions
    [2, ..., m] of the data sequence (`Sig`) using the default parameters:
    embedding dimension = 2, time delay = 1, 
    angular thresholds = [5 45],  logarithm = base 2 

        Slop = SlopEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, Lvls::AbstractArray{T,1} where T<:Real=[5, 45], Logx::Real=2, Norm::Bool=true)

    Returns the slope entropy (`Slop`) estimate of the data sequence (`Sig`)  
    using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, an integer > 1   	\n
              SlopEn returns estimates for each dimension [2,...,m]
    `tau`   - Time Delay, a positive integer    \n
    `Lvls`  - Angular thresolds, a vector of monotonically increasing   
              values in the range [0 90] degrees.\n
    `Logx`  - Logarithm base, a positive scalar (enter 0 for natural log)   \n
    `Norm`  - Normalisation of SlopEn value, a boolean operator: \n
              [false]  no normalisation
              [true]   normalises w.r.t. the number of patterns found (default)

    # See also `PhasEn`, `GridEn`, `MSEn`, `CoSiEn`, `SampEn`, `ApEn`

    # References:
        [1] David Cuesta-Frau,
            "Slope Entropy: A New Time Series Complexity Estimator Based on
            Both Symbolic Patterns and Amplitude Information." 
            Entropy 
            21.12 (2019): 1167.

    """
    function SlopEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1,
         Lvls::AbstractArray{T,1} where T<:Real=[5, 45], Logx::Real=2, Norm::Bool=true)
        
    Logx == 0 ? Logx = exp(1) : nothing
    (size(Sig,1) >10) ? nothing :  error("Sig:   must be a numeric vector")
    (m > 1) ? nothing :  error("m:     must be an integer > 1")
    (tau>0) ? nothing :  error("tau:   must be an integer > 0")
    (length(Lvls)>1 && all(diff(Lvls).>0) && all(0 .< Lvls .< 90)) ? nothing :
        error("Lvls:    must be a vector of 2 or more monotonically increasing
                        values in the range [0 90] degrees")
    (Logx>0) ? nothing : error("Logx:  must be a positive number > 0")

    m = m-1;
    Tx = atand.(Sig[1+tau:end] .- Sig[1:end-tau])
    N = size(Tx,1)
    Sx = zeros(Int,N,m)
    Symbx = zeros(Int,size(Tx));
    Slop = zeros(m)
    sort!(Lvls)

    for q = 2:length(Lvls)
        Symbx[(Tx.<= Lvls[q]) .& (Tx .> Lvls[q-1])] .= q-1
        Symbx[(Tx.>=-Lvls[q]) .& (Tx .<-Lvls[q-1])] .= -(q-1)
        
        if q == length(Lvls)
            Symbx[Tx.> Lvls[q]] .= q
            Symbx[Tx.<-Lvls[q]] .= -q
        end
    end

    for k = 1:m
        Sx[1:N-k+1,k] = Symbx[k:N]
        Locs = groupslices(Sx[1:N-k+1,1:k])
        p = []
        [push!(p, sum(Locs.==n)) for n in unique(Locs)]
            
        Norm ? p ./=(N-k+1) : p./= length(p)

        if Norm && round(sum(p)) != 1
            @warn("Potential Error: Some permutations not accounted for!")
            print(round(sum(p)))
        end
        
        Slop[k] = -sum(p.*log.(Logx, p))
    end

    return Slop
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