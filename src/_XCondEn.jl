module _XCondEn
export XCondEn
using StatsBase: fit, Histogram
using Statistics: mean, std
    """
        XCond, SEw, SEz = XCondEn(Sig1, Sig2) 

    Returns the corrected cross-conditional entropy estimates (`XCond`) and the
    corresponding Shannon entropies (m: SEw, m+1: SEz) for m = [1,2] 
    estimated for the data sequences contained in `Sig1` and `Sig2` using the default
    parameters:  embedding dimension = 2, time delay = 1, number of symbols = 6, 
    logarithm = natural
    ** Note: XCondEn is direction-dependent. Therefore, the order of the
    data sequences `Sig1` and `Sig2` matters. If `Sig1` is the
    sequence 'y', and `Sig2` is the second sequence 'u', the `XCond` is
    the amount of information carried by y(i) when the pattern u(i) is found.**

        XCond, SEw, SEz = XCondEn(Sig1::Union{AbstractMatrix{T}, AbstractVector{T}} where T<:Real, Sig2::Union{AbstractVector{T} where T<:Real, Nothing} = nothing; m::Int=2, tau::Int=1, c::Int=6, Logx::Real=exp(1), Norm::Bool=false)

    Returns the corrected cross-conditional entropy estimates (`XCond`) for
    the data sequences contained in `Sig1` and `Sig2` using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, an integer > 1        [default: 2]  \n
    `tau`   - Time Delay, a positive integer             [default: 1]  \n
    `c`     - Number of symbols, an integer > 1          [default: 6]  \n
    `Logx`  - Logarithm base, a positive scalar          [default: natural] \n
    `Norm`  - Normalisation of `XCond` values: 
                [false]  no normalisation                  [default]\n
                [true]   normalises w.r.t cross-Shannon entropy.  \n

    # See also `XFuzzEn`, `XSampEn`, `XApEn`, `XPermEn`, `CondEn`, `XMSEn`

    # References:
        [1] Alberto Porta, et al.,
            "Conditional entropy approach for the evaluation of the 
            coupling strength." 
            Biological cybernetics 
            81.2 (1999): 119-129.
                       
    """
    function XCondEn(Sig1::Union{AbstractMatrix{T}, AbstractVector{T}} where T<:Real, Sig2::Union{AbstractVector{T} where T<:Real, Nothing} = nothing;
             m::Int=2, tau::Int=1,  c::Int=6, Logx::Real=exp(1), Norm::Bool=false)

    if all(isa.((Sig1,Sig2), AbstractVector))
        N = size(Sig1,1);   
        S1 = (Sig1 .- mean(Sig1))/std(Sig1,corrected=false)
        S2 = (Sig2 .- mean(Sig2))/std(Sig2,corrected=false) 
    elseif (minimum(size(Sig1))==2 && (Sig2 isa Nothing)) 
        argmin(size(Sig1)) == 2 ? nothing : Sig1 = Sig1'
        S1 = (Sig1[:,1] .- mean(Sig1[:,1]))/std(Sig1[:,1],corrected=false)
        S2 = (Sig2[:,2] .- mean(Sig2[:,2]))/std(Sig2[:,2],corrected=false)
        N = maximum(size(Sig1)); 
    else   error("""Sig1 and Sig2 must be 2 separate vectors 
        \t\t\t - OR - 
        Sig1 must be 2-column matrix and Sig2 nothing""")
    end

    length(S2)==N ? nothing : error("Sig1 and Sig2 must be the same length!")
    (N>=10) ? nothing :  error("Sig1/Sig2:   sequences must have >= 10 values")
    (m > 1) ? nothing : error("m:     must be an integer > 1")
    (tau>0) ? nothing : error("tau:   must be an integer > 0")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")
    (c>1) ? nothing :  error("c:     must be an integer > 1")

    Edges = range(minimum(S1),maximum(S1),length=c+1)
    Sx1 = map(x -> sum(Edges[1:c].<=x), S1)
    Edges = range(minimum(S2),maximum(S2),length=c+1)
    Sx2 = map(x -> sum(Edges[1:c].<=x), S2)

    SEw = zeros(m-1)
    SEz = zeros(m-1)
    Prcm = zeros(m-1)
    Xi = zeros(Int,N,m)
    for k = 1:m-1
        Nx = N-(k-1)*tau
        Xi[1:Nx,m-(k-1)] = Sx1[(k-1)*tau+1:N]
        #Wi = transpose(c.^(k-1:-1:0))*transpose(Xi[1:Nx,m-k+1:m])
        Wi = Xi[1:Nx,m-k+1:m]*(c.^(k-1:-1:0))
        #Zi = (c^k)*transpose(Sx2[(k-1)*tau+1:N]) .+ Wi 
        Zi = (c^k)*Sx2[(k-1)*tau+1:N]  + Wi 
        Pw = fit(Histogram, Wi[:], minimum(Wi)-.5:maximum(Wi)+.5).weights
        Pz = fit(Histogram, Zi[:], minimum(Zi)-.5:maximum(Zi)+.5).weights
        Prcm[k] = sum(Pz.==1)/Nx
        
        (sum(Pw)!= Nx || sum(Pz)!= Nx) ? @warn("Potential error estimating probabilities.") : nothing
        
        Pw = Pw[Pw.!=0];     Pw = Pw/N
        Pz = Pz[Pz.!=0];     Pz = Pz/N
        SEw[k] = -transpose(Pw)*log.(Logx, Pw)
        SEz[k] = -transpose(Pz)*log.(Logx, Pz)
    end

    Temp = fit(Histogram,Sx2,nbins=c).weights
    Temp = Temp[Temp.!=0]./N
    Sy = -transpose(Temp)*log.(Logx, Temp)
    XC = SEz - SEw + Prcm*Sy
    XC = vcat(Sy, XC)
    Norm ? XC = XC/Sy : nothing

    return XC, SEw, SEz
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