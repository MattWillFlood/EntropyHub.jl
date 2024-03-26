module _CondEn
export CondEn
using Statistics: std, mean
using StatsBase: Histogram, fit
    """
        Cond, SEw, SEz = CondEn(Sig) 

    Returns the corrected conditional entropy estimates (`Cond`) and the
    corresponding Shannon entropies (m: `SEw`, m+1: `SEz`) for m = [1,2] 
    estimated from the data sequence (`Sig`) using the default  parameters:
    embedding dimension = 2, time delay = 1, symbols = 6, logarithm = natural,
    normalisation = false
    *Note: CondEn(m=1) returns the Shannon entropy of `Sig`.*

        Cond, SEw, SEz = CondEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, c::Int=6, Logx::Real=exp(1), Norm::Bool=false)

    Returns the corrected conditional entropy estimates (`Cond`) from the data
    sequence (`Sig`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, an integer > 1  \n
    `tau`   - Time Delay, a positive integer  \n
    `c`     - # of symbols, an integer > 1  \n
    `Logx`  - Logarithm base, a positive scalar  \n 
    `Norm`  - Normalisation of CondEn value:  \n
              [false]  no normalisation - default
              [true]   normalises w.r.t Shannon entropy of data sequence `Sig`  

    # See also `XCondEn`, `MSEn`, `PermEn`, `DistEn`, `XPermEn`
  
    # References:
        [1] Alberto Porta, et al.,
            "Measuring regularity by means of a corrected conditional
            entropy in sympathetic outflow." 
            Biological cybernetics 
            78.1 (1998): 71-78.
            
    """
    function CondEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, 
        c::Int=6, Logx::Real=exp(1), Norm::Bool=false)
        
    (size(Sig)[1] > 10) ? nothing :   error("Sig:   must be a numeric vector")
    (m > 1) ? nothing :  error("m:     must be an integer > 1")
    (tau > 0) ? nothing :  error("tau:   must be an integer > 0")
    (Logx > 0) ? nothing : error("Logx:     must be a positive number > 0")
    (c > 1) ? nothing :  error("c:     must be an integer > 1")

    Edges = range(minimum(Sig),maximum(Sig),length=c+1)
    Sx = map(x -> sum(Edges[1:c].<=x), Sig)

    N = size(Sx)[1]
    SEw = zeros(m-1)
    SEz = zeros(m-1)
    Prcm = zeros(m-1)
    Xi = zeros(N,m)
    Xi[:,m] = Sx
    for k = 1:m-1
        Nx = N-(k*tau)
        Xi[1:Nx,end-k] = Sx[(k*tau)+1:N]
        Wi = Xi[1:Nx,m-k+1:m] * (c.^collect(k-1:-1:0)) # Maybe dot notation here???
        Zi = Xi[1:Nx,m-k:m] * (c.^collect(k:-1:0))

        Pw = fit(Histogram, Wi, minimum(Wi)-.5:maximum(Wi)+.5).weights
        Pz = fit(Histogram, Zi, minimum(Zi)-.5:maximum(Zi)+.5).weights
        Prcm[k] = sum(Pw.==1)/Nx
        
        if sum(Pw)!= Nx || sum(Pz)!= Nx
            @warn("Potential error estimating probabilities.")
        end
        
        Pw = Pw[Pw.!=0];    Pw /= N;
        Pz = Pz[Pz.!=0];    Pz /= N;
        SEw[k] = -transpose(Pw)*log.(Logx, Pw)
        SEz[k] = -transpose(Pz)*log.(Logx, Pz)
        
    end

    Temp = fit(Histogram,Sx,.5:c+.5).weights/N;
    Temp = Temp[Temp.!=0]
    S1 = -transpose(Temp)*log.(Logx, Temp)
    Cond = SEz - SEw + Prcm*S1;
    Cond = vcat(S1, Cond);
    if Norm
        Cond = Cond/S1;
    end

    return Cond, SEw, SEz
    end

end

    #= Sig = (Sig.-mean(Sig))./std(Sig,corrected=false);
    Edges = range(minimum(Sig),maximum(Sig),length=c+1)
    #Edges[1] -= .1; Edges[end] += .1
    Sx = map(x -> searchsortedfirst(Edges,x), Sig) .- 1
    Sx[Sx.==0] .= 1 =#

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