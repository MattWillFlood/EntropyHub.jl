module _SyDyEn
export SyDyEn
using Clustering: kmeans, assignments
using StatsBase: Histogram, fit
using Statistics: mean
    """
        SyDy, Zt = SyDyEn(Sig) 

    Returns the symbolic dynamic entropy (`SyDy`) and the symbolic sequence
    (`Zt`) of the data sequence (`Sig`) using the default parameters: 
    embedding dimension = 2, time delay = 1, symbols = 3, logarithm = natural,
    symbolic partition type = maximum entropy partitioning (`MEP`), 
    normalisation = normalises w.r.t # possible vector permutations (c^m) 

        SyDy, Zt = SyDyEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, c::Int=3, Typex::String="MEP", Logx::Real=exp(1), Norm::Bool=true)

    Returns the symbolic dynamic entropy (`SyDy`) and the symbolic sequence
    (`Zt`) of the data sequence (`Sig`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer  \n
    `tau`   - Time Delay, a positive integer  \n
    `c`     - Number of symbols, an integer > 1  \n
    `Typex` - Type of symbolic sequnce partitioning method, one of the following:  \n
              {"linear","uniform","MEP"(default),"kmeans"}    
    `Logx`  - Logarithm base, a positive scalar    \n
    `Norm`  - Normalisation of SyDyEn value: \n 
              [false]  no normalisation 
              [true]   normalises w.r.t # possible vector permutations (c^m+1) - default

    See the EntropyHub guide for more info on these parameters.

    # See also `DispEn`, `PermEn`, `CondEn`, `SampEn`, `MSEn`
  
    # References:
        [1] Yongbo Li, et al.,
            "A fault diagnosis scheme for planetary gearboxes using 
            modified multi-scale symbolic dynamic entropy and mRMR feature 
            selection." 
            Mechanical Systems and Signal Processing 
            91 (2017): 295-312. 
  
        [2] Jian Wang, et al.,
            "Fault feature extraction for multiple electrical faults of 
            aviation electro-mechanical actuator based on symbolic dynamics
            entropy." 
            IEEE International Conference on Signal Processing, 
            Communications and Computing (ICSPCC), 2015.
  
        [3] Venkatesh Rajagopalan and Asok Ray,
            "Symbolic time series analysis via wavelet-based partitioning."
            Signal processing 
            86.11 (2006): 3309-3320.
  
    """
    function SyDyEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, 
        c::Int=3, Typex::String="MEP", Logx::Real=exp(1), Norm::Bool=true)

    N = size(Sig,1)
    (N > 10) ?  nothing : error("Sig:   must be a numeric vector")
    (m > 0) ?  nothing :  error("m:     must be an integer > 0")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (c > 1) ? nothing :   error("c:     must be an integer > 1")
    (lowercase(Typex) in ["linear", "kmeans", "uniform","mep"]) ? nothing :
        error("Typex:    must be one of the following strings - 'linear','kmeans','uniform','MEP'")    
    (Logx>0) ? nothing :  error("Logx:  must be a positive number > 0")

    Nx = N-((m-1)*tau)
    Zt = zeros(N)
    if lowercase(Typex) ==  "linear"
        Edges = range(minimum(Sig),maximum(Sig),length=c+1)
        Zt = map(x -> sum(Edges[1:c].<=x), Sig)
    elseif lowercase(Typex) == "uniform"
        Ix = sortperm(Sig)
        Edges = range(1,N,length=c+1)
        z = map(x -> sum(Edges[1:c].<=x), 1:N)
        Zt[Ix] .= z 
            
    elseif lowercase(Typex) == "kmeans"
        Tx = kmeans(transpose(Sig), c; maxiter=200)
        z = Int.(assignments(Tx))
        ix = sortperm(Tx.centers[:]);   Zt = zeros(Int,N)
        for k = 1:c
            Zt[z.==ix[k]] .= k;
        end
            
    else 
        Tx = sort(Sig)
        Edges = Tx[vcat(1, Int.(ceil.((1:c-1)*N/c)),N)]
        Zt = map(x -> sum(Edges[1:end-1].<=x), Sig)
    end

    Zm = zeros(Int,Nx,m);
    for n = 1:m
        Zm[:,n] = Zt[(n-1)*tau + 1:Nx+(n-1)*tau]
    end

    T = unique(Zm,dims=1)
    Counter = zeros(size(T,1))
    Counter2 = zeros(size(T,1),c)
    Bins = range(0.5,c+.5,step=1)
    for n = 1:size(T,1)
        Ordx = all(Zm .- transpose(T[n,:]).==0,dims=2)
        Counter[n] = mean(Ordx)
        Temp = Zm[vcat(falses(m*tau), Ordx[1:end-(m*tau)]),1]
        Counter2[n,:] = fit(Histogram,Temp,Bins).weights
    end
    Counter2 ./= sum(Counter2,dims=2)
    Counter2[isnan.(Counter2)] .= 0

    P1 = -sum(Counter.*log.(Logx, Counter))
    P2 = log.(Logx, repeat(Counter,outer=(1,c)).*Counter2)
    P2[isinf.(P2)] .= 0
    SyDy = P1 .- transpose(Counter)*(sum(P2,dims=2))

    if round(sum(Counter),digits=4) != 1 || maximum(round.(sum(Counter2,dims=2),digits=4)) != 1
        print(maximum(sum(Counter2,dims=2)))
        @warn("Potential Error calculating probabilities")
    end
    if Norm
        SyDy = SyDy/(log(Logx, c^(m+1)))
    end

    return SyDy, Zt
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