module _IncrEn
export IncrEn
using Statistics: std
    """
        Incr = IncrEn(Sig) 

    Returns the increment entropy (`Incr`) estimate of the data sequence 
    (`Sig`) using the default parameters: 
    embedding dimension = 2, time delay = 1, quantifying resolution = 4,
    logarithm = base 2,

        Incr = IncrEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, R::Int=4, Logx::Real=2, Norm::Bool=false)

    Returns the increment entropy (`Incr`) estimate of the data sequence
    (`Sig`) using the specified 'keyword' arguments:

    # Arguments:    
    `m`     - Embedding Dimension, an integer > 1   \n
    `tau`   - Time Delay, a positive integer    \n
    `R`     - Quantifying resolution, a positive scalar    \n
    `Logx`  - Logarithm base, a positive scalar (enter 0 for natural log) \n
    `Norm`  - Normalisation of IncrEn value: \n
              [false]  no normalisation - default
              [true]   normalises w.r.t embedding dimension (m-1). 

    # See also `PermEn`, `SyDyEn`, `MSEn`

    # References:
        [1] Xiaofeng Liu, et al.,
            "Increment entropy as a measure of complexity for time series."
            Entropy
            18.1 (2016): 22.1.

        ***   "Correction on Liu, X.; Jiang, A.; Xu, N.; Xue, J. - Increment 
            Entropy as a Measure of Complexity for Time Series,
            Entropy 2016, 18, 22." 
            Entropy 
            18.4 (2016): 133.

        [2] Xiaofeng Liu, et al.,
            "Appropriate use of the increment entropy for 
            electrophysiological time series." 
            Computers in biology and medicine 
            95 (2018): 13-23.


    """
    function IncrEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, 
        R::Int=4, Logx::Real=2, Norm::Bool=false)

    Logx == 0  ? Logx = exp(1) : nothing
        
    (size(Sig,1) > 10) ? nothing :  error("Sig:   must be a numeric vector")
    (m > 1) ? nothing :  error("m:     must be an integer > 1")
    (tau>0) ? nothing :  error("tau:   must be an integer > 0")
    (R > 0) ? nothing :  error("R:     must be a positive integer > 0")
    (Logx>0) ? nothing : error("Logx:  must be a positive number > 0")
    
    Vi = diff(Sig)
    N = size(Vi,1)-((m-1)*tau)
    Vk = zeros(N,m)
    for k = 1:m
        Vk[:,k] = Vi[1+(k-1)*tau:N+(k-1)*tau]
    end

    Sk = sign.(Vk)
    Temp = std(Vk,dims=2)[:]
    Qk = min.(R, floor.((abs.(Vk)*R)./repeat(Temp,outer=(1,m))))
    Qk[any(Temp.==0,dims=2), :] .= 0  #should that be all()
    Wk = Sk.*Qk  
    Wk[Wk.==-0] .= 0
    Px = unique(Wk,dims=1)
    Counter = zeros(Int,size(Px,1));
    for k = 1:size(Px,1) 
        Counter[k] = sum(all(Wk .- transpose(Px[k,:]) .==0 ,dims=2))
    end
    Ppi = Counter/N

    if size(Px,1) > (2*R + 1)^m
        @warn("Error with probability estimation'")
    elseif round(sum(Ppi),digits=5) != 1
        @warn("Error with probability estimation")
    end
    Incr = -sum(Ppi.*(log.(Logx, Ppi)))
    if Norm
        Incr = Incr/(m-1);
    end
   
    return Incr
    end

end

"""
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