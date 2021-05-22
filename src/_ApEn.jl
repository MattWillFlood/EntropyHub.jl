module _ApEn
export ApEn
using Statistics: mean, std
    """
        Ap, Phi = ApEn(Sig)

    Returns the approximate entropy estimates `Ap` and the log-average number of 
    matched vectors `Phi` for `m` = [0,1,2], estimated from the data sequence `Sig`
    using the default parameters: embedding dimension = 2, time delay = 1,
    radius distance threshold = 0.2*SD(`Sig`), logarithm = natural

        Ap, Phi = ApEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1))

    Returns the approximate entropy estimates `Ap` of the data sequence `Sig`
    for dimensions = [0,1,...,`m`] using the specified keyword arguments:

    # Arguments:        
    `m`      - Embedding Dimension, a positive integer\n
    `tau`    - Time Delay, a positive integer\n
    `r`      - Radius Distance Threshold, a positive scalar  \n
    `Logx`   - Logarithm base, a positive scalar\n  

    # See also `XApEn`, `SampEn`, `MSEn`, `FuzzEn`, `PermEn`, `CondEn`, `DispEn`
  
    # References:
        [1] Steven M. Pincus, 
            "Approximate entropy as a measure of system complexity." 
            Proceedings of the National Academy of Sciences 
            88.6 (1991): 2297-2301.
  
    """
    function ApEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1))

    N = length(Sig)
    (N>10) ? nothing : error("Sig:   must be a numeric vector with >10 samples")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (r>=0) ? nothing :  error("r:     must be a positive value")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")

    Counter = 1*(abs.(Sig .- transpose(Sig)) .<= r)  
    M = Int.([m*ones(N-m*tau); repeat(collect(m-1:-1:1),inner=tau)])
    Ap = zeros(m+1)
    Phi = zeros(m+2)

    for n = 1:N-tau
        ix = findall(Counter[n, :] .== 1)
        
        for k = 1:M[n]
            ix = ix[(ix .+ (k*tau)) .<= N]     
            p1 = repeat(transpose(Sig[n:tau: n+(tau*k)]), length(ix))                       
            p2 = Sig[ix .+ transpose(collect(0:tau:(k*tau)))]           
            ix = ix[findall(maximum(abs.(p1 - p2),dims=2) .<= r)] 
            Counter[n, ix] .+= 1
        end
    end

    Phi[1] = log(Logx, N)/N
    Phi[2] = mean(log.(Logx, sum(Counter.>0,dims=2)/N))
    Ap[1]  = Phi[1] - Phi[2]

    for k = 0:m-1
        ai = sum(Counter.>k+1,dims=2)/(N-(k+1)*tau)
        bi = sum(Counter.>k,dims=2)/(N-(k*tau))
        ai = ai[ai.>0]
        bi = bi[bi.>0]
        Phi[k+3] = sum(log.(Logx,ai))/(N-(k+1)*tau)
        Ap[k+2]  = sum(log.(Logx,bi))/(N-(k*tau)) - Phi[k+3]
    end

    return Ap, Phi
    end

end

"""Copyright 2021 Matthew W. Flood, EntropyHub
  
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