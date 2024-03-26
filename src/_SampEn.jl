module _SampEn
export SampEn
using Statistics: mean, std
using LinearAlgebra: UpperTriangular, I
"""
    Samp, A, B = SampEn(Sig) 

 Returns the sample entropy estimates `Samp` and the number of matched state 
 vectors (`m`:B, `m+1`:A) for `m` = [0,1,2] estimated from the data sequence `Sig`
 using the default parameters: embedding dimension = 2, time delay = 1, 
 radius threshold = 0.2*SD(`Sig`), logarithm = natural
 
    Samp, A, B, (Vcp, Ka, Kb) = SampEn(Sig, ..., Vcp = true) 
    
 If `Vcp == true`, an additional tuple `(Vcp, Ka, Kb)` is returned with    
 the sample entropy estimates (`Samp`) and the number of matched state
 vectors (`m: B`, `m+1: A`). `(Vcp, Ka, Kb)`  contains the variance of the conditional
 probabilities (`Vcp`), i.e. CP = A/B,  and the number of **overlapping**
 matching vector pairs of lengths m+1 (`Ka`) and m (`Kb`),
 respectively.  Note `Vcp` is undefined for the zeroth embedding dimension (m = 0) 
 and due to the computational demand, **will take substantially more time to return function outputs.**
 See Appendix B in [2] for more info.

    Samp, A, B = SampEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1), Vcp::Bool=false)
    
 Returns the sample entropy estimates `Samp` for dimensions = [0,1,...,`m`]
 estimated from the data sequence `Sig` using the specified keyword arguments:
    
# Arguments:
 `m`     - Embedding Dimension, a positive integer\n
 `tau`   - Time Delay, a positive integer\n
 `r`     - Radius Distance Threshold, a positive scalar  \n
 `Logx`  - Logarithm base, a positive scalar \n 
 
# See also `ApEn`, `FuzzEn`, `PermEn`, `CondEn`, `XSampEn`, `SampEn2D`, `MSEn`
  
# References:
    [1] Joshua S Richman and J. Randall Moorman. 
        "Physiological time-series analysis using approximate entropy
        and sample entropy." 
        American Journal of Physiology-Heart and Circulatory Physiology (2000).

    [2] Douglas E Lake, Joshua S Richman, M.P. Griffin, J. Randall Moorman
        "Sample entropy analysis of neonatal heart rate variability."
        American Journal of Physiology-Regulatory, Integrative and Comparative Physiology
        283, no. 3 (2002): R789-R797.

"""
    function SampEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1), Vcp::Bool=false)

    N = length(Sig)
    (N>10) ? nothing : error("Sig:   must be a numeric vector")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (r>=0) ? nothing : error("r:     must be a positive scalar value")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")

    Counter = 1*(abs.(Sig .- transpose(Sig)) .<= r).*UpperTriangular(ones(N,N)) - I(N)
    M = Int.([m*ones(N-m*tau); repeat(collect(m-1:-1:1),inner=tau)])
    A = zeros(m + 1)
    B = zeros(m + 1)
    A[1] = sum(Counter)
    B[1] = N*(N-1)/2

    for n = 1:N-tau
        ix = findall(Counter[n, :] .== 1)    
        for k = 1:M[n]            
            ix = ix[ix .+ (k*tau) .<= N]
            p1 = repeat(transpose(Sig[n:tau:n+(tau*k)]), length(ix))  
            p2 = Sig[ix .+ transpose(collect(0:tau:(k*tau)))]           
            ix = ix[findall(maximum(abs.(p1 - p2),dims=2) .<= r)]  
            if length(ix)>0
                Counter[n, ix] .+= 1
            else
                break
            end      
        end
    end

    for k = 1:m
        A[k+1] = sum(Counter.>k)
        B[k+1] = sum(Counter[:,1:N-(k*tau)].>=k)
    end

    Samp = -log.(Logx, A./B)

    if Vcp
        Temp = hcat(getindex.(findall(Counter.>m),1), getindex.(findall(Counter.>m),2))
        if length(Temp[:,1])>1
            Ka = zeros(Int, length(Temp[:,1]) -1)             
            for k = 1:size(Temp,1)-1 # (length(Temp[:,1])-1)
                TF = (abs.(Temp[k+1:end,:] .- Temp[k,1]) .<= m*tau) .+ (abs.(Temp[k+1:end,:] .- Temp[k,2]) .<= m*tau)
                Ka[k] = sum(any(TF.>0, dims=2))
            end
        else
            Ka = 0
        end
        
        Temp = hcat(getindex.(findall(Counter[:,1:end-m*tau].>=m),1), getindex.(findall(Counter[:,1:end-m*tau].>=m),2))       
        if length(Temp[:,1])>1
            Kb = zeros(Int, length(Temp[:,1]) -1)   
            for k = 1:size(Temp,1)-1  # (length(Temp[:,1]) -1)
                TF = (abs.(Temp[k+1:end,:] .- Temp[k,1]) .<= (m-1)*tau) + (abs.(Temp[k+1:end,:] .- Temp[k,2]) .<= (m-1)*tau)
                Kb[k] = sum(any(TF.>0, dims=2))      
            end
        else
            Kb = 0
        end

        Ka = sum(Ka)
        Kb = sum(Kb)
        CP = A[end]/B[end]
        Vcp = (CP*(1-CP)/B[end]) + (Ka - Kb*(CP^2))/(B[end]^2)

        return Samp, A, B, (Vcp, Ka, Kb)
    
    else
        return Samp, A, B 

    # return Samp, A, B
    end

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