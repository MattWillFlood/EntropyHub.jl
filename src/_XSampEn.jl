module _XSampEn
export XSampEn
using Statistics: std, var
using LinearAlgebra: UpperTriangular, I

    """
        XSamp, A, B = XSampEn(Sig1, Sig2) 

    Returns the cross-sample entropy estimates (`XSamp`) and the number of 
    matched vectors (m:B, m+1:A) for m = [0,1,2] estimated for the two 
    univariate data sequences contained in `Sig1` and  `Sig2` using the default parameters:
    embedding dimension = 2, time delay = 1, 
    radius distance threshold= 0.2*SDpooled(`Sig1`,`Sig2`),  logarithm = natural

        XSamp, A, B, (Vcp, Ka, Kb) = XSampEn(Sig1, Sig2, ..., Vcp = true) 
    
    If `Vcp == true`, an additional tuple `(Vcp, Ka, Kb)` is returned with    
    the cross-sample entropy estimates (`XSamp`) and the number of matched state
    vectors (`m: B`, `m+1: A`). `(Vcp, Ka, Kb)` contains the variance of the conditional
    probabilities (`Vcp`), i.e. CP = A/B,  and the number of **overlapping**
    matching vector pairs of lengths m+1 (`Ka`) and m (`Kb`),
    respectively.  Note `Vcp` is undefined for the zeroth embedding dimension (m = 0) 
    and due to the computational demand, **will take substantially more time to return function outputs.**
    See Appendix B in [2] for more info.

        XSamp, A, B = XSampEn(Sig1::Union{AbstractMatrix{T}, AbstractVector{T}} where T<:Real, Sig2::Union{AbstractVector{T} where T<:Real, Nothing} = nothing; m::Int=2, tau::Int=1, r::Union{Real,Nothing}=nothing, Logx::Real=exp(1), Vcp::Bool=false)

    Returns the cross-sample entropy estimates (`XSamp`) for dimensions [0,1,...,m]
    estimated between the data sequences in `Sig1` and `Sig2` using the specified 
    'keyword'  arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer  [default: 2]   \n
    `tau`   - Time Delay, a positive integer         [default: 1]   \n
    `r`     - Radius Distance Threshold, a positive scalar [default: 0.2*SDpooled(`Sig1`,`Sig2`)] \n
    `Logx`  - Logarithm base, a positive scalar      [default: natural] \n

    See also `XFuzzEn`, `XApEn`, `SampEn`, `SampEn2D`, `XMSEn`, `ApEn`
    
    # References:
    	[1] Joshua S Richman and J. Randall Moorman. 
            "Physiological time-series analysis using approximate entropy
            and sample entropy." 
            American Journal of Physiology-Heart and Circulatory Physiology
            (2000)
        
        [2] Douglas E Lake, Joshua S Richman, M.P. Griffin, J. Randall Moorman
            "Sample entropy analysis of neonatal heart rate variability."
            American Journal of Physiology-Regulatory, Integrative and Comparative Physiology
            283, no. 3 (2002): R789-R797.

    """
    function XSampEn(Sig1::Union{AbstractMatrix{T}, AbstractVector{T}} where T<:Real, Sig2::Union{AbstractVector{T} where T<:Real, Nothing} = nothing;
                        m::Int=2, tau::Int=1, r::Union{Nothing,Real}=nothing, Logx::Real=exp(1), Vcp::Bool=false)

    if all(isa.((Sig1,Sig2), AbstractVector))
        N1 = size(Sig1,1);  N2 = size(Sig2,1)  
        S1 = copy(Sig1); S2 = copy(Sig2)
    elseif (minimum(size(Sig1))==2 && (Sig2 isa Nothing)) 
        argmin(size(Sig1)) == 2 ? nothing : Sig1 = Sig1'
        S1 = Sig1[:,1]; S2 = Sig1[:,2];
        N1 = maximum(size(Sig1)); N2 = maximum(size(Sig1));
    else   error("""Sig1 and Sig2 must be 2 separate vectors 
                \t\t\t - OR - 
                Sig1 must be 2-column matrix and Sig2 nothing""")
    end

    r isa Nothing ? r = 0.2*sqrt((var(S1,corrected=false)*(N1-1) + var(S2,corrected=false)*(N2-1))/(N1+N2-1)) : nothing

    (N1>=10 && N2>=10) ? nothing :  error("Sig1/Sig2:   sequences must have >= 10 values")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau>0) ? nothing : error("tau:   must be an integer > 0")
    (r >=0) ? nothing : error("r:     must be a positive value")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")

    Counter = 1*(abs.(S1 .- transpose(S2)) .<= r)
    M = vcat(m*ones(Int,N1-(m*tau)), repeat((m-1):-1:1,inner=tau))
    A = zeros(Int,m+1)
    B = zeros(Int,m+1)
    A[1] = sum(Counter);  B[1] = N1*N2;

    for n = 1:N1 - tau 
        ix = findall(Counter[n, :] .== 1)    
        for k = 1:M[n]            
            ix = ix[ix .+ (k*tau) .<= N2]
            isempty(ix) ? break : nothing
            p1 = repeat(transpose(S1[n:tau:n+(tau*k)]), size(ix,1))  
            p2 = S2[ix .+ transpose(collect(0:tau:(k*tau)))]           
            ix = ix[findall(maximum(abs.(p1 - p2),dims=2) .<= r)]  
            Counter[n, ix] .+= 1
        end
    end

    for k = 1:m
        A[k+1] = sum(Counter.>k)
        B[k+1] = sum(Counter.>=k)
    end

    XSamp = -log.(Logx, A./B)
    #return XSamp, A, B
    
    if Vcp
        T1 = getindex.(findall(Counter.>m),1)
        T2 = getindex.(findall(Counter.>m),2)
        Ka = UpperTriangular(((abs.(T1.-T1').<=m*tau) .+ (abs.(T2.-T2').<=m*tau))[1:end-1,2:end])

        T1 = getindex.(findall(Counter[:,1:end-m*tau].>=m),1)
        T2 = getindex.(findall(Counter[:,1:end-m*tau].>=m),2)
        Kb = UpperTriangular(((abs.(T1.-T1').<=(m-1)*tau) .+ (abs.(T2.-T2').<=(m-1)*tau))[1:end-1,2:end])
                    
        Ka = sum(Ka.>0)
        Kb = sum(Kb.>0)
        CP = A[end]/B[end]
        Vcp = (CP*(1-CP)/B[end]) + (Ka - Kb*(CP^2))/(B[end]^2)

        return XSamp, A, B, (Vcp, Ka, Kb)
    
    else
        return XSamp, A, B 
    end

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