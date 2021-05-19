module _XSampEn
export XSampEn
using Statistics: std
    """
    # XSamp, A, B = XSampEn(`Sig`) 

    Returns the cross-sample entropy estimates (`XSamp`) and the number of 
    matched vectors (m:B, m+1:A) for m = [0,1,2] estimated for the two 
    univariate data sequences contained in 'Sig' using the default parameters:
    embedding dimension = 2, time delay = 1, 
    radius distance threshold= 0.2*SD(`Sig`),  logarithm = natural

    # XSamp, A, B = XSampEn(`Sig`, 'keyword' = value, ...)

    Returns the cross-sample entropy estimates (`XSamp`) for dimensions [0,1,...,m]
    estimated between the data sequences in `Sig` using the specified 
    'keyword'  arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer  [default: 2]   \n
    `tau`   - Time Delay, a positive integer         [default: 1]   \n
    `r`     - Radius Distance Threshold, a positive scalar [default: 0.2*SD(`Sig`)] \n
    `Logx`  - Logarithm base, a positive scalar      [default: natural] \n

    See also `XFuzzEn`, `XApEn`, `SampEn`, `SampEn2D`, `XMSEn`, `ApEn`
    
    # References:
    	[1] Joshua S Richman and J. Randall Moorman. 
            "Physiological time-series analysis using approximate entropy
            and sample entropy." 
            American Journal of Physiology-Heart and Circulatory Physiology
            (2000)

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
    function XSampEn(Sig::AbstractArray{T,2} where T<:Real; m::Int=2, 
        tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1))

    (size(Sig,2) > size(Sig,1)) ? Sig = transpose(Sig) : nothing

    N = size(Sig,1)
    (N>10 && size(Sig,2)==2) ? nothing :  error("Sig:   must be a 2-columns matrix")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau>0) ? nothing : error("tau:   must be an integer > 0")
    (r >=0) ? nothing : error("r:     must be a positive value")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")

    S1 = Sig[:,1]; S2 = Sig[:,2]
    Counter = 1*(abs.(S1 .- transpose(S2)) .<= r)
    M = vcat(m*ones(Int,N-(m*tau)), repeat((m-1):-1:1,inner=tau))
    A = zeros(Int,m+1)
    B = zeros(Int,m+1)
    A[1] = sum(Counter);  B[1] = N*N;

    for n = 1:N - tau 
        ix = findall(Counter[n, :] .== 1)    
        for k = 1:M[n]            
            ix = ix[ix .+ (k*tau) .<= N]
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
    return XSamp, A, B
    end

end