module _XApEn
export XApEn
using Statistics: mean, std, var
    """
        XAp, Phi = XApEn(Sig1, Sig2)

    Returns the cross-approximate entropy estimates (`XAp`) and the average
    number of matched vectors (`Phi`) for m = [0,1,2], estimated for the data
    sequences contained in `Sig1` and `Sig2` using the default parameters:
    embedding dimension = 2, time delay = 1, 
    radius distance threshold= 0.2*SDpooled(`Sig1`,`Sig2`), logarithm = natural

    **NOTE**: XApEn is direction-dependent. Thus, `Sig1` is used as the template data sequence,
    and `Sig2` is the matching sequence.``

        XAp, Phi = XApEn(Sig1::Union{AbstractMatrix{T}, AbstractVector{T}} where T<:Real, Sig2::Union{AbstractVector{T} where T<:Real, Nothing} = nothing; m::Int=2, tau::Int=1, r::Union{Real,Nothing}=nothing, Logx::Real=exp(1))
        
    Returns the cross-approximate entropy estimates (`XAp`) between the data
    sequences contained in `Sig1` and `Sig2` using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer  [default: 2] \n
    `tau`   - Time Delay, a positive integer        [default: 1]  \n
    `r`     - Radius Distance Threshold, a positive scalar [default: 0.2*SDpooled(`Sig1`,`Sig2`)] \n
    `Logx`  - Logarithm base, a positive scalar     [default: natural]  \n

    # See also `XSampEn`, `XFuzzEn`, `XMSEn`, `ApEn`, `SampEn`, `MSEn`
  
    # References:
        [1] Steven Pincus and Burton H. Singer,
            "Randomness and degrees of irregularity." 
            Proceedings of the National Academy of Sciences 
            93.5 (1996): 2083-2088.
  
        [2] Steven Pincus,
            "Assessing serial irregularity and its implications for health."
            Annals of the New York Academy of Sciences 
            954.1 (2001): 245-267.
  
    """
    function XApEn(Sig1::Union{AbstractMatrix{T}, AbstractVector{T}} where T<:Real, Sig2::Union{AbstractVector{T} where T<:Real, Nothing} = nothing; 
                        m::Int=2, tau::Int=1, r::Union{Real,Nothing}=nothing, Logx::Real=exp(1))

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
    (Logx>0) ? nothing : error("Logx: must be a positive number > 0")

    Counter = 1*(abs.(S1 .- transpose(S2)) .<= r)
    M = vcat(m*ones(Int,N1-(m*tau)), repeat((m-1):-1:1,inner=tau))
    XAp = zeros(m+1)
    Phi = zeros(m+2)
    for n = 1:N1 - tau 
        ix = findall(Counter[n, :] .== 1)    
        for k = 1:M[n]            
            ix = ix[ix .+ (k*tau) .<= N2]
            isempty(ix) ? break : nothing
            p1 = repeat(transpose(S1[n:tau:n+(tau*k)]), outer=length(ix))  
            p2 = S2[ix .+ transpose(collect(0:tau:(k*tau)))]           
            ix = ix[findall(maximum(abs.(p1 - p2),dims=2) .<= r)]  
            Counter[n, ix] .+= 1
        end
    end

    #Phi[1] = log(Logx, N1)/N1
    #Phi[2] = mean(log.(Logx, sum(Counter.>0,dims=1)/N)) 
    Temp = sum(Counter.>0,dims=1); Temp = Temp[Temp.!=0]
    Phi[2] = mean(log.(Logx, Temp/N1)) 
    XAp[1] = Phi[1] - Phi[2]
    for k = 0:m-1
        ai = sum(Counter.>(k+1),dims=1)/(N1-(k+1)*tau)
        bi = sum(Counter.>k,dims=1)/(N1-(k*tau))
        ai = ai[ai.!=0]
        bi = bi[bi.!=0]
        Phi[k+3] = sum(log.(Logx, ai))/(N1-(k+1)*tau)
        XAp[k+2]= sum(log.(Logx, bi))/(N1-(k*tau)) - Phi[k+3]
    end

    return XAp, Phi
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