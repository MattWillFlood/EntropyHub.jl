module _XApEn
export XApEn
using Statistics: mean, std
    """
        XAp, Phi = XApEn(Sig)

    Returns the cross-approximate entropy estimates (`XAp`) and the average
    number of matched vectors (`Phi`) for m = [0,1,2], estimated for the data
    sequences contained in 'Sig' using the default parameters:
    embedding dimension = 2, time delay = 1, 
    radius distance threshold= 0.2*SD(`Sig`), logarithm = natural

    **NOTE**: XApEn is direction-dependent. Thus, the first column of
    `Sig` is used as the template data sequence, and the second
    column is the matching sequence.``

        XAp, Phi = XApEn(Sig::AbstractArray{T,2} where T<:Real; m::Int=2, tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1))

    Returns the cross-approximate entropy estimates (`XAp`) between the data
    sequences contained in 'Sig' using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer  [default: 2] \n
    `tau`   - Time Delay, a positive integer        [default: 1]  \n
    `r`     - Radius Distance Threshold, a positive scalar [default: 0.2*SD(`Sig`)] \n
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
    function XApEn(Sig::AbstractArray{T,2} where T<:Real; m::Int=2, 
        tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1))

    size(Sig,2) > size(Sig,1) ? Sig = transpose(Sig) : nothing

    N = size(Sig,1)
    (N>10 && size(Sig,2)==2) ? nothing :  error("Sig:   must be a 2-columns matrix")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau>0) ? nothing : error("tau:   must be an integer > 0")
    (r >=0) ? nothing : error("r:     must be a positive value")
    (Logx>0) ? nothing : error("Logx: must be a positive number > 0")

    S1 = Sig[:,1]; S2 = Sig[:,2]
    Counter = 1*(abs.(S1 .- transpose(S2)) .<= r)
    M = vcat(m*ones(Int,N-(m*tau)), repeat((m-1):-1:1,inner=tau))
    XAp = zeros(m+1)
    Phi = zeros(m+2)
    for n = 1:N - tau 
        ix = findall(Counter[n, :] .== 1)    
        for k = 1:M[n]            
            ix = ix[ix .+ (k*tau) .<= N]
            isempty(ix) ? break : nothing
            p1 = repeat(transpose(S1[n:tau:n+(tau*k)]), outer=length(ix))  
            p2 = S2[ix .+ transpose(collect(0:tau:(k*tau)))]           
            ix = ix[findall(maximum(abs.(p1 - p2),dims=2) .<= r)]  
            Counter[n, ix] .+= 1
        end
    end

    Phi[1] = log(Logx, N)/N
    #Phi[2] = mean(log.(Logx, sum(Counter.>0,dims=1)/N)) 
    Temp = sum(Counter.>0,dims=1); Temp = Temp[Temp.!=0]
    Phi[2] = mean(log.(Logx, Temp/N)) 
    XAp[1] = Phi[1] - Phi[2]
    for k = 0:m-1
        ai = sum(Counter.>(k+1),dims=1)/(N-(k+1)*tau)
        bi = sum(Counter.>k,dims=1)/(N-(k*tau))
        ai = ai[ai.!=0]
        bi = bi[bi.!=0]
        Phi[k+3] = sum(log.(Logx, ai))/(N-(k+1)*tau)
        XAp[k+2]= sum(log.(Logx, bi))/(N-(k*tau)) - Phi[k+3]
    end

    return XAp, Phi
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