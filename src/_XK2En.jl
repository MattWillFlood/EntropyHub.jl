module _XK2En
export XK2En
using Statistics: std, mean
    """
        XK2, Ci = XK2En(Sig) 

    Returns the cross-Kolmogorov entropy estimates (`XK2`) and the correlation
    integrals (`Ci`) for m = [1,2] estimated between the data sequences 
    contained in `Sig` using the default parameters: 
    embedding dimension = 2, time delay = 1, distance threshold (r) = 0.2*SD(Sig),
    logarithm = natural

        XK2, Ci = XK2En(Sig::AbstractArray{T,2} where T<:Real; m::Int=2, tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1))

    Returns the cross-Kolmogorov entropy estimates (`XK2`) estimated between
    the data sequences contained in `Sig` using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer [default: 2]  \n
    `tau`   - Time Delay, a positive integer         [default: 1]  \n
    `r`     - Radius Distance Threshold, a positive scalar  [default: 0.2*SD(`Sig`)]  \n
    `Logx`  - Logarithm base, a positive scalar      [default: natural]  \n

    # See also `XSampEn`, `XFuzzEn`, `XApEn`, `K2En`, `XMSEn`, `XDistEn`

    # References:
        [1]  Matthew W. Flood,
             "XK2En - EntropyHub Project"
             (2021) https://github.com/MattWillFlood/EntropyHub

    """
    function XK2En(Sig::AbstractArray{T,2} where T<:Real; m::Int=2, 
        tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1))

    size(Sig,2) > size(Sig,1) ? Sig = transpose(Sig) : nothing

    N = size(Sig,1)
    (N>10 && size(Sig,2)==2) ? nothing :  error("Sig:   must be a 2-columns matrix")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau>0) ? nothing : error("tau:   must be an integer > 0")
    (Logx>0) ? nothing : error("Logx: must be a positive number > 0")
    (r>0) ? nothing : error("r:  must be 2 element tuple of positive values")

    m   += 1
    Zm1 = zeros(N,m)
    Zm2 = zeros(N,m)
    Ci = zeros(m)
    for n = 1:m
        N2 = N-(n-1)*tau
        Zm1[1:N2,n] = Sig[(n-1)*tau + 1:N,1]
        Zm2[1:N2,n] = Sig[(n-1)*tau + 1:N,2]   
        Norm = zeros(N2,N2)    
        for k = 1:N2
            Temp = repeat(transpose(Zm1[k,1:n]),outer=N2) .- Zm2[1:N2,1:n]
            Norm[k,:] = sqrt.(sum(Temp.*Temp,dims=2))
        end
        Ci[n] = mean(Norm .< r)
    end
    
    XK2 = log.(Logx, Ci[1:m-1]./Ci[2:m])/tau
    XK2[isinf.(XK2)] .= NaN

    return XK2, Ci

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