module _K2En
export K2En
using Statistics: std
    """
        K2, Ci = K2En(Sig) 
        
    Returns the Kolmogorov entropy estimates `K2` and the correlation
    integrals `Ci` for `m` = [1,2] estimated from the data sequence `Sig`
    using the default parameters: embedding dimension = 2, time delay = 1, 
    r = 0.2*SD(`Sig`), logarithm = natural

        K2, Ci = K2En(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1))
        
    Returns the Kolmogorov entropy estimates `K2` for dimensions = [1,...,`m`]
    estimated from the data sequence `Sig` using the 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer\n
    `tau`   - Time Delay, a positive integer\n
    `r`     - Radius, a positive scalar \n 
    `Logx`  - Logarithm base, a positive scalar\n  

    # See also `DistEn`, `XK2En`, `MSEn`

    # References:
        [1] Peter Grassberger and Itamar Procaccia,
            "Estimation of the Kolmogorov entropy from a chaotic signal." 
            Physical review A 28.4 (1983): 2591.

        [2] Lin Gao, Jue Wang  and Longwei Chen
            "Event-related desynchronization and synchronization 
            quantification in motor-related EEG by Kolmogorov entropy"
            J Neural Eng. 2013 Jun;10(3):03602

    """
    function K2En(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, r::Real=0.2*std(Sig,corrected=false), Logx::Real=exp(1))

    N = size(Sig)[1]
    (N>10) ? nothing : error("Sig:   must be a numeric vector")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (r>0) ? nothing :  error("r:     must be a positive scalar > 0")
    (Logx>0) ? nothing :  error("Logx:     must be a positive number > 0")
        
    m = m+1
    Zm = zeros(N,m)
    Ci = zeros(m)
    for n = 1:m
        N2 = N-(n-1)*tau
        Zm[1:N2,n] = Sig[(n-1)*tau + 1:N]   
        Norm = Inf*ones(N2-1,N2-1)
        for k = 1:N2-1
            Temp = repeat(transpose(Zm[k,1:n]),N2-k,1) - Zm[k+1:N2,1:n]
            Norm[k,k:N2-1] = sqrt.(sum(Temp.^2, dims=2)) 
        end
        Ci[n] = 2*sum(Norm .< r)/(N2*(N2-1))    
    end
    
    K2 = log.(Logx, Ci[1:m-1]./Ci[2:m])/tau
    K2[isinf.(K2)] .= NaN

    return K2, Ci
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