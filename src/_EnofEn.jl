module _EnofEn
export EnofEn
using StatsBase: countmap, Histogram, fit
    """
        EoE, AvEn, S2 = EnofEn(Sig) 

    Returns the entropy of entropy (`EoE`), the average Shannon entropy (`AvEn`), and
    the number of levels (`S2`) across all windows estimated from the data sequence (`Sig`) using
    the default parameters: 
    window length (samples) = 10, slices = 10, logarithm = natural, 
    heartbeat interval range (xmin, xmax) = (min(Sig), max(Sig))


        EoE, AvEn, S2 = EnofEn(Sig::AbstractArray{T,1} where T<:Real; tau::Int=10, S::Int=10, Xrange::Tuple{Real,REal}, Logx::Real=exp(1))

    Returns the entropy of entropy (`EoE`) estimated from the data sequence 
    (`Sig`)  using the specified 'keyword' arguments:

    # Arguments:
    `tau`    - Window length, an integer > 1 \n
    `S`      - Number of slices (s1,s2), a two-element tuple of integers > 2 \n
    `Xrange` - The min and max heartbeat interval, a two-element tuple where X[1] <= X[2]\n
    `Logx`   - Logarithm base, a positive scalar  \n

    # See also `SampEn`, `MSEn`, `ApEn`

    # References:
        [1] Chang Francis Hsu, et al.,
            "Entropy of entropy: Measurement of dynamical complexity for
            biological systems." 
            Entropy 
            19.10 (2017): 550.

    """
    function EnofEn(Sig::AbstractArray{T,1} where T<:Real; tau::Int=10, S::Int=10,
                        Xrange::Tuple{Real,Real}=(minimum(Sig),maximum(Sig)), Logx::Real=exp(1))
        
    N = size(Sig,1)
    (N > 10) ? nothing :  error("Sig:   must be a numeric vector")
    (tau > 1 && tau < length(Sig)) ? nothing : error("tau:   must be an integer > 1")
    (S > 1) ? nothing : error("S:  must be an integer > 1")
    (length(Xrange)==2 && (Xrange[1]<=Xrange[2])) ? nothing : error("Xrange:   must be a two-element numeric tuple where Xrange[1]<Xrange[2]")
    (Logx>0) ? nothing :  error("Logx:     must be a positive number > 0")

    Wn = Int(floor(N/tau))
    Wj = transpose(reshape(Sig[1:Wn*tau],tau,Wn))
    Yj = zeros(Wn)
    #Edges = collect(range(minimum(Sig),maximum(Sig),length=(S[1]+1)))
    Edges = collect(range(Xrange[1],Xrange[2],length=(S+1)))
    Edges[1] -= .1;  Edges[end] += .1
    for n = 1:Wn
        Temp = fit(Histogram,Wj[n,:],Edges).weights/tau
        Temp = Temp[Temp.>0]
        Yj[n] = -sum(Temp.*log.(Logx, Temp))
    end

    AvEn = sum(Yj)/Wn
    #Edges = collect(range(minimum(Yj),maximum(Yj),length=(S[2]+1)))
    #Edges[1] -= .1;  Edges[end] += .1
    #Pjl = fit(Histogram,Yj,Edges).weights/Wn
    #Pjl = Pjl[Pjl.>0]
    Pjl = collect(values(countmap(round.(Yj,digits=12))))./Wn    
    S2 = length(Pjl)
    if round(sum(Pjl),digits=5) != 1
        @warn("Possible error estimating probabilities")
    end

    EoE = -sum(Pjl.*log.(Logx, Pjl))
    return EoE, AvEn, S2

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