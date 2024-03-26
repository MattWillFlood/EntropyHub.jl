module _XK2En
export XK2En
using Statistics: std, mean, var
    """
        XK2, Ci = XK2En(Sig1, Sig2) 

    Returns the cross-Kolmogorov entropy estimates (`XK2`) and the correlation
    integrals (`Ci`) for m = [1,2] estimated between the data sequences 
    contained in `Sig1` and `Sig2` using the default parameters: 
    embedding dimension = 2, time delay = 1, distance threshold (r) = 0.2*SDpooled(Sig1, Sig2),
    logarithm = natural

        XK2, Ci = XK2En(Sig1::Union{AbstractMatrix{T}, AbstractVector{T}} where T<:Real, Sig2::Union{AbstractVector{T} where T<:Real, Nothing} = nothing; m::Int=2, tau::Int=1, r::Union{Real,Nothing}=nothing, Logx::Real=exp(1))

    Returns the cross-Kolmogorov entropy estimates (`XK2`) estimated between
    the data sequences contained in `Sig1` and `Sig2` using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer [default: 2]  \n
    `tau`   - Time Delay, a positive integer         [default: 1]  \n
    `r`     - Radius Distance Threshold, a positive scalar  [default: 0.2*SDpooled(`Sig1`,`Sig2`)]  \n
    `Logx`  - Logarithm base, a positive scalar      [default: natural]  \n

    # See also `XSampEn`, `XFuzzEn`, `XApEn`, `K2En`, `XMSEn`, `XDistEn`

    # References:
        [1]  Matthew W. Flood,
             "XK2En - EntropyHub Project"
             (2021) https://github.com/MattWillFlood/EntropyHub

    """
    function XK2En(Sig1::Union{AbstractMatrix{T}, AbstractVector{T}} where T<:Real, Sig2::Union{AbstractVector{T} where T<:Real, Nothing} = nothing;
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
    (Logx>0) ? nothing : error("Logx: must be a positive number > 0")
    (r>0) ? nothing : error("r:  must be 2 element tuple of positive values")

    m   += 1
    Zm1 = zeros(N1,m)
    Zm2 = zeros(N2,m)
    Ci = zeros(m)
    for n = 1:m
        Nx = N1-(n-1)*tau
        Zm1[1:Nx,n] = S1[(n-1)*tau + 1:N1]
        Ny = N2-(n-1)*tau
        Zm2[1:Ny,n] = S2[(n-1)*tau + 1:N2]   
        Norm = zeros(Nx,Ny)    
        for k = 1:Nx
            Temp = repeat(transpose(Zm1[k,1:n]),outer=Ny) .- Zm2[1:Ny,1:n]
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