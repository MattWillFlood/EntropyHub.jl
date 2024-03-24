module _CoSiEn
export CoSiEn
using Statistics: std, mean, median
using LinearAlgebra: Diagonal, UpperTriangular
    """
        CoSi, Bm = CoSiEn(Sig) 

    Returns the cosine similarity entropy (`CoSi`) and the corresponding
    global probabilities estimated from the data sequence (`Sig`) using the
    default parameters:   embedding dimension = 2, time delay = 1, 
    angular threshold = .1,  logarithm = base 2,

        CoSi, Bm = CoSiEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, r::Real=.1, Logx::Real=2, Norm::Int=0)

    Returns the cosine similarity entropy (`CoSi`) estimated from the data
    sequence (`Sig`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, an integer > 1   \n
    `tau`   - Time Delay, a positive integer    \n
    `r`     - Angular threshold, a value in range [0 < r < 1]   \n
    `Logx`  - Logarithm base, a positive scalar (enter 0 for natural log) \n
    `Norm`  - Normalisation of `Sig`, one of the following integers:    \n
            [0]  no normalisation - default
            [1]  normalises `Sig` by removing median(`Sig`)
            [2]  normalises `Sig` by removing mean(`Sig`)
            [3]  normalises `Sig` w.r.t. SD(`Sig`)
            [4]  normalises `Sig` values to range [-1 1]

    # See also `PhasEn`, `SlopEn`, `GridEn`, `MSEn`, `cMSEn`
    
    # References:
        [1] Theerasak Chanwimalueang and Danilo Mandic,
            "Cosine similarity entropy: Self-correlation-based complexity
            analysis of dynamical systems."
            Entropy 
            19.12 (2017): 652.

 
    """
    function CoSiEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, 
        r::Real=.1, Logx::Real=2, Norm::Int=0)

    Logx == 0  ? Logx = exp(1) : nothing

    N = size(Sig,1)
    (N > 10) ? nothing : error("Sig:   must be a numeric vector")
    (m > 1) ? nothing :  error("m:     must be an integer > 1")
    (tau>0) ? nothing :  error("tau:   must be an integer > 0")
    (0<r<1) ? nothing :  error("r:     must be a scalar in range [0 1]")
    (Logx>0) ? nothing : error("Logx:  must be a positive number > 0")
    (Norm in collect(0:4)) ? nothing : error("Norm:   must be an integer in range [0 4]")

    if Norm == 1
        Xi = Sig .- median(Sig);
    elseif Norm == 2
        Xi = Sig .- mean(Sig);
    elseif Norm == 3
        Xi = (Sig .- mean(Sig))/std(Sig,corrected=false)
    elseif Norm == 4
        Xi = (2*(Sig .- minimum(Sig))/(maximum(Sig)-minimum(Sig))) .- 1;
    else
        Xi = Sig;
    end

    Nx = N-((m-1)*tau);
    Zm = zeros(Nx,m);
    for n = 1:m
        Zm[:,n] = Xi[(n-1)*tau+1:Nx+(n-1)*tau]
    end

    Num = Zm*transpose(Zm); 
    Mag = sqrt.(sum(Diagonal(Num),dims=1))[:]
    Den = Mag*transpose(Mag)
    AngDis = round.(acos.(round.(Num./Den,digits=6))/pi,digits=6)
    if maximum(imag.(AngDis)) < (10^-5)
        Bm = (sum(UpperTriangular(AngDis .< r))-Nx)/(Nx*(Nx-1)/2)
    else
        Bm = (sum(UpperTriangular(real.(AngDis) .< r))-Nx)/(Nx*(Nx-1)/2)
        @warn("Complex values ignored.")
    end
    if Bm == 1 || Bm == 0
        CoSi = 0
    else
        CoSi = -(Bm*log(Logx, Bm)) - ((1-Bm)*log(Logx, 1-Bm))
    end

    return CoSi, Bm
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