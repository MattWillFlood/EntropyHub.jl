module _DivEn
export DivEn
using StatsBase: Histogram, fit

    """
        Div, CDEn, Bm = DivEn(Sig) 

    Returns the diversity entropy (`Div`), the cumulative diversity entropy (`CDEn`),
    and the corresponding probabilities (`Bm`) estimated from the data sequence (`Sig`) 
    using the default parameters:   embedding dimension = 2, time delay = 1, #bins = 5,  logarithm = natural,
    
        Div, CDEn, Bm = DivEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, r::Int=5, Logx::Real=exp(1))
    
    Returns the diversity entropy (`Div`) estimated from the data sequence (`Sig`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, an integer > 1   \n
    `tau`   - Time Delay, a positive integer    \n
    `r`     - Histogram bins #: either \n
                * an integer [1 < `r`] representing the number of bins
                * a list/numpy array of 3 or more increasing values in range [-1 1] representing the bin edges including the rightmost edge.\n
    `Logx`  - Logarithm base, a positive scalar  (Enter 0 for natural logarithm)

    # See also  `CoSiEn`, `PhasEn`, `SlopEn`, `GridEn`, `MSEn`
    
    # References:     
        [1] X. Wang, S. Si and Y. Li, 
            "Multiscale Diversity Entropy: A Novel Dynamical Measure for Fault 
            Diagnosis of Rotating Machinery," 
            IEEE Transactions on Industrial Informatics,
            vol. 17, no. 8, pp. 5419-5429, Aug. 2021
            
        [2] Y. Wang, M. Liu, Y. Guo, F. Shu, C. Chen and W. Chen, 
            "Cumulative Diversity Pattern Entropy (CDEn): A High-Performance, 
            Almost-Parameter-Free Complexity Estimator for Nonstationary Time Series,"
            IEEE Transactions on Industrial Informatics
            vol. 19, no. 9, pp. 9642-9653, Sept. 2023      
     
    """
    function DivEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, 
        r::Union{Int, Vector, StepRangeLen, Tuple}=5, Logx::Real=exp(1))

    Logx == 0  ? Logx = exp(1) : nothing
    N = size(Sig,1)
    (N > 10) ? nothing : error("Sig:   must be a numeric vector")
    (m > 1) ? nothing :  error("m:     must be an integer > 1")
    (tau>0) ? nothing :  error("tau:   must be an integer > 0")
    (Logx>0) ? nothing : error("Logx:  must be a positive number > 0")

    r isa Union{Vector, StepRangeLen, Tuple} ? r = collect(r) : nothing
    if r isa Int
        (r<1) ? error("r:    must be an int > 1 or a vector of 3 or more increasing values in range [-1 1]") : r = collect(LinRange(-1,1,r+1))
    elseif r isa Vector
        (r isa Vector && length(r) > 2 && minimum(r) >= -1 && maximum(r)<= 1 && minimum(diff(r))>0) ? nothing : 
        error("r:    must be an int > 1 or a vector of 3 or more increasing values in range [-1 1]") 
    else
        error("r:    must be an int > 1 or a vector of 3 or more increasing values in range [-1 1]") 
    end

    Nx = N - (m-1)*tau
    Zm = zeros((Nx,m))
    for n = 1:m
        Zm[:,n] = Sig[1 + (n-1)*tau:Nx+((n-1)*tau)]
    end

    Num = sum(Zm[1:end-1,:].*Zm[2:end,:],dims=2)
    Den = sqrt.(sum(Zm[2:end,:].^2,dims=2)).*sqrt.(sum(Zm[1:end-1,:].^2,dims=2))
    Di = (Num./Den)[:]
    Bm = fit(Histogram, Di, r).weights
    Bm = Bm[Bm.>0]/sum(Bm)
    
    round(sum(Bm),digits = 6) != 1.0 ? (@warn "Warning: Potential error is probability estimation! Sum(Pi) == " round(sum(Bm),digits=6)) : nothing
    
    r = length(r)-1
    Pj = 1 .- cumsum(Bm)
    Pj = (Pj./sum(Pj))[1:end-1]
    CDEn = -sum(Pj.*log.(Pj)./log(Logx))./(log(r)/log(Logx))
    Div = -sum(Bm.*log.(Bm)./log(Logx))./(log(r)/log(Logx))

    return Div, CDEn, Bm                 
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