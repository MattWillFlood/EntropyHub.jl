module _XPermEn
export XPermEn
using StatsBase: fit, Histogram
using Combinatorics: permutations
    """
    # XPerm = XPermEn(`Sig`) 

    Returns the cross-permuation entropy estimates (`XPerm`) estimated betweeen
    the data sequences contained in `Sig` using the default parameters:
    embedding dimension = 3, time delay = 1, logarithm = base 2, 

    # XPerm = XPermEn(`Sig`, 'keyword' = value, ...)

    Returns the permutation entropy estimates (`XPerm`) estimated between the 
    data sequences contained in `Sig` using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, an integer > 2   [default: 3]    
            ``**Note: XPerm is undefined for embedding dimensions < 3.``\n
    `tau`   - Time Delay, a positive integer        [default: 1]    \n
    `Logx`  - Logarithm base, a positive scalar     [default: 2]    
            ``* enter 0 for natural log.``    \n

    # See also `PermEn`, `XApEn`, `XSampEn`, `XFuzzEn`, `XMSEn` 

    # References:
        [1] Wenbin Shi, Pengjian Shang, and Aijing Lin,
            "The coupling analysis of stock market indices based on 
            cross-permutation entropy."
            Nonlinear Dynamics
            79.4 (2015): 2439-2447.

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
    function XPermEn(Sig::AbstractArray{T,2} where T<:Real; m::Int=3, tau::Int=1, Logx::Real=exp(1))

    (size(Sig,2) > size(Sig,1)) ? Sig = transpose(Sig) : nothing

    N = size(Sig,1)
    (N>10 && size(Sig,2)==2) ? nothing :  error("Sig:   must be a 2-columns matrix")
    (m > 2) ? nothing : error("m:     must be an integer > 1")
    (tau>0) ? nothing : error("tau:   must be an integer > 0")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")

    S1 = Sig[:,1]; S2 = Sig[:,2];
    N = length(S1)-(m-1)*tau
    Sx1 = zeros(N,m)
    Sx2 = zeros(N,m)
    for k = 1:m    
        Sx1[:,k] = S1[1+(k-1)*tau:N+(k-1)*tau]
        Sx2[:,k] = S2[1+(k-1)*tau:N+(k-1)*tau]
    end

    Temp = sortind(Sx1[1:N,1:m])
    Gx = zeros(N,m)
    for k = 1:N
        Gx[k,:] = Sx2[k,Temp[k,:]]   
    end

    Kt = zeros(m-2,m-2,N)
    for k = 1:m-2           
        for j = k+1:m-1           
            G1 = Gx[:,j+1] .- Gx[:,k]
            G2 = Gx[:,k] .- Gx[:,j]       
            Kt[k,j-1,:] = (G1.*G2 .> 0)       
        end      
    end

    Di = sum(Kt,dims=(1,2))[:]
    Ppi = fit(Histogram, Di, -.5:((m-2)*(m-1) + 1)/2).weights/N
    Ppi = Ppi[Ppi.!=0]
    XPerm = -sum(Ppi.*log.(Logx,Ppi))
    if round(sum(Ppi),digits=6)!=1
        @warn("Potential error with probability calculation")
    end

    return XPerm
    end

    function sortind(X)
        Y = zeros(Int, size(X))
        for k = 1:length(X[:,1])
            Y[k,:] = sortperm(X[k,:])
        end
        return Y
    end
end