module _XDistEn
export XDistEn
using StatsBase: fit, Histogram, skewness
using Statistics: mean, std
    """
    # XDist, Ppi = XDistEn(`Sig`) 

    Returns the cross-distribution entropy estimate (`XDist`) and the
    corresponding distribution probabilities (`Ppi`) estimated between the data 
    sequences contained in `Sig` using the default parameters: 
    embedding dimension = 2, time delay = 1, binning method = 'Sturges',
    logarithm = base 2, normalisation = w.r.t # of histogram bins

    # XDist, Ppi = XDistEn(`Sig`, 'keyword' = value, ...)

    Returns the cross-distribution entropy estimate (`XDist`) estimated between the 
    data sequences contained in `Sig` using the specified 'keyword' = arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer   [default: 2]    \n
    `tau`   - Time Delay, a positive integer            [default: 1]    \n
    `Bins`  - Histogram bin selection method for distance distribution, 
              an integer > 1 indicating the number of bins, or one of the 
              following strings {'sturges','sqrt','rice','doanes'}
              [default: 'sturges'] \n
    `Logx`  - Logarithm base, a positive scalar         [default: 2]    
            `` * enter 0 for natural log``\n
    `Norm`  - Normalisation of DistEn value:
              [false]  no normalisation.
              [true]   normalises w.r.t # of histogram bins [default] \n

    # See also `XSampEn`, `XApEn`, `XPermEn`, `XCondEn`, `DistEn`, `DistEn2D`, `XMSEn`
  
    # References:
      [1] Yuanyuan Wang and Pengjian Shang,
            "Analysis of financial stock markets through the multiscale
            cross-distribution entropy based on the Tsallis entropy."
            Nonlinear Dynamics 
            94.2 (2018): 1361-1376.
  
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
    function XDistEn(Sig::AbstractArray{T,2} where T<:Real; m::Int=2, tau::Int=1, 
        Bins::Union{Int,String}="Sturges", Logx::Real=2, Norm::Bool=true)

    (size(Sig,2) > size(Sig,1)) ? Sig = transpose(Sig) : nothing
    Logx == 0  ? Logx = exp(1) : nothing

    N = size(Sig,1)
    (N>10 && size(Sig,2)==2) ? nothing :  error("Sig:   must be a 2-columns matrix")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau>0) ? nothing : error("tau:   must be an integer > 0")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")
    if typeof(Bins)<:Int
        (Bins>1) ? nothing : error("Bins:  must be an integer > 1 (or name of binning method)")        
    elseif typeof(Bins)<:String
        (lowercase(Bins) in ["sturges","sqrt","rice","doanes"]) ? nothing : 
        error("Bins:    must be one of the following strings
                'sturges', 'sqrt', 'rice', 'doanes' (or an integer >1)")
    end

    S1 = Sig[:,1]; S2 = Sig[:,2];
    Nx = size(S1,1) - ((m-1)*tau)
    Zm1 = zeros(Nx,m)
    Zm2 = zeros(Nx,m)
    for n = 1:m
        Zm1[:,n] = S1[(n-1)*tau + 1:Nx+(n-1)*tau]
        Zm2[:,n] = S2[(n-1)*tau + 1:Nx+(n-1)*tau]
    end

    DistMat = zeros(Nx,Nx);
    for k = 1:Nx
        DistMat[k,:] = maximum(abs.(repeat(transpose(Zm1[k,:]),outer=Nx) - Zm2),dims=2)
    end

    Ny = Nx*Nx
    DistMat = reshape(DistMat,1,Ny)
    if eltype(Bins)<:Char
        if lowercase(Bins) == "sturges"
            Bx = ceil(log2(Ny) + 1)
        elseif lowercase(Bins) == "rice"
            Bx = ceil(2*(Ny^(1/3)));
        elseif lowercase(Bins) == "sqrt"
            Bx = ceil(sqrt(Ny));
        elseif lowercase(Bins) == "doanes"
            sigma = sqrt(6*(Ny-2)/((Ny+1)*(Ny+3)));
            Bx = ceil(1+log2(Ny)+log2(1+abs(skewness(DistMat)/sigma)));
        else 
            error("Please enter a valid binning method")
        end
    else    
        Bx = Bins
    end

    By = collect(range(minimum(DistMat),maximum(DistMat),length=Int(Bx+1)))
    By[end] += 1; By[1] -= 1
#    Ppi = fit(Histogram, transpose(DistMat), By).weights/Ny
    Ppi = fit(Histogram, DistMat[:], By).weights/Ny
 
    if round(sum(Ppi),digits=6) != 1
        @warn("Potential error estimating probabilities (p = $(sum(Ppi))")
        Ppi = Ppi[Ppi.!=0]
    elseif any(Ppi.==0)
        print("Note: $(sum(Ppi.==0))/$(length(Ppi)) bins were empty")
        Ppi = Ppi[Ppi.!=0]
    end
    XDist = -sum(Ppi.*log.(Logx, Ppi))
    Norm ? XDist = XDist/log(Logx, Bx) : nothing

    return XDist, Ppi
    end
end