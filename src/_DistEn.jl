module _DistEn
export DistEn
using StatsBase: Histogram, fit, skewness
    """ 
    # Dist, Ppi = DistEn(`Sig`) 

    Returns the distribution entropy estimate (`Dist`) and the
    corresponding distribution probabilities (`Ppi`) estimated from 
    the data sequence (`Sig`) using the default  parameters: 
    embedding dimension = 2, time delay = 1, binning method = 'Sturges',
    logarithm = base 2, normalisation = w.r.t # of histogram bins

    # Dist, Ppi = DistEn(`Sig`, 'keyword' = value, ...)

    Returns the distribution entropy estimate (`Dist`) estimated from 
    the data sequence (`Sig`) using the specified 'keyword' arguments:

    #Arguments:
    `m`     - Embedding Dimension, a positive integer  \n
    `tau`   - Time Delay, a positive integer  \n
    `Bins`  - Histogram bin selection method for distance distribution, 
              one of the following:  \n
              an integer > 1 indicating the number of bins, or one of the 
              following strings {'sturges','sqrt','rice','doanes'}
              [default: 'sturges']
    `Logx`  - Logarithm base, a positive scalar (enter 0 for natural log)   \n
    `Norm`  - Normalisation of DistEn value:  \n
              [false]  no normalisation.
              [true]   normalises w.r.t # of histogram bins - default

    # See also `XDistEn`, `DistEn2D`, `MSEn`, `K2En`

    # References:
        [1] Li, Peng, et al.,
            "Assessing the complexity of short-term heartbeat interval 
            series by distribution entropy." 
            Medical & biological engineering & computing 
            53.1 (2015): 77-87. 
    
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
    function DistEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, 
        Bins::Union{Int,String}="Sturges", Logx::Real=2, Norm::Bool=true)

    Logx == 0  ? Logx = exp(1) : nothing
    
    N = size(Sig)[1]
    (N>10) ? nothing : error("Sig:  must be a numeric vector")
    (m > 0) ? nothing : error("m:   must be an integer > 0")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (Logx>0) ? nothing :  error("Logx:  must be a positive number > 0")
    if typeof(Bins)<:Int
        (Bins>1) ? nothing : 
        error("Bins:  must be an integer > 1 (or name of binning method)")        
    elseif typeof(Bins)<:String
        (lowercase(Bins) in ["sturges","sqrt","rice","doanes"]) ? nothing : 
        error("Bins:    must be one of the following strings
                'sturges', 'sqrt', 'rice', 'doanes' (or an integer >1)")
    end

    Nx = size(Sig)[1] - ((m-1)*tau)
    Zm = zeros(Nx,m)
    for n = 1:m
        Zm[:,n] = Sig[(n-1)*tau + 1:Nx+(n-1)*tau]
    end

    DistMat = zeros(Int(Nx*(Nx-1)/2))
    for k = 1:Nx-1
        Ix = [Int((k-1)*(Nx - k/2)+1), Int(k*(Nx-((k+1)/2)))]
        DistMat[Ix[1]:Ix[2]] = maximum(abs.(transpose(Zm[k,:]) .- Zm[k+1:end,:]),dims=2)
        # DistMat[Ix[1]:Ix[2]] = maximum(abs.(repeat(Zm[k,:],outer=(1,Nx-k)) .- Zm[k+1:end,:]'),dims=2)
    end

    Ny = size(DistMat)[1]
    if eltype(Bins)<:Char
        if lowercase(Bins) == "sturges"
            Bx = ceil(log2(Ny) + 1)
        elseif lowercase(Bins) == "rice"
            Bx = ceil(2*(Ny^(1/3)))
        elseif lowercase(Bins) == "sqrt"
            Bx = ceil(sqrt(Ny))
        elseif lowercase(Bins) == "doanes"
            sigma = sqrt(6*(Ny-2)/((Ny+1)*(Ny+3)))
            Bx = ceil(1+log2(Ny)+log2(1+abs(skewness(DistMat)/sigma)))
        else 
            error("Please enter a valid binning method")
        end
    else    
        Bx = Bins
    end
    By = collect(range(minimum(DistMat),maximum(DistMat),length=Int(Bx+1)))
    By[end] += 1; By[1]-= 1
    Ppi = fit(Histogram, DistMat, By).weights/Ny

    if round(sum(Ppi),digits=6) != 1
        @warn("Potential error estimating probabilities (p = $(sum(Ppi))")
        Ppi = Ppi[Ppi.!=0]
    elseif any(Ppi.==0)
        print("Note: $(sum(Ppi.==0))/$(length(Ppi)) bins were empty \n")
        Ppi = Ppi[Ppi.!=0]
    end
    Dist = -sum(Ppi.*log.(Logx, Ppi))
    if Norm
        Dist = Dist/(log(Logx, Bx));
    end

    return Dist, Ppi
    end

end