module _DistEn2D
export DistEn2D
using StatsBase: fit, Histogram, skewness
    """
        Dist2D = DistEn2D(Mat) 

    Returns the bidimensional distribution entropy estimate (`Dist2D`)
    estimated for the data matrix (`Mat`) using the default parameters:
    time delay = 1, histogram binning method = "sturges", logarithm = natural, 
    template matrix size = [floor(H/10) floor(W/10)], (where H and W represent
    the height (rows) and width (columns) of the data matrix `Mat`) \n
    ** The minimum number of rows and columns of Mat must be > 10.**

        Dist2D = DistEn2D(Mat::AbstractArray{T,2} where T<:Real; m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10), tau::Int=1,
                            Bins::Union{Int,String}="Sturges", Logx::Real=2, Norm::Int=2, Lock::Bool=true)

    Returns the bidimensional distribution entropy (`Dist2D`) estimate for 
    the data matrix (`Mat`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Template submatrix dimensions, an integer scalar (i.e. the same
              height and width) or a two-element tuple of integers 
              [height, width] with a minimum value > 1.
              [default: [floor(H/10) floor(W/10)]] \n
    `tau`   - Time Delay, a positive integer       [default: 1]  \n
    `Bins`  - Histogram bin selection method for distance distribution,
              an integer > 1 indicating the number of bins, or one of the 
              following strings {`"sturges", "sqrt", "rice", "doanes"``}
              [default: 'sturges']  \n
    `Logx`  - Logarithm base, a positive scalar        [default: natural]\n
              ** enter 0 for natural logarithm.**\n
    `Norm`  - Normalisation of `Dist2D` value, one of the following integers:
              [0]  no normalisation.
              [1]  normalises values of data matrix (`Mat`) to range [0 1].
              [2]  normalises values of data matrix (`Mat`) to range [0 1],
                   and normalises the distribution entropy value (`Dist2D`)
                   w.r.t the number of histogram bins.  [default]
              [3]  normalises the distribution entropy value
                   w.r.t the number of histogram bins, without normalising
                   data matrix values. \n
    `Lock`  - By default, DistEn2D only permits matrices with a maximum
              size of 128 x 128 to prevent memory errors when storing data on RAM. 
              e.g. For Mat = [200 x 200], m = 3, and tau = 1, DistEn2D 
              creates a vector of 753049836 elements. To enable matrices
              greater than [128 x 128] elements, set `Lock` to false.
              [default: 'true']
              `WARNING: unlocking the permitted matrix size may cause your Julia
               IDE to crash.`

    # See also `DistEn`, `XDistEn`, `SampEn2D`, `FuzzEn2D`, `MSEn`

    # References:
        [1] Hamed Azami, Javier Escudero and Anne Humeau-Heurtier,
            "Bidimensional distribution entropy to analyze the irregularity
            of small-sized textures."
            IEEE Signal Processing Letters 
            24.9 (2017): 1338-1342.


    """
    function DistEn2D(Mat::AbstractArray{T,2} where T<:Real; 
        m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10), tau::Int=1,
        Bins::Union{Int,String}="Sturges", Logx::Real=2, Norm::Int=2, Lock::Bool=true)

    Logx == 0  ? Logx = exp(1) : nothing

    NL, NW = size(Mat)
    ((NL > 128 || NW > 128) && Lock) ? 
        error("To prevent memory errors, matrix width & length must have <= 128 elements.
            To estimate DistEn2D  for the current matrix ($NL,$NW) change Lock to 'false'.
            Caution: unlocking the safe matrix size may cause the Julia IDE to crash.") :
            nothing

    length(m)==1 ?  (mL = m; mW = m) : (mL = m[1];  mW = m[2])
    (NL > 10 && NW > 10) ? nothing : 
    error("Number of rows and columns in Mat must be > 10")
    (minimum(m) > 1) ? nothing : error("m:   must be an integer > 1")
    (tau >0) ? nothing : error("tau:   must be an integer > 0")
    (Logx>0) ? nothing : error("Logx:  must be a positive number > 0")
    if typeof(Bins)<:Int
        (Bins>1) ? nothing : 
        error("Bins:  must be an integer > 1 (or name of binning method)")        
    elseif typeof(Bins)<:String
        (lowercase(Bins) in ["sturges","sqrt","rice","doanes"]) ? nothing : 
        error("Bins:    must be one of the following strings
                'sturges', 'sqrt', 'rice', 'doanes' (or an integer >1)")
    end
    (Norm in [0, 1, 2, 3]) ? nothing :  error("Norm:  must be an integer in the range [0 3]")

    Norm in [1, 2] ? Mat = (Mat.-minimum(Mat))./maximum(Mat.-minimum(Mat)) : nothing
    NL = NL - (mL-1)*tau
    NW = NW - (mW-1)*tau
    X = zeros(mL,mW,NL*NW)
    p = 0
    for k = 1:NL
        for n = 1:NW
            p += 1
            X[:,:,p] = Mat[k:tau:(mL-1)*tau+k,n:tau:(mW-1)*tau+n]
        end
    end

    p = size(X,3)
    p != NL*NW ? @warn("Potential error with submatrix division.") : nothing
    Ny = Int(p*(p-1)/2)
    Ny > 300000000 ? @warn("Number of pairwise distance calculations is $Ny") : nothing

    Y = zeros(Ny)
    for k = 1:p-1
        Ix = Int.([(k-1)*(p - k/2)+1, k*(p-((k+1)/2))])
        Y[Ix[1]:Ix[2]] = maximum(abs.(X[:,:,k+1:end] .- X[:,:,k]),dims=(1,2))
    end
    
    if eltype(Bins)<:Char
        if lowercase(Bins) == "sturges"
            Bx = ceil(log2(Ny) + 1)
        elseif lowercase(Bins) == "rice"
            Bx = ceil(2*(Ny^(1/3)))
        elseif lowercase(Bins) == "sqrt"
            Bx = ceil(sqrt(Ny))
        elseif lowercase(Bins) == "doanes"
            sigma = sqrt(6*(Ny-2)/((Ny+1)*(Ny+3)))
            Bx = ceil(1+log2(Ny)+log2(1+abs(skewness(convert(Array{Float64,1},Y))/sigma)))
        else 
            error("Please enter a valid binning method")
        end
    else    
        Bx = Bins
    end

    By = collect(range(minimum(Y),maximum(Y),length=Int(Bx+1)))
    By[end] += 1; By[1] -= 1

    Ppi = fit(Histogram, Y[:], By).weights/Ny
    if round(sum(Ppi),digits=6) != 1
        @warn("Potential error estimating probabilities (p = $(sum(Ppi)))")
        Ppi = Ppi[Ppi.!=0]
    elseif any(Ppi.==0)
        print("Note: $(sum(Ppi.==0))/$(length(Ppi)) bins were empty \n")
        Ppi = Ppi[Ppi.!=0]
    end
    Dist2D = -sum(Ppi.*log.(Logx, Ppi))
    Norm >= 2  ? Dist2D /= (log(Logx, Bx)) : nothing

    return Dist2D
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



#= Y = []
for k = 1:p-1
    append!(Y,maximum(abs.(X[:,:,k+1:end] .- X[:,:,k]),dims=(1,2)))
end =#