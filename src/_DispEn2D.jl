module _DispEn2D
export DispEn2D
using Clustering: kmeans, assignments
using Statistics: std, mean
using StatsFuns: normcdf
    """
        Disp2D, RDE = DispEn2D(Mat) 

    Returns the bidimensional dispersion entropy estimate (`Disp2D`) and reverse
    bidimensional dispersion entropy (`RDE`) estimated for the data matrix (`Mat`) 
    using the default parameters:  time delay = 1, symbols = 3, logarithm = natural, 
    data transform = normalised cumulative density function (`'ncdf'`), logarithm = natural, 
    template matrix size = [floor(H/10) floor(W/10)], (where H and W represent
    the height (rows) and width (columns) of the data matrix `Mat`) \n
    ** The minimum number of rows and columns of Mat must be > 10.**

        Disp2D, RDE = DispEn2D(Mat::AbstractArray{T,2} where T<:Real; 
                            m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10), tau::Int=1,
                                c::Int=3, Typex::String="ncdf", Logx::Real=exp(1), Norm::Bool=false, Lock::Bool=true)

    Returns the bidimensional dispersion entropy (`Disp2D`) and reverse 
    bidimensional distribution entropy (`RDE`) estimate for the data matrix (`Mat`) 
    using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Template submatrix dimensions, an integer scalar (i.e. the same
              height and width) or a two-element tuple of integers 
              [height, width] with a minimum value > 1.
              [default: [floor(H/10) floor(W/10)]] \n
    `tau`   - Time Delay, a positive integer       [default: 1]  \n
    `c`     - Number of symbols, an integer > 1
    `Typex` - Type of symbolic mapping transform, one of the following:
              {`linear`, `kmeans`, `ncdf`, `equal`}
              See the `EntropyHub Guide` for more info on these transforms.    
    `Logx`  - Logarithm base, a positive scalar        [default: natural]\n
              ** enter 0 for natural logarithm.**\n
    `Norm`  - Normalisation of `Disp2D` value, a boolean:
                - [false]   no normalisation - default
                - [true]    normalises w.r.t number of possible dispersion patterns.    
    `Lock`  - By default, DispEn2D only permits matrices with a maximum
              size of 128 x 128 to prevent memory errors when storing data on RAM. 
              e.g. For Mat = [200 x 200], m = 3, and tau = 1, DispEn2D 
              creates a vector of 753049836 elements. To enable matrices
              greater than [128 x 128] elements, set `Lock` to false.
              [default: 'true']
              `WARNING: unlocking the permitted matrix size may cause your Julia
               IDE to crash.`

    # See also `DispEn`, `DistEn2D`, `SampEn2D`, `FuzzEn2D`, `MSEn`

    # References:
        [1] Hamed Azami, et al.,
            "Two-dimensional dispersion entropy: An information-theoretic 
            method for irregularity analysis of images."
            Signal Processing: Image Communication, 
            75 (2019): 178-187.

    """
    function DispEn2D(Mat::AbstractArray{T,2} where T<:Real; m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10), tau::Int=1,
                        c::Int=3, Typex::String="ncdf", Logx::Real=exp(1), Norm::Bool=false, Lock::Bool=true)

    Logx == 0  ? Logx = exp(1) : nothing

    NL, NW = size(Mat)
    ((NL > 128 || NW > 128) && Lock) ? 
        error("To prevent memory errors, matrix width & length must have <= 128 elements.
            To estimatDispEn2D  for the current matrix ($NL,$NW) change Lock to 'false'.
            Caution: unlocking the safe matrix size may cause the Julia IDE to crash.") :
            nothing

    length(m)==1 ?  (mL = m; mW = m) : (mL = m[1];  mW = m[2])
    (NL > 10 && NW > 10) ? nothing : 
    error("Number of rows and columns in Mat must be > 10")
    (minimum(m) > 1) ? nothing : error("m:   must be an integer > 1, or 2-element integer tuple w/ values > 1")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (c > 1) ? nothing :   error("c:     must be an integer > 1")
    (lowercase(Typex) in ["linear", "kmeans", "ncdf", "equal"]) ? nothing :
        error("Typex:    must be one of the following strings - 'linear','kmeans','ncdf','equal'")    
    (Logx>0) ? nothing :  error("Logx:  must be a positive number > 0")

    if lowercase(Typex) == "linear"    
        Edges = range(minimum(Mat),maximum(Mat),length=c+1)
        Zi = map(x -> sum(Edges[1:c].<=x), Mat)

    elseif lowercase(Typex) == "kmeans"
        Temp = kmeans(transpose(Mat[:]), c; maxiter=200)
        Zx = assignments(Temp)
        Clux = Temp.centers
        xx = sortperm(Clux[:]);   Zi = zeros(Int,length(Mat))
        for k = 1:c
            Zi[Zx.==xx[k]] .= k;
        end
        Zi = reshape(Zi,size(Mat))
            
    elseif lowercase(Typex) == "ncdf"
        Zx = normcdf.(mean(Mat),std(Mat,corrected=false),Mat);
        Zi = map(x -> sum(range(0,1,length=c+1)[1:c].<=x), Zx);

    elseif lowercase(Typex) == "equal"
        ix = sortperm(Mat[:]);
        xx = Int.(round.(range(0,length(Mat),length=c+1)))
        Zi = zeros(Int, length(Mat))
        for k = 1:c
            Zi[ix[xx[k]+1:xx[k+1]]] .= k
        end
        Zi = reshape(Zi,size(Mat))
    end

    NL = NL - (mL-1)*tau
    NW = NW - (mW-1)*tau
    X = zeros(Int,NL*NW,mL*mW)
    p = 0
    for k = 1:NL
        for n = 1:NW
            p += 1
            X[p,:] = Zi[k:tau:(mL-1)*tau+k,n:tau:(mW-1)*tau+n][:]
        end
    end

    p != NL*NW ? @warn("Potential error with submatrix division.") : nothing

    T = unique(X,dims=1)
    Nx = size(T)[1]
    Counter = zeros(Nx)
    for n = 1:Nx
        Counter[n] = sum(all(X .- transpose(T[n,:]) .==0, dims=2))
    end    
    
    Ppi = Counter[Counter.!= 0]/size(X)[1]
    
    big(c)^(mL*mW) > 10^16 ? error("RDE cannot be estimated with c = $c and 
    a submatrix of size $mL x $mW. Required floating point precision exceeds 10^16.
    Consider reducing the template submatrix size (m) or the number of symbols (c).") : nothing
    
    #RDE = sum(Ppi.^2) - (1/Nx)
    RDE = sum((Ppi .- (1/(c^(mL*mW)))).^2);

    if round(sum(Ppi),digits=4) != 1
        @warn("Potential Error calculating probabilities")
    end

    Disp2D = -sum(Ppi.*log.(Logx, Ppi))
    if Norm
        Disp2D = Disp2D/log(Logx, c^(mL*mW))
        RDE = RDE./(1 - (1/(c^(mL*mW))))
    end

    return Disp2D, RDE
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