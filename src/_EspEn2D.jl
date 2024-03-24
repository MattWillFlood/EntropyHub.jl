module _EspEn2D   
export EspEn2D
    """
        Esp2D,  = EspEn2D(Mat) 

    Returns the bidimensional Espinosa entropy estimate (`Esp2D`) 
    estimated for the data matrix (`Mat`) using the default parameters: 
    time delay = 1, tolerance threshold = 20, percentage similarity = 0.7
    logarithm = natural, matrix template size = [floor(H/10) floor(W/10)], 
    (where H and W represent the height (rows) and width (columns) of 
    the data matrix `Mat`) 
    ** The minimum number of rows and columns of `Mat` must be > 10.
      

        Esp2D = EspEn2D(Mat::AbstractArray{T,2} where T<:Real; m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10),             
                                        tau::Int=1, r::Real=20, ps::Float=.7, Logx::Real=exp(1), Lock::Bool=true)

    Returns the bidimensional Espinosa entropy (`Esp2D`) estimates for the data
    matrix (`Mat`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Template submatrix dimensions, an integer scalar (i.e. the same 
              height and width) or a two-element vector of integers [height, width] with a minimum value > 1.
              (default: [floor(H/10) floor(W/10)]) \n
    `tau`   - Time Delay, a positive integer       (default: 1) \n
    `r`     - Tolerance threshold, a positive scalar  (default: 20) \n
    `ps`    - Percentage similarity, a value in range [0 1],  (default: 0.7) \n
    `Logx`  - Logarithm base, a positive scalar    (default: natural)  \n
    `Lock`  - By default, EspEn2D only permits matrices with a maximum
              size of 128 x 128 to prevent memory errors when storing data
              on RAM. e.g. For Mat = [200 x 200], m = 3, and tau = 1, 
              EspEn2D creates a vector of 753049836 elements. 
              To enable matrices greater than [128 x 128] elements,
              set `Lock` to false.  (default: true)  \n
              `WARNING: unlocking the permitted matrix size may cause your Julia
              IDE to crash.`

    #  See also `SampEn2D`, `FuzzEn2D`, `DispEn2D`, `DistEn2D`, `PermEn2D` 

    # References:
        [1] Ricardo Espinosa, et al.,
            "Two-dimensional EspEn: A New Approach to Analyze Image Texture 
            by Irregularity." 
            Entropy,
            23:1261 (2021)

    """
    function EspEn2D(Mat::AbstractArray{T,2} where T<:Real; 
        m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10), 
        tau::Int=1, r::Real=20, ps::Real=0.7, Logx::Real=exp(1), Lock::Bool=true)
        
    NL, NW = size(Mat)
    ((NL > 128 || NW > 128) && Lock) ? 
        error("To prevent memory errors, matrix width & length must have <= 128 elements.
            To estimate EspEn2D  for the current matrix ($NL,$NW) change Lock to 'false'.
            Caution: unlocking the safe matrix size may cause the Julia IDE to crash.") :
            nothing

    length(m)==1 ?  (mL = m; mW = m) : (mL = m[1];  mW = m[2])
    (NL > 10 && NW > 10) ? nothing : 
        error("Number of rows and columns in Mat must be > 10")
    (minimum(m)>1) ? nothing :
        error("m:     must be an integer > 1, or a 2 element tuple of integer values > 1")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (r >= 0) ?  nothing : error("r:     must be a positive value")
    ((ps >= 0) && (ps <= 1)) ?  nothing : error("ps:   must be a value in range [0 1]")
    (Logx>0) ? nothing :  error("Logx:     must be a positive number > 0")
    
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
    Ny = p*(p-1)/2
    Ny > 300000000 ? @warn("Number of pairwise distance calculations is $Ny") : nothing
       
    Cij = -ones(p-1,p-1)
    for k = 1:p-1
        Temp = abs.(X[:,:,k+1:p] .- X[:,:,k]) .<= r
        Cij[1:(end-k+1),k] = sum(Temp,dims=(1,2))    
    end

    Dm = sum((Cij[:]/(mL*mW)).>=ps)/(p*(p-1)/2)
    Esp2D = -log(Logx, Dm)

    return Esp2D
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