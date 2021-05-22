module _SampEn2D   
export SampEn2D
using Statistics: mean, std
    """
        SE2D, Phi1, Phi2 = SampEn2D(Mat) 

    Returns the bidimensional sample entropy estimate (`SE2D`) and the number
    of matched sub-matricess (m:Phi1, m+1:Phi2) estimated for the data 
    matrix (`Mat`) using the default parameters: time delay = 1,
    radius distance threshold = 0.2*SD(`Mat`), logarithm = natural
    matrix template size = [floor(H/10) floor(W/10)]  (where H and W represent
    the height (rows) and width (columns) of the data matrix `Mat`) 
    ** The minimum dimension size of Mat must be > 10.**

        SE2D, Phi1, Phi2 = SampEn2D(Mat::AbstractArray{T,2} where T<:Real; m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10),             
                                        tau::Int=1, r::Real=0.2*std(Mat,corrected=false), Logx::Real=exp(1), Lock::Bool=true)

    Returns the bidimensional sample entropy (`SE2D`) estimates for the data
    matrix (`Mat`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Template submatrix dimensions, an integer scalar (i.e. the same 
              height and width) or a two-element vector of integers 
              [height, width] with a minimum value > 1.
              (default: [floor(H/10) floor(W/10)]) \n
    `tau`   - Time Delay, a positive integer       (default: 1) \n
    `r`     - Distance Threshold, a positive scalar  (default: 0.2*SD(Mat)) \n
    `Logx`  - Logarithm base, a positive scalar    (default: natural)  \n
    `Lock`  - By default, SampEn2D only permits matrices with a maximum
              size of 128 x 128 to prevent memory errors when storing data
              on RAM. e.g. For Mat = [200 x 200], m = 3, and tau = 1, 
              SampEn2D creates a vector of 753049836 elements. 
              To enable matrices greater than [128 x 128] elements,
              set `Lock` to false.  (default: true)  \n
              `WARNING: unlocking the permitted matrix size may cause your Julia
              IDE to crash.`

    #  See also `SampEn`, `FuzzEn2D`, `XSampEn`, `MSEn` 

    # References:
        [1] Luiz Eduardo Virgili Silva, et al.,
             "Two-dimensional sample entropy: Assessing image texture 
             through irregularity." 
             Biomedical Physics & Engineering Express
             2.4 (2016): 045002.


    """
    function SampEn2D(Mat::AbstractArray{T,2} where T<:Real; 
        m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10), 
        tau::Int=1, r::Real=0.2*std(Mat,corrected=false), Logx::Real=exp(1), Lock::Bool=true)
        
    NL, NW = size(Mat)
    ((NL > 128 || NW > 128) && Lock) ? 
        error("To prevent memory errors, matrix width & length must have <= 128 elements.
            To estimate SampEn2D  for the current matrix ($NL,$NW) change Lock to 'false'.
            Caution: unlocking the safe matrix size may cause the Julia IDE to crash.") :
            nothing

    length(m)==1 ?  (mL = m; mW = m) : (mL = m[1];  mW = m[2])
    (NL > 10 && NW > 10) ? nothing : 
        error("Number of rows and columns in Mat must be > 10")
    (minimum(m)>1) ? nothing :
        error("m:     must be an integer > 1, or a 2 element tuple of integer values > 1")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (r >= 0) ?  nothing : error("r:     must be a positive value")
    (Logx>0) ? nothing :  error("Logx:     must be a positive number > 0")
    
    NL = NL - mL*tau
    NW = NW - mW*tau
    X = zeros(mL+1,mW+1,NL*NW)
    p = 0
    for k = 1:NL
        for n = 1:NW
            p += 1
            X[:,:,p] = Mat[k:tau:(mL)*tau+k,n:tau:(mW)*tau+n]
        end
    end
    p = size(X,3)
    p != NL*NW ? @warn("Potential error with submatrix division.") : nothing
    Ny = p*(p-1)/2
    Ny > 300000000 ? @warn("Number of pairwise distance calculations is $Ny") : nothing

    Y1 = zeros(p-1)
    Y2 = zeros(p-1)
    for k = 1:p-1
        Temp = maximum(abs.(X[1:mL,1:mW,k+1:p] .- X[1:mL,1:mW,k]), dims=(1,2))[:] .< r
        Y1[k] = sum(Temp);    Temp = findall(Temp.>0) .+ k
        Y2[k] = sum(maximum(abs.(X[:,:,Temp] .- X[:,:,k]),dims=(1,2)) .< r)
    end
    Phi1 = sum(Y1)/Ny
    Phi2 = sum(Y2)/Ny
    SE2D = -log(Logx, Phi2/Phi1)

    return SE2D, Phi1, Phi2
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