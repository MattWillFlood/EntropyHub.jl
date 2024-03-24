module _PermEn2D   
export PermEn2D
    """
        Perm2D = PermEn2D(Mat) 

    Returns the bidimensional permutation entropy estimate (`Perm2D`) estimated for 
    the data matrix (`Mat`) using the default parameters: time delay = 1,
    logarithm = natural, template matrix size = [floor(H/10) floor(W/10)]  
    (where H and W represent the height (rows) and width (columns) of the data matrix `Mat`) \n
    ** The minimum dimension size of Mat must be > 10.**

        Perm2D = PermEn2D(Mat::AbstractArray{T,2} where T<:Real; m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10),             
                                        tau::Int=1, Norm::Bool=true, Logx::Real=exp(1), Lock::Bool=true)

    Returns the bidimensional permutation entropy (`Perm2D`) estimates for the data
    matrix (`Mat`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Template submatrix dimensions, an integer scalar (i.e. the same 
              height and width) or a two-element vector of integers 
              [height, width] with a minimum value > 1.
              (default: [floor(H/10) floor(W/10)]) \n
    `tau`   - Time Delay, a positive integer       (default: 1) \n
    `Norm`  - Normalization of permutation entropy estimate, a boolean  (default: true) \n
    `Logx`  - Logarithm base, a positive scalar    (default: natural)  \n
    `Lock`  - By default, PermEn2D only permits matrices with a maximum
              size of 128 x 128 to prevent memory errors when storing data
              on RAM. e.g. For Mat = [200 x 200], m = 3, and tau = 1, 
              SampEn2D creates a vector of 753049836 elements. 
              To enable matrices greater than [128 x 128] elements,
              set `Lock` to false.  (default: true)  \n
              `WARNING: unlocking the permitted matrix size may cause your Julia
              IDE to crash.`\n\n

    **NOTE** - The original bidimensional permutation entropy algorithms 
              [1][2] do not account for equal-valued elements of the embedding
              matrices. 
              To overcome this, `PermEn2D` uses the lowest common rank for
              such instances. For example, given an embedding matrix A where,
              A = [3.4  5.5  7.3]
                  |2.1  6    9.9|
                  [7.3  1.1  2.1]
              would normally be mapped to an ordinal pattern like so,
              [3.4  5.5  7.3  2.1  6  9.9  7.3  1.1  2.1] =>
              [ 8    4    9    1   2   5    3    7    6 ]
              However, indices 4 & 9, and 3 & 7 have the same values, 2.1
              and 7.3 respectively. Instead, PermEn2D uses the ordinal pattern
              [ 8    4    4    1   2   5    3    3    6 ]
              where the lowest rank (4 & 3) are used instead (of 9 & 7). 
              Therefore, the number of possible permutations is no longer 
              (mx*my)!, but (mx*my)^(mx*my). Here, the PermEn2D value is 
              normalized by the maximum Shannon entropy (Smax = log((mx*my)!) 
              ``assuming that no equal values are found in the permutation
              motif matrices``, as presented in [1].
 

    #  See also `SampEn2D`, `FuzzEn2D`, `DispEn2D`, `DistEn2D` 

    # References:
        [1] Haroldo Ribeiro et al.,
                "Complexity-Entropy Causality Plane as a Complexity Measure 
                for Two-Dimensional Patterns"
                PLoS ONE (2012), 7(8):e40689, 

        [2] Luciano Zunino and Haroldo Ribeiro,
                "Discriminating image textures with the multiscale
                two-dimensional complexity-entropy causality plane"
                Chaos, Solitons and Fractals,  91:679-688 (2016)

        [3] Matthew Flood and Bernd Grimm,
                "EntropyHub: An Open-source Toolkit for Entropic Time Series Analysis"
                PLoS ONE (2021) 16(11): e0259448.
    """
    function PermEn2D(Mat::AbstractArray{T,2} where T<:Real; m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10), 
        tau::Int=1, Norm::Bool=true, Logx::Real=exp(1), Lock::Bool=true)
        
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
    (Logx>0) ? nothing :  error("Logx:  must be a positive number > 0")
    
    NL = NL - (mL-1)*tau
    NW = NW - (mW-1)*tau
    
    Temp = Mat[1:tau:(mL-1)*tau+1,1:tau:(mW-1)*tau+1]
    Sord = sort(Temp[:])
    Dix = sortperm(Temp[:])

    if any(diff(Sord).==0)
        for x in findall(diff(Sord).==0).+1
            Dix[x] = Dix[x-1]        
        end
    end
    Counter = [0]

    for k in 1:NL
        for n in 1:NW        
            Temp = Mat[k:tau:(mL-1)*tau+k,n:tau:(mW-1)*tau+n]
            Sord = sort(Temp[:])
            Dx = sortperm(Temp[:])

            if any(diff(Sord).==0)
                for x in findall(diff(Sord).==0).+1
                    Dx[x] = Dx[x-1]        
                end
            end     
                    
            if any(all((Dix .- Dx).==0,dims=1))
                Counter .+= all((Dix .- Dx).==0,dims=1).*1         
            else
                Dix = [Dix Dx]
                Counter = [Counter 1]            
            end                   
        end
    end

    sum(Counter) != NL*NW ? @warn("Potential error with permutation comparisons.") : nothing
   
    Pi = Counter/sum(Counter)
    Perm2D = -sum(Pi.*log.(Logx, Pi))
    Norm ? Perm2D /= log(Logx, factorial(big(mL*mW))) : nothing

    return Perm2D
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