module _MvCoSiEn
export MvCoSiEn
using Statistics: std, mean, median
using LinearAlgebra: Diagonal, UpperTriangular
    """
        MCoSi, Bm = MvCoSiEn(Data) 

     Returns the multivariate cosine similarity entropy estimate (`MCoSi`)
     and the corresponding global probabilities (`Bm`) estimated for the 
     M multivariate sequences in `Data` using the default parameters: 
     embedding dimension = 2*ones(M), time delay = ones(M), 
     angular threshold = 0.1, logarithm = 2, data normalization = none, 
    

    !!! note
    
        To maximize the number of points in the embedding process, this algorithm 
        uses N-max(m * tau) delay vectors and _*not*_ N-max(m) * max(tau) as employed 
        in [1][2].
        
    -------------------------------------------------------------

        MCoSi, Bm = MvCoSiEn(Data::AbstractArray{T,2} where T<:Real; m::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, tau::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, r::Real=.1, Logx::Real=2, Norm::Int=0)

     Returns the multivariate cosine similarity entropy estimates (`MSamp`) 
     estimated from the M multivariate data sequences in `Data` using the 
     specified keyword arguments:

     # Arguments:
     `Data`    - Multivariate dataset, NxM matrix of N (>10) observations (rows) and M (cols) univariate data sequences\n
     `m`       - Embedding Dimension, a vector of M positive integers\n
     `tau`     - Time Delay, a vector of M positive integers\n
     `r`     - Angular threshold, a value in range [0 < r < 1]   \n
     `Logx`    - Logarithm base, a positive scalar (enter 0 for natural log) \n
     `Norm`    - Normalisation of `Data`, one of the following integers:\n
                *  [0]  no normalisation - default
                *  [1]  remove median(`Data`) to get zero-median series
                *  [2]  remove mean(`Data`) to get zero-mean series
                *  [3]  normalises each sequence in `Data` to unit variance and zero mean
                *  [4]  normalises each sequence in `Data` values to range [-1 1]
                    
     # See also  `CoSiEn`, `MvDispEn`, `MvSampEn`, `MvFuzzEn`, `MvPermEn`, `MSEn`
            
     # References:
        [1] H. Xiao, T. Chanwimalueang and D. P. Mandic, 
            "Multivariate Multiscale Cosine Similarity Entropy" 
            IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP),
            pp. 5997-6001, doi: 10.1109/ICASSP43922.2022.9747282.

        [2] Xiao, H.; Chanwimalueang, T.; Mandic, D.P., 
            "Multivariate Multiscale Cosine Similarity Entropy and Its 
            Application to Examine Circularity Properties in Division Algebras."
            Entropy 2022, 24, 1287. 

        [3] Ahmed Mosabber Uddin, Danilo P. Mandic
            "Multivariate multiscale entropy: A tool for complexity
            analysis of multichannel data."
            Physical Review E 84.6 (2011): 061918.
    
        [4] Theerasak Chanwimalueang and Danilo Mandic,
            "Cosine similarity entropy: Self-correlation-based complexity
            analysis of dynamical systems."
            Entropy 
            19.12 (2017): 652.
    """
    function MvCoSiEn(Data::AbstractArray{T,2} where T<:Real; m::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, 
        tau::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, r::Real=.1, Logx::Real=2, Norm::Int=0)

    N, Dn = size(Data)
    isnothing(m) ? m = 2*ones(Int, Dn) : nothing
    isnothing(tau) ? tau = ones(Int, Dn) : nothing
    Logx==0 ? Logx = exp(1) : nothing

    (N>10) && (Dn>1) ? nothing : error("Data:   must be an NxM matrix where N>10 and M>1")
    ndims(m)==1  && length(m)==Dn  && all(m.>0) && eltype(m)<:Int ? nothing : error("m:  vector of M positive integers")
    ndims(tau)==1  && length(tau)==Dn  && all(tau.>0) && eltype(tau)<:Int ? nothing : error("tau:  vector of M positive integers")
    (0<r<1) ? nothing :  error("r:     must be a scalar in range 0 < r < 1")
    (Logx>0) ? nothing : error("Logx:  must be a positive number > 0")
    (Norm in collect(0:4)) ? nothing : error("Norm:   must be an integer in range [0 4]")

    if Norm == 1
        Xi = Data .- median(Data,dims=1)
    elseif Norm == 2
        Xi = Data .- mean(Data,dims=1)
    elseif Norm == 3
        Xi = (Data .- mean(Data, dims=1))./std(Data,corrected=false, dims=1)
    elseif Norm == 4
        Xi = (2*(Data .- minimum(Data, dims=1))./(maximum(Data,dims=1).-minimum(Data, dims=1))) .- 1;
    else
        Xi = Data;
    end

    Nx = N-maximum((m.-1).*tau)
    Zm = zeros(Nx,sum(m));
    q=1
    for k = 1:Dn
        for p = 1:m[k]
            Zm[:,q] = Xi[1+(p-1)*tau[k]:Nx+(p-1)*tau[k],  k]
            q += 1
        end    
    end


    Num = Zm*transpose(Zm); 
    Mag = sqrt.(sum(Diagonal(Num),dims=1))[:]
    Den = Mag*transpose(Mag)
    AngDis = round.(acos.(round.(Num./Den,digits=8))/pi,digits=6)
    if maximum(imag.(AngDis)) < (10^-5)
        Bm = (sum(UpperTriangular(AngDis .< r))-Nx)/(Nx*(Nx-1)/2)
    else
        Bm = (sum(UpperTriangular(real.(AngDis) .< r))-Nx)/(Nx*(Nx-1)/2)
        @warn("Complex values ignored.")
    end
    Bm == 1 || Bm == 0  ? MCoSi = NaN :  MCoSi = -(Bm*log(Logx, Bm)) - ((1-Bm)*log(Logx, 1-Bm))
    
    return MCoSi, Bm
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