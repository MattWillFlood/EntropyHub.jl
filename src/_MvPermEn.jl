module _MvPermEn
export MvPermEn
using Combinatorics: permutations
using Statistics: var, mean
using DSP: hilbert, unwrap, angle
"""
    MPerm, MPnorm = MvPermEn(Data) 
    
 Returns the multivariate permutation entropy estimate (`MPerm`) and
 the normalized permutation entropy for the M multivariate sequences in
 `Data` using the default parameters:
 embedding dimension = 2*ones(M,1), time delay = ones(M,1), 
 logarithm = 2, normalisation = w.r.t #symbols (sum(`m-1`))

 ________________________________________________________________________ 
 
!!! note 

    The multivariate permutation entropy algorithm implemented here uses
    multivariate embedding based on Takens' embedding theorem, and follows
    the methods for multivariate entropy estimation through shared spatial 
    reconstruction as originally presented by Ahmed & Mandic [1]. 

    This function does _*NOT*_ use the multivariate permutation entropy 
    algorithm of Morabito et al. (Entropy, 2012) where the entropy values of 
    individual univariate sequences are averaged because such methods do not
    follow the definition of multivariate embedding and therefore do not
    consider cross-channel statistical complexity.

    To maximize the number of points in the embedding process, this
    algorithm uses N- max(tau * m) delay vectors and _*not*_ N-max(m) * max(tau)
    as employed in [1].

 ________________________________________________________________________ 
      
    MPerm, MPnorm = MvPermEn(Data::AbstractArray{T} where T<:Real; m::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, tau::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, Typex::String="none", tpx::Union{Int,Nothing}=nothing, Norm::Bool=false, Logx::Real=2)
                
 Returns the multivariate permutation entropy estimate (`MPerm`) for
 the M multivariate data sequences in `Data` using the specified keyword arguments:

 # Arguments:
 `Data`  - Multivariate dataset, NxM matrix of N (>10) observations (rows) and M (cols) univariate data sequences\n
 `m`     - Embedding Dimension, a vector of M positive integers\n
 `tau`   - Time Delay, a vector of M positive integers\n
 `Typex` - Permutation entropy variation, can be one of the following strings:\n
            {`'modified'`, `'ampaware'`, `'weighted'`, `'edge'`, `'phase'`}
            See the `EntropyHub guide <https://github.com/MattWillFlood/EntropyHub/blob/main/EntropyHub%20Guide.pdf>`_ for more info on MvPermEn variants.    
 `tpx`   - Tuning parameter for associated permutation entropy variation.
            *   [ampaware]  `tpx` is the A parameter, a value in range [0 1] (default = 0.5)
            *   [edge]      `tpx` is the r sensitivity parameter, a scalar > 0 (default = 1)
            *   [phase]     `tpx` is the option to unwrap the phase angle of Hilbert-transformed signal, either [] or 1 (default = 0)\n
 `Norm`  - Normalisation of MPnorm value, a boolean operator:\n
            * false -  normalises w.r.t log(# of permutation symbols [sum(m)-1]) - default
            * true  -  normalises w.r.t log(# of all possible permutations [sum(m)!])
 `Logx`   - Logarithm base, a positive scalar \n 
        
 # See also    `PermEn`, `PermEn2D`, `XPermEn`, `MSEn`, `MvFuzzEn`, `MvSampEn`

 # References:
    [1] Ahmed Mosabber Uddin, Danilo P. Mandic
        "Multivariate multiscale entropy: A tool for complexity
        analysis of multichannel data."
        Physical Review E 84.6 (2011): 061918.

    [2] Christoph Bandt and Bernd Pompe, 
        "Permutation entropy: A natural complexity measure for time series." 
        Physical Review Letters,
        88.17 (2002): 174102.

    [3] Chunhua Bian, et al.,
        "Modified permutation-entropy analysis of heartbeat dynamics."
        Physical Review E
        85.2 (2012) : 021906

    [4] Bilal Fadlallah, et al.,
        "Weighted-permutation entropy: A complexity measure for time 
        series incorporating amplitude information." 
        Physical Review E 
        87.2 (2013): 022911.

    [5] Hamed Azami and Javier Escudero,
        "Amplitude-aware permutation entropy: Illustration in spike 
        detection and signal segmentation." 
        Computer methods and programs in biomedicine,
        128 (2016): 40-51.

    [6] Zhiqiang Huo, et al.,
        "Edge Permutation Entropy: An Improved Entropy Measure for 
        Time-Series Analysis," 
        45th Annual Conference of the IEEE Industrial Electronics Soc,
        (2019), 5998-6003

    [7] Maik Riedl, Andreas Müller, and Niels Wessel,
        "Practical considerations of permutation entropy." 
        The European Physical Journal Special Topics 
        222.2 (2013): 249-262.

    [8] Kang Huan, Xiaofeng Zhang, and Guangbin Zhang,
        "Phase permutation entropy: A complexity measure for nonlinear time
        series incorporating phase information."
        Physica A: Statistical Mechanics and its Applications
        568 (2021): 125686.

"""
    function MvPermEn(Data::AbstractArray{T,2} where T<:Real; m::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, tau::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, 
                Typex::String="none", tpx::Union{Real,Nothing}=nothing, Norm::Bool=false, Logx::Real=2)

    N, Dn = size(Data)
    isnothing(m) ? m = 2*ones(Int, Dn) : nothing
    isnothing(tau) ? tau = ones(Int, Dn) : nothing
    Logx==0 ? Logx = exp(1) : nothing

    (N>10) && (Dn>1) ? nothing : error("Data:   must be an NxM matrix where N>10 and M>1")
    ndims(m)==1  && length(m)==Dn  && all(m.>0) && eltype(m)<:Int ? nothing : error("m:  vector of M positive integers")
    ndims(tau)==1  && length(tau)==Dn  && all(tau.>0) && eltype(tau)<:Int ? nothing : error("tau:  vector of M positive integers")
    (lowercase(Typex) in ["none","modified","ampaware","weighted","edge","phase"]) ?
    nothing : error("Typex:    must be one of the following strings - 'modified','ampaware','weighted','edge','phase'")
    (isnothing(tpx) || tpx>0) ? nothing : error("tpx:   the value of tpx relates to 'Type'.
       See the EntropyHub guide for further info on the 'tpx' value.")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")

    if lowercase(Typex) == "phase"
        Data = angle.(hilbert(Data))
        tpx == 1 ? Data = unwrap(Data) : nothing
    end

    Nx = N-maximum((m.-1).*tau)
    Sx = zeros(Nx,sum(m))
    q=1
    for k = 1:Dn
        for p = 1:m[k]
            Sx[:,q] = Data[1+(p-1)*tau[k]:Nx+(p-1)*tau[k],  k]
            q += 1
        end    
    end

    Temp = sortind(Sx)
    #Px = collect(permutations(collect(1:sum(m))))
    Px = unique(Temp,dims=1)
    #Counter = zeros(Int, size(Px,1))      
    if  lowercase(Typex) == "modified"
        Tx = (diff(sort(Sx,dims=2),dims=2).==0)
        for km = 1:sum(m)-1
            Temp[Tx[:,km],km+1] = Temp[Tx[:,km],km];
        end            
        Px = unique(Temp,dims=1)
        Counter = map(n -> sum(all(Temp .- transpose(Px[n,:]) .==0,dims=2)), 1:size(Px,1))
        Counter = Counter[Counter.!=0]
        Ppi = Counter/sum(Counter)
                
    elseif lowercase(Typex) == "weighted"
        Wj = var(Sx,corrected=false,dims=2)
        Counter = map(n -> sum(Wj[all(Temp .- transpose(Px[n,:]) .==0,dims=2)]), 1:size(Px,1))
        Counter = Counter[Counter.!=0]
        Ppi = Counter/sum(Wj)
        #=for n = 1:size(Px,1)
            Counter[n] = sum(Wj[all(Temp .- transpose(Px[n]) .==0,dims=2)])
        end=#

    elseif lowercase(Typex) == "ampaware"
        isnothing(tpx) ?  tpx = 0.5 : nothing
        tpx<0 || tpx>1 ?  error("When Typex = 'ampaware', the A parameter must be in the range [0 1]") : nothing            
        AA = sum(abs.(Sx),dims=2)
        AB = sum(abs.(diff(Sx,dims=2)),dims=2)
        Ax = (tpx*AA/sum(m)) + ((1-tpx)*AB/(sum(m)-1));                
        #= for n = 1:size(Px,1)
            Counter[n] = sum(Ax[all(Temp.-transpose(Px[n]).==0,dims=2)])
        end =#
        Counter = map(n -> sum(Ax[all(Temp .- transpose(Px[n,:]) .==0,dims=2)]), 1:size(Px,1))
        Counter = Counter[Counter.!=0]
        Ppi = Counter/sum(Ax);
                
    elseif lowercase(Typex) == "edge"
        isnothing(tpx) ? tpx = 1 : nothing
        tpx <=0 ? error("When Typex = 'Edge', the r sensitivity parameter (tpx) must be > 0") : nothing
        Counter = zeros(size(Px,1))      
        for n in eachindex(1:size(Px,1))
            Tx = diff(Sx[all(Temp .- transpose(Px[n,:]) .==0,dims=2)[:],:],dims=2)
            Counter[n] = sum(mean(hypot.(Tx,1),dims=2).^tpx)
        end
        Counter = Counter[Counter.!=0]
        Ppi = Counter/sum(Counter)   
    
    else
        Counter = map(n -> sum(all(Temp .- transpose(Px[n,:]) .==0,dims=2)), 1:size(Px,1))
        Counter = Counter[Counter.!=0]
        Ppi = Counter/sum(Counter)
    end
        
    if round(sum(Ppi),digits=5) != 1
        @warn ("Potential error with probability calculation")
    end
            
    MPerm = -sum(Ppi.*(log.(Logx, Ppi)));
    Norm ? Pnorm = MPerm/(log(Logx, factorial(sum(m)))) : Pnorm = MPerm/(sum(m)-1)
    
    return MPerm, Pnorm
    end
    
    function sortind(X)
        Y = zeros(Int, size(X))
        for k = 1:length(X[:,1])
            Y[k,:] = sortperm(X[k,:])
        end
        return Y
    end
end

"""Copyright 2024 Matthew W. Flood, EntropyHub
  
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

For Terms of Use see https://github.com/MattWillFlood/EntropyHub"""