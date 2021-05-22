module _PermEn
export PermEn
using Combinatorics: permutations
using Statistics: std, var, mean
    """
        Perm, Pnorm, cPE = PermEn(Sig) 

    Returns the permuation entropy estimates `Perm`, the normalised
    permutation entropy `Pnorm` and the conditional permutation entropy `cPE`
    for `m` = [1,2] estimated from the data sequence `Sig`  using the default 
    parameters: embedding dimension = 2, time delay = 1, logarithm = base 2, 
    normalisation = w.r.t #symbols (`m`-1)
    Note: using the standard PermEn estimation, `Perm` = 0 when `m` = 1.

        Perm, Pnorm, cPE = PermEn(Sig, m)

    Returns the permutation entropy estimates `Perm` estimated from the data
    sequence `Sig` using the specified embedding dimensions = [1,...,`m`] 
    with other default parameters as listed above.

        Perm, Pnorm, cPE = PermEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, Typex::String="none", tpx::Union{Real,Nothing}=nothing, Logx::Real=2, Norm::Bool=false)

    Returns the permutation entropy estimates `Perm` for dimensions = [1,...,`m`]
    estimated from the data sequence `Sig` using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, an integer > 1 \n
    `tau`   - Time Delay, a positive integer\n
    `Logx`  - Logarithm base, a positive scalar (enter 0 for natural log) \n
    `Norm`  - Normalisation of PermEn value:\n
              false -  normalises w.r.t log(# of permutation symbols [m-1]) - default
              true  -  normalises w.r.t log(# of all possible permutations [m!])
              * Note: Normalised permutation entropy is undefined for m = 1.
              ** Note: When Typex = 'uniquant' and Norm = true, normalisation
              is calculated w.r.t. log(tpx^m)\n
    `Typex`  - Permutation entropy variation, one of the following:
              {'none','uniquant','finegrain','modified','ampaware','weighted','edge'}
              See the EntropyHub guide for more info on PermEn variations.    \n
    `tpx`   - Tuning parameter for associated permutation entropy variation.\n
              [uniquant]  'tpx' is the L parameter, an integer > 1 (default = 4).           
              [finegrain] 'tpx' is the alpha parameter, a positive scalar (default = 1)
              [ampaware]  'tpx' is the A parameter, a value in range [0 1] (default = 0.5)
              [edge]      'tpx' is the r sensitivity parameter, a scalar > 0 (default = 1)
              See the EntropyHub guide for more info on PermEn variations.    \n

    # See also `XPermEn`, `MSEn`, `XMSEn`, `SampEn`, `ApEn`, `CondEn`

    # References:
        [1] Christoph Bandt and Bernd Pompe, 
                "Permutation entropy: A natural complexity measure for time series." 
                Physical Review Letters,
                88.17 (2002): 174102.

        [2] Xiao-Feng Liu, and Wang Yue,
                "Fine-grained permutation entropy as a measure of natural 
                complexity for time series." 
                Chinese Physics B 
                18.7 (2009): 2690.

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

        [8] Zhe Chen, et al. 
                "Improved permutation entropy for measuring complexity of time
                series under noisy condition." 
                Complexity 
                1403829 (2019).

        [9] Maik Riedl, Andreas Müller, and Niels Wessel,
                "Practical considerations of permutation entropy." 
                The European Physical Journal Special Topics 
                222.2 (2013): 249-262.


    """
    function PermEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, 
        Typex::String="none", tpx::Union{Real,Nothing}=nothing, Logx::Real=2, Norm::Bool=false)

    Logx == 0  ? Logx = exp(1) : nothing

    N = size(Sig)[1]
    (N>10) ? nothing :   error("Sig:   must be a numeric vector")
    (m > 1) ? nothing :  error("m:     must be an integer > 1")
    (tau > 0) ? nothing :  error("tau:   must be an integer > 0")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")
    (lowercase(Typex) in ["none","uniquant","finegrain","modified","ampaware","weighted","edge"]) ?
     nothing : error("Typex:    must be one of the following strings - 
        'uniquant','finegrain','modified','ampaware','weighted','edge'")
    (isnothing(tpx) || tpx>0) ? nothing : error("tpx:   the value of tpx relates to 'Type'.
        See the EntropyHub guide for further info on the 'tpx' value.")

    Sx = zeros(N,m)
    Perm = zeros(m)
    Pnorm = zeros(m)
    for k = 1:m
        Nx = N-(k-1)*tau
        Sx[1:Nx,k] = Sig[1+(k-1)*tau:N]
        Temp = sortind(Sx[1:Nx,1:k])
        Px = collect(permutations(collect(1:k)))
        Counter = zeros(length(Px))
            
        if lowercase(Typex) == "uniquant"
            Temp = sort(Sx[1:Nx,1:k], dims=2)
            S = zeros(size(Temp))
            if isnothing(tpx)
                tpx = 4;
            elseif tpx <= 1 || typeof(tpx)!=Int
                error("When Typex = 'Uniquant', L parameter (tpx) must be an integer > 1")
            end
            delta = (maximum(Sig)-minimum(Sig))/tpx
            S[:,1] = map(x -> searchsortedfirst(
                minimum(Sig):delta:maximum(Sig),x), Temp[:,1]) .- 1
            #S[findall(S[:,1].==0),1] .+= 1
            S[S[:,1].==0,1] .+= 1        
            if k > 1
                S[:,2:k] = S[:,1] .+ floor.((Temp[:,2:k] .- Temp[:,1])/delta)
            end
            Px = unique(S,dims=1)
            Counter = zeros(size(Px,1))
            for n = 1:size(Px,1)
                Counter[n] = sum(all(S .- transpose(Px[n,:]).==0,dims=2))
            end
            Counter = Counter[Counter.!=0]
            Ppi = Counter/sum(Counter)
            Norm = true  # Might need to change back to 3         

        elseif lowercase(Typex) == "finegrain"
            if k > 1
                if isnothing(tpx)
                    tpx = 1;
                elseif tpx <= 0
                    error("When Typex = 'finegrain', Alpha parameter (tpx) must be greater than 0")
                end
                q =  floor.(maximum(abs.(diff(Sx[1:Nx,1:k],dims=2)),dims=2)./
                                (tpx*std(abs.(diff(Sig)),corrected=false)))
                Temp = hcat(Temp, q)
                Px = unique(Temp,dims=1)
                Counter = zeros(size(Px,1))
                for n = 1:size(Px,1)
                    Counter[n] = sum(all(Temp .- transpose(Px[n,:]).== 0,dims=2))
                end
                Counter = Counter[Counter.!=0]
                Ppi = Counter./sum(Counter)
                #clear q n qt
            else
                Ppi = 1
            end
                    
        elseif lowercase(Typex) == "modified"
            Tx = (diff(sort(Sx[1:Nx,1:k],dims=2),dims=2).==0)
            for km = 1:k-1
                Temp[Tx[:,km],km+1] = Temp[Tx[:,km],km];
            end            
            Px = unique(Temp,dims=1)
            Counter = zeros(size(Px,1))
            for n = 1:size(Px,1)
                Counter[n] = sum(all(Temp .- transpose(Px[n,:]) .==0,dims=2))
            end
            Counter = Counter[Counter.!=0]
            Ppi = Counter/sum(Counter)
            #clear Tx km
                    
        elseif lowercase(Typex) == "weighted"
            if k > 1
                Wj = var(Sx[1:Nx,1:k],corrected=false,dims=2)
                for n = 1:size(Px,1)
                    Counter[n] = sum(Wj[all(Temp .- transpose(Px[n]) .==0,dims=2)])
                end
                Counter = Counter[Counter.!=0]
                Ppi = Counter/sum(Wj)
                #clear Wj n
            else
                Ppi = 1;
            end
                    
        elseif lowercase(Typex) == "ampaware"
            if k > 1
                if isnothing(tpx)
                    tpx = 0.5;
                elseif tpx<0 || tpx>1
                    error("When Typex = 'ampaware', the A parameter must be in the range [0 1]")
                end
                AA = sum(abs.(Sx[1:Nx,1:k]),dims=2)
                AB = sum(abs.(diff(Sx[1:Nx,1:k],dims=2)),dims=2)
                Ax = (tpx*AA/k) + ((1-tpx)*AB/(k-1));                
                for n = 1:size(Px,1)
                    Counter[n] = sum(Ax[all(Temp.-transpose(Px[n]).==0,dims=2)])
                end
                Counter = Counter[Counter.!=0]
                Ppi = Counter/sum(Ax);
            else
                Ppi = 1;
            end
            #clear AA AB Ax
                    
        elseif lowercase(Typex) == "edge"
            if isnothing(tpx)
                tpx = 1;
            elseif tpx <=0 
                error("When Typex = 'Edge', the r sensitivity parameter (tpx) must be > 0")
            end
            if k > 1
                for n = 1:size(Px,1)
                    Sy = Sx[1:Nx,1:k]
                    Tx = diff(Sy[all(Temp .- transpose(Px[n]) .==0,dims=2)[:],:],dims=2)
                    Counter[n] = sum(mean(hypot.(Tx,1),dims=2).^tpx)
                end
                Counter = Counter[Counter.!=0]
                Ppi = Counter/sum(Counter)
            else Ppi = 1;
            end    
        
        else
            for n = 1:size(Px,1)
                Counter[n] = sum(all(Temp .- transpose(Px[n]).==0,dims=2));
                #sum(all(Temp .- transpose(Px[n]).==0,dims=2));
            end
            Counter = Counter[Counter.!=0]
            Ppi = Counter/sum(Counter)
        end
            
        if round(sum(Ppi),digits=3) != 1
            @warn ("Potential error with probability calculation")
        end
            
        Perm[k] = -sum(Ppi.*(log.(Logx, Ppi)));
        if Norm
            if Norm && lowercase(Typex)=="uniquant"
                Pnorm[k] = Perm[k]/(log(Logx, tpx^k));
            else
                Pnorm[k] = Perm[k]/(log(Logx, factorial(k)));
            end
        else
            Pnorm[k] = Perm[k]/(k-1);
        end     
    end
    cPE = diff(Perm);

    return Perm, Pnorm, cPE
    end
    
    
    function sortind(X)
        Y = zeros(Int, size(X))
        for k = 1:length(X[:,1])
            Y[k,:] = sortperm(X[k,:])
        end
        return Y
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