module _DispEn
export DispEn
using Clustering: kmeans, assignments
using Statistics: std, mean
using StatsFuns: normcdf
    """
        Dispx, RDE = DispEn(Sig)

    Returns the dispersion entropy (`Dispx`) and the reverse dispersion entropy
    (`RDE`) estimated from the data sequence (`Sig`) using the default parameters:
    embedding dimension = 2, time delay = 1, symbols = 3, logarithm = natural,
    data transform = normalised cumulative density function (ncdf)

        Dispx, RDE = DispEn(Sig::AbstractArray{T,1} where T<:Real; c::Int=3, m::Int=2, tau::Int=1, Typex::String="ncdf", Logx::Real=exp(1), Fluct::Bool=false, Norm::Bool=false, rho::Real=1)

    Returns the dispersion entropy (`Dispx`) and the reverse dispersion entropy (`RDE`)
    estimated from the data sequence (`Sig`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer\n
    `tau`   - Time Delay, a positive integer\n
    `c`     - Number of symbols, an integer > 1\n
    `Typex` - Type of data-to-symbolic sequence transform, one of the following:
              {`"linear", "kmeans" ,"ncdf", "finesort", "equal"`}\n
              See the EntropyHub guide for more info on these transforms.\n
    `Logx`  - Logarithm base, a positive scalar\n
    `Fluct` - When Fluct == true, DispEn returns the fluctuation-based
              Dispersion entropy.   [default: false]\n
    `Norm`  - Normalisation of Dispx and RDE value:
              [false]  no normalisation - default
              [true]   normalises w.r.t number of possible dispersion patterns 
                       (c^m  or (2c -1)^m-1 if Fluct == true).\n
    `rho`   - *If Typex == 'finesort', rho is the tuning parameter* (default: 1)\n

    # See also `PermEn`, `SyDyEn`, `MSEn`

    # References:
        [1] Mostafa Rostaghi and Hamed Azami,
            "Dispersion entropy: A measure for time-series analysis." 
            IEEE Signal Processing Letters 
            23.5 (2016): 610-614.

        [2] Hamed Azami and Javier Escudero,
            "Amplitude-and fluctuation-based dispersion entropy." 
            Entropy 
            20.3 (2018): 210.

        [3] Li Yuxing, Xiang Gao and Long Wang,
            "Reverse dispersion entropy: A new complexity measure for 
            sensor signal." 
            Sensors 
            19.23 (2019): 5203.

        [4] Wenlong Fu, et al.,
            "Fault diagnosis for rolling bearings based on fine-sorted 
            dispersion entropy and SVM optimized with mutation SCA-PSO."
            Entropy
            21.4 (2019): 404.

    """
    function DispEn(Sig::AbstractArray{T,1} where T<:Real; c::Int=3, m::Int=2, tau::Int=1, Typex::String="ncdf", 
        Logx::Real=exp(1), Fluct::Bool=false, Norm::Bool=false, rho::Real=1)

    N = size(Sig)[1]
    (N > 10) ? nothing :  error("Sig:   must be a numeric vector")
    (c > 1) ? nothing :   error("c:     must be an integer > 1")
    (m > 0) ? nothing :   error("m:     must be an integer > 0")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (lowercase(Typex) in ["linear", "kmeans", "ncdf", "finesort","equal"]) ? nothing :
        error("Typex:    must be one of the following strings - 'linear','kmeans','ncdf','finesort','equal'")    
    (Logx>0) ? nothing :  error("Logx:  must be a positive number > 0")
    (rho>=0) ? nothing :  error("rho:   must be a positive scalar.")

    if lowercase(Typex) == "linear"    
        Edges = range(minimum(Sig),maximum(Sig),length=c+1)
        Zi = map(x -> sum(Edges[1:c].<=x), Sig)

    elseif lowercase(Typex) == "kmeans"
        Temp = kmeans(transpose(Sig), c; maxiter=200)
        Zx = assignments(Temp)
        Clux = Temp.centers
        xx = sortperm(Clux[:]);   Zi = zeros(N)
        for k = 1:c
            Zi[Zx.==xx[k]] .= k;
        end
            
    elseif lowercase(Typex) == "ncdf"
        Zx = normcdf.(mean(Sig),std(Sig,corrected=false),Sig);
        #= Zi= map(x -> searchsortedfirst(range(0,1,length=c+1),x), Zx) .- 1
        Zi[Zi.==0] .= 1 =#
        Zi = map(x -> sum(range(0,1,length=c+1)[1:c].<=x), Zx)

    elseif lowercase(Typex)  == "finesort"
        Zx = normcdf.(mean(Sig),std(Sig,corrected=false),Sig)
        Zi = map(x -> sum(range(0,1,length=c+1)[1:c].<=x), Zx)
        Ym = zeros(N-(m-1)*tau, m)
        for n = 1:m
            Ym[:,n] = Zx[1+(n-1)*tau:N-((m-n)*tau)]
        end
        Yi = floor.(maximum(abs.(diff(Ym,dims=2)),dims=2)./(rho*std(abs.(diff(Sig)),corrected=false)))

    elseif lowercase(Typex) == "equal"
        ix = sortperm(Sig,alg=MergeSort);
        xx = Int.(round.(range(0,N,length=c+1)))
        Zi = zeros(N)
        for k = 1:c
            Zi[ix[xx[k]+1:xx[k+1]]] .= k
        end
    end

    Zm = zeros(N-(m-1)*tau, m)
    for n = 1:m
        Zm[:,n] = Zi[1+(n-1)*tau:N-((m-n)*tau)]
    end

    (lowercase(Typex) == "finesort") ? Zm = hcat(Zm, Yi) : nothing
    if Fluct
        Zm = diff(Zm,dims=2)
        (m < 2) ? @warn(["Fluctuation-based Dispersion Entropy is undefined for m = 1. "...
                "An embedding dimension (m) > 1 should be used."]) : nothing
    end

    T = unique(Zm,dims=1)
    Nx = size(T)[1]
    Counter = zeros(Nx)
    for n = 1:Nx
        Counter[n] = sum(all(Zm .- transpose(T[n,:]) .==0, dims=2))
    end    
    
    Ppi = Counter[Counter.!= 0]/size(Zm)[1]

    if Fluct
        RDE = sum((Ppi .- (1/((2*c - 1)^(m-1)))).^2)
    else
        RDE = sum((Ppi .- (1/(c^m))).^2)
    end
    #RDE = sum(Ppi.^2) - (1/Nx)
    
    if round(sum(Ppi)) != 1
        @warn("Potential Error calculating probabilities")
    end

    Dispx = -sum(Ppi.*log.(Logx, Ppi))
    if Norm
        #Dispx = Dispx/log(Logx, Nx)
        #RDE = RDE/(1 - (1/Nx))
        if Fluct
            Dispx = Dispx/(log(Logx, (2*c - 1)^(m-1)))
            RDE = RDE/(1 - (1/((2*c - 1)^(m-1))))
        else
            Dispx = Dispx/(log(Logx, c^m))
            RDE = RDE/(1 - (1/(c^m)))
        end
    end

    return Dispx, RDE
    end

end

#= for n = 1:m
    Zm[:,n] = Zi[1+(n-1)*tau:N-((m-n)*tau)]
    T[:,n] = repeat(1:c,inner=(c^(n-1),1),outer=(c^(m-n),1))
end

if lowercase(Typex) == "finesort"
    Zm = hcat(Zm, Yi)
    temp = sort(unique(Yi))
    T = repeat(T,outer=(length(temp),1))
    T = hcat(T,repeat(temp,inner=Nx))
    Counter = repeat(Counter,outer=length(Yi))
    Nx = length(Counter) 
end


if Fluct
    Zm = diff(Zm,dims=2)
    T  = unique(Zm,dims=1)
    Nx = size(T,1)
    Counter = zeros(Nx)
end =#

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