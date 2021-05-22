module _PhasEn
export PhasEn
using Plots
using Random: randperm
    """
        Phas = PhasEn(Sig) 
   
    Returns the phase entropy (`Phas`) estimate of the data sequence (`Sig`)
    using the default parameters: 
    angular partitions = 4, time delay = 1, logarithm = natural,
 
        Phas = PhasEn(Sig::AbstractArray{T,1} where T<:Real; K::Int=4, tau::Int=1, Logx::Real=exp(1), Norm::Bool=true, Plotx::Bool=false)
   
    Returns the phase entropy (`Phas`) estimate of the data sequence (`Sig`)  
    using the specified 'keyword' arguments:
   
    # Arguments:
    `K`     - Angular partitions (coarse graining), an integer > 1  \n
    `tau`   - Time Delay, a positive integer    \n
    `Logx`  - Logarithm base, a positive scalar \n 
    `Norm`  - Normalisation of `Phas` value:    \n
              [false]  no normalisation
              [true]   normalises w.r.t. the number of partitions Log(`K`)
    `Plotx` - When `Plotx` == true, returns Poicaré plot (default: false)  \n
 
    # See also `SampEn`, `ApEn`, `GridEn`, `MSEn`, `SlopEn`, `CoSiEn`, `BubbEn`
  
    # References:
        [1] Ashish Rohila and Ambalika Sharma,
            "Phase entropy: a new complexity measure for heart rate
            variability." 
            Physiological measurement
            40.10 (2019): 105006.
  
  
    """
    function PhasEn(Sig::AbstractArray{T,1} where T<:Real; K::Int=4, tau::Int=1, 
        Logx::Real=exp(1), Norm::Bool=true, Plotx::Bool=false)

    N = size(Sig,1)
    (N > 10) ? nothing : error("Sig:   must be a numeric vector")
    (K > 1) ? nothing :  error("K:     must be an integer > 1")
    (tau >0) ? nothing : error("tau:   must be an integer > 0")
    (Logx>0) ? nothing : error("Logx:  must be a positive number > 0")

    Yn = Sig[1+2*tau:end] .- Sig[tau+1:end-tau]
    Xn = Sig[tau+1:end-tau] .- Sig[1:end-2*tau]
    Theta_r = atan.(Yn./Xn)
    Theta_r[(Yn.<0) .& (Xn.<0)] .+= pi
    Theta_r[(Yn.<0) .& (Xn.>0)] .+= 2*pi
    Theta_r[(Yn.>0) .& (Xn.<0)] .+= pi

    Limx = ceil.(maximum(abs.(vcat(Yn, Xn))))
    Angs = range(0,2*pi,length=K+1)
    Tx = zeros(Int, K, length(Theta_r))
    Si = zeros(K)
    for n = 1:K
        Temp = (Theta_r .> Angs[n]) .& (Theta_r .< Angs[n+1]);
        Tx[n,Temp] .= 1
        Si[n] = sum(Theta_r[Temp])
    end
    Si = Si[Si.!=0]
    Phas = -sum((Si./sum(Si)).*(log.(Logx, Si./sum(Si))))
    if Norm
        Phas = Phas/(log(Logx, K))
    end

    if Plotx
        Ys = sin.(Angs)*Limx*sqrt(2);
        Xs = cos.(Angs)*Limx*sqrt(2);
        Cols = hcat(zeros(K), repeat(randperm(K)/K,outer=(1,2)))
        Tx = Bool.(Tx)
        xx = plot()
        for n = 1:K
            plot!(Xn[Tx[n,:]], Yn[Tx[n,:]], seriestype =:scatter, 
                markersize = 2,
                markercolor=RGB(Cols[n,:]...,),
                markerstrokecolor=RGB(Cols[n,:]...,))
        end
        #plot!(hcat(zeros(K+1), Xs), hcat(Ys,zeros(K+1)),c=:magenta,
        plot!(hcat(Xs,zeros(K+1))', hcat(Ys,zeros(K+1))',c=:magenta,
            xlim = (-Limx, Limx), ylim = (-Limx, Limx),
            size = (400, 400), legend=false,
            xticks = [-Limx, 0, Limx],  yticks = [-Limx, 0 ,Limx],
            grid = false)
        xlabel!("X ₙ₊ₜ - X ₙ")
        ylabel!("X ₙ₊₂ₜ - X ₙ₊ₜ")
        display(xx)
    end

    return Phas
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