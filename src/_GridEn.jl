module _GridEn 
export GridEn
using Plots
using Plots.PlotMeasures
using StatsBase: Histogram, fit 
    """
        GDE, GDR, _ = GridEn(Sig) 

    Returns the gridded distribution entropy (`GDE`) and the gridded 
    distribution rate (`GDR`) estimated from the data sequence (`Sig`) using 
    the default  parameters:
    grid coarse-grain = 3, time delay = 1, logarithm = base 2
   
        GDE, GDR, PIx, GIx, SIx, AIx = GridEn(Sig)

    In addition to GDE and GDR, GridEn returns the following indices 
    estimated for the data sequence (`Sig`) using the default  parameters:
    [`PIx`]   - Percentage of points below the line of identity (LI)
    [`GIx`]   - Proportion of point distances above the LI
    [`SIx`]   - Ratio of phase angles (w.r.t. LI) of the points above the LI
    [`AIx`]   - Ratio of the cumulative area of sectors of points above the LI

        GDE, GDR, ..., = GridEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=3, tau::Int=1, Logx::Real=exp(1), Plotx::Bool=false)
    
    Returns the gridded distribution entropy (`GDE`) estimated from the data 
    sequence (`Sig`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Grid coarse-grain (m x m sectors), an integer > 1 \n
    `tau`   - Time Delay, a positive integer    \n
    `Logx`  - Logarithm base, a positive scalar \n
    `Plotx` - When Plotx == true, returns gridded Poicaré plot and a bivariate 
              histogram of the grid point distribution (default: false)  \n
 
    # See also `PhasEn`, `CoSiEn`, `SlopEn`, `BubbEn`, `MSEn`

    # References:
        [1] Chang Yan, et al.,
                "Novel gridded descriptors of poincaré plot for analyzing 
                heartbeat interval time-series." 
                Computers in biology and medicine 
                109 (2019): 280-289.

        [2] Chang Yan, et al. 
                "Area asymmetry of heart rate variability signal." 
                Biomedical engineering online 
                16.1 (2017): 1-14.

        [3] Alberto Porta, et al.,
                "Temporal asymmetries of short-term heart period variability 
                are linked to autonomic regulation." 
                American Journal of Physiology-Regulatory, Integrative and 
                Comparative Physiology 
                295.2 (2008): R550-R557.

        [4] C.K. Karmakar, A.H. Khandoker and M. Palaniswami,
                "Phase asymmetry of heart rate variability signal." 
                Physiological measurement 
                36.2 (2015): 303.


    """
    function GridEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=3, tau::Int=1, 
        Logx::Real=exp(1), Plotx::Bool=false)
        
    N = size(Sig,1)
    (N > 10) ? nothing : error("Sig:   must be a numeric vector")
    (m > 1) ? nothing :  error("m:     must be an integer > 1")
    (tau >0) ? nothing : error("tau:   must be an integer > 0")
    (Logx>0) ? nothing : error("Logx:  must be a positive number > 0")

    Sig_n = (Sig.-minimum(Sig))/(maximum(Sig)-minimum(Sig))
    Temp = hcat(Sig_n[1:end-tau], Sig_n[1+tau:end])
    Edges = collect(0:1/m:1); Edges[end] = 1.1
    N = fit(Histogram,(Temp[:,1],Temp[:,2]),(Edges,Edges)).weights
    Pj = reverse(transpose(N),dims=1)/size(Temp,1); Ppi = Pj[Pj.>0];
    if round(sum(Ppi),digits=4) != 1
        @warn("Potential error of estimated probabilities: P = $(sum(Ppi))")
    end
    GDE = -sum(Ppi.*(log.(Ppi)/log(Logx)));
    GDR = sum(N.!=0)/(m*m);

    T2   = atand.(Temp[:,2]./Temp[:,1])
    Dup  = sum(abs.(diff(Temp[T2.>45,:],dims=2)))
    Dtot = sum(abs.(diff(Temp[T2.!=45,:],dims=2)))
    Sup  = sum((T2[T2.>45].-45))
    Stot = sum(abs.(T2[T2.!=45].-45))
    Aup  = sum(abs.((T2[T2.>45].-45).*sqrt.(sum(Temp[T2.>45,:].^2,dims=2))))
    Atot = sum(abs.((T2[T2.!=45].-45).*sqrt.(sum(Temp[T2.!=45,:].^2,dims=2))))

    PIx = 100*sum(T2.<45)/sum(T2.!=45)
    GIx = 100*Dup/Dtot
    SIx = 100*Sup/Stot
    AIx = 100*Aup/Atot

    if Plotx
        ntrvl = range(0,1,length=m+1)
        p1 = plot(Sig_n[1:end-tau],Sig_n[tau+1:end],seriestype =:scatter,
            c=:magenta, markerstrokecolor=:magenta, markersize = 2)
        plot!(hcat(ntrvl,ntrvl)',hcat(zeros(m+1),ones(m+1))',c=:cyan,lw=3)
        plot!(hcat(zeros(m+1),ones(m+1))',hcat(ntrvl,ntrvl)',c=:cyan,lw=3)
        plot!([0, 1],[0, 1],c=:black, legend=false, aspect_ratio=:equal,
            xlim = (0, 1), ylim = (0, 1),
            xticks=[0, 1], yticks=[0, 1],  xlabel= "X ₙ" , ylabel="X ₙ₊ₜ")
    
        p2 = histogram2d(Temp[:,1],Temp[:,2], nbins=ntrvl,
            xticks=[0, 1], yticks=[0, 1], xlabel= "X ₙ" , ylabel="X ₙ₊ₜ",
            xlim = (0, 1), ylim = (0, 1), aspect_ratio=:equal, 
            colorbar_entry=false, seriescolor=:cool,show_empty_bins=true)

        xx = plot(p1,p2,layout=(1,2),legend=false, framestyle=:grid, margin=7mm)
        display(xx)
    end

    return GDE, GDR, PIx, GIx, SIx, AIx
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