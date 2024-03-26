module _MSEn
export MSEn, EMD
using Statistics: std, mean, median, var
using DataInterpolations: CubicSpline
#using Dierckx: Spline1D
using Plots
using DSP: conv

#=function __init__()
 @warn("\n\n Methodx option IMF (Intrinisic Mode Function) is not stable.
 Random or highly aperiodic signals may not decompose fully.
 Access to the IMFs decomposed by the empirical mode decomposition (EMD) function
 can be found by calling EntropyHub.EMD(`Sig`,`MaxIMFs`).
 A stable EMD function will be included in future releases.\n\n")
end=#

    """
        MSx, CI = MSEn(Sig, Mobj) 

    Returns a vector of multiscale entropy values `MSx` and the complexity 
    index `CI` of the data sequence `Sig` using the parameters specified 
    by the multiscale object `Mobj` over 3 temporal scales with coarse-
    graining (default). 

        MSx, CI = MSEn(Sig::AbstractArray{T,1} where T<:Real, Mobj::NamedTuple; Scales::Int=3, 
                             Methodx::String="coarse", RadNew::Int=0, Plotx::Bool=false)

    Returns a vector of multiscale entropy values `MSx` and the complexity 
    index `CI` of the data sequence `Sig` using the parameters specified by
    the multiscale object `Mobj` and the following 'keyword' arguments:

    # Arguments:
    `Scales`   - Number of temporal scales, an integer > 1   (default: 3) \n
    `Method`   - Graining method, one of the following: 
                 {`coarse`,`modified`,`imf`,`timeshift`,`generalized`} [default = `coarse`]  
                 For further info on these graining procedures, see the EntropyHub guide.  \n
    `RadNew`   - Radius rescaling method, an integer in the range [1 4].
                 When the entropy specified by `Mobj` is `SampEn` or `ApEn`, 
                 RadNew allows the radius threshold to be updated at each 
                 time scale (Xt). If a radius value is specified by `Mobj` (`r`),
                 this becomes the rescaling coefficient, otherwise it is set
                 to 0.2 (default). The value of RadNew specifies one of the 
                 following methods:\n
                 [1]    Standard Deviation          - r*std(Xt)\n
                 [2]    Variance                    - r*var(Xt) \n
                 [3]    Mean Absolute Deviation     - r*mean_ad(Xt) \n
                 [4]    Median Absolute Deviation   - r*med_ad(Xt)\n
    `Plotx`    - When Plotx == true, returns a plot of the entropy value at each
                 time scale (i.e. the multiscale entropy curve) [default: false]\n


   `For further info on these graining procedures see the EntropyHub guide.`

    # See also `MSobject`, `rMSEn`, `cMSEn`, `hMSEn`, `SampEn`, `ApEn`, `XMSEn`

    # References:
        [1] Madalena Costa, Ary Goldberger, and C-K. Peng,
                "Multiscale entropy analysis of complex physiologic time series."
                Physical review letters
                89.6 (2002): 068102.

        [2] Vadim V. Nikulin, and Tom Brismar,
                "Comment on “Multiscale entropy analysis of complex physiologic
                time series”." 
                Physical review letters 
                92.8 (2004): 089803.

        [3] Madalena Costa, Ary L. Goldberger, and C-K. Peng. 
                "Costa, Goldberger, and Peng reply." 
                Physical Review Letters
                92.8 (2004): 089804.

        [4] Madalena Costa, Ary L. Goldberger and C-K. Peng,
                "Multiscale entropy analysis of biological signals." 
                Physical review E 
                71.2 (2005): 021906.

        [5] Ranjit A. Thuraisingham and Georg A. Gottwald,
                "On multiscale entropy analysis for physiological data."
                Physica A: Statistical Mechanics and its Applications
                366 (2006): 323-332.

        [6] Meng Hu and Hualou Liang,
                "Intrinsic mode entropy based on multivariate empirical mode
                decomposition and its application to neural data analysis." 
                Cognitive neurodynamics
                5.3 (2011): 277-284.

        [7] Anne Humeau-Heurtier 
                "The multiscale entropy algorithm and its variants: A review." 
                Entropy 
                17.5 (2015): 3110-3123.

        [8] Jianbo Gao, et al.,
                "Multiscale entropy analysis of biological signals: a 
                fundamental bi-scaling law." 
                Frontiers in computational neuroscience 
                9 (2015): 64.

        [9]  Paolo Castiglioni, et al.,
                "Multiscale Sample Entropy of cardiovascular signals: Does the
                choice between fixed-or varying-tolerance among scales 
                influence its evaluation and interpretation?." 
                Entropy
                19.11 (2017): 590.

        [10] Tuan D Pham,
                "Time-shift multiscale entropy analysis of physiological signals." 
                Entropy 
                19.6 (2017): 257.

        [11] Hamed Azami and Javier Escudero,
                "Coarse-graining approaches in univariate multiscale sample 
                and dispersion entropy." 
                Entropy 20.2 (2018): 138.

        [12] Madalena Costa and Ary L. Goldberger,
                "Generalized multiscale entropy analysis: Application to quantifying 
                the complex volatility of human heartbeat time series." 
                Entropy 17.3 (2015): 1197-1203.

    """
    function MSEn(Sig::AbstractArray{T,1} where T<:Real, Mobj::NamedTuple; Scales::Int=3, 
        Methodx::String="coarse", RadNew::Int=0, Plotx::Bool=false)

    (size(Sig,1)>10) ? nothing : error("Sig:   must be a numeric vector" )
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
                with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (lowercase(Methodx) in ["coarse","modified","imf","timeshift","generalized"]) ? nothing :
        error("Method:  must be one of the following string names - 'coarse','modified','imf','timeshift','generalized'")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("SampEn","ApEn"))) ? nothing :
        error("RadNew:     must be 0, or an integer in range [1 4] with 
                entropy function 'SampEn' or 'ApEn'")
        
    lowercase(String(Symbol(Mobj.Func))[1]) == 'x' ? error("Base entropy estimator is a cross-entropy method. 
    To perform multiscale CROSS-entropy estimation, use rXMSEn.") : nothing
    
    String(Symbol(Mobj.Func))=="SampEn" ? Mobj = merge(Mobj,(Vcp=false,)) : nothing

    if lowercase(Methodx)=="imf"
        Sig, _  = EMD(Sig,Scales-1)
       # sum(all(Sig.==0,dims=2))==0 ? nothing : Sig = Sig[all(Sig.!=0,dims=2),:]
        #Scales >= size(Sig,1) ? nothing : 
        #@warn("Max number of IMF's decomposed from EMD is less than number of Scales.
        #    MSEn evaluated over $(size(Sig,1)) scales instead of $Scales.")
        #Scales = size(Sig,1)
    
        sum(all(Sig.==0,dims=2))==0 ? Tp1 = ones(Int,size(Sig,1)) : Tp1 = vec(collect(all(Sig .!= 0, dims=2)));
        if !all(Bool.(Tp1)) || (size(Sig,1)<Scales) 
            Sig = Sig[Bool.(Tp1),:]
            @warn("Max number of IMF's decomposed from EMD is less than number of Scales.
                MSEn evaluated over $(size(Sig,1)) scales instead of $Scales.")
            Scales = size(Sig,1);
        end  
    end

    MSx = zeros(Scales)
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
    Func2 = getfield(_MSEn,Symbol(lowercase(Methodx)))
    if RadNew > 0
        if RadNew == 1
            Rnew = x -> std(x, corrected=false)
        elseif RadNew == 2
            Rnew = x -> var(x, corrected=false)
        elseif RadNew == 3
            Rnew = x -> mean(abs.(x .- mean(x)))
        elseif RadNew == 4
            Rnew = x -> median(abs.(x .- median(x)))    
        end

        if haskey(Mobj,:r)
            Cx = Mobj.r
        else
            Cy = ("Standard Deviation","Variance","Mean Abs Deviation",
                    "Median Abs Deviation")
            @warn("No radius value provided in Mobj.
                Default set to 0.2*$(Cy[RadNew]) of each new time-series.")            
            Cx = .2
        end
    end

    for T = 1:Scales
        print(". ")
        Temp = Func2(Sig,T) 

        if lowercase(Methodx) == "timeshift"
            Tempx = zeros(T)
            for k = 1:T
                RadNew > 0 ? Args = (Args..., r=Cx*Rnew(Temp[k,:])) : nothing
                Tempy = Mobj.Func(Temp[k,:]; Args...)
                typeof(Tempy)<:Tuple ? Tempx[k] = Tempy[1][end] : Tempx[k] = Tempy[end]
                # Tempx[k] = Tempy[end]  
            end
            Temp2 = mean(Tempx)
        else
            RadNew > 0 ? Args = (Args..., r=Cx*Rnew(Temp[:])) : nothing   
            Tempx = Mobj.Func(Temp[:]; Args...)
            typeof(Tempx)<:Tuple ? Temp2 = Tempx[1][end] : Temp2 = Tempx[end]
            #Temp2 = Tempx[1][end]  
        end

        MSx[T] = Temp2
    end
    CI = sum(MSx)
    print("\n")
    if any(isnan.(MSx))
        println("Some entropy values may be undefined.")
    end
   

    if Plotx

        p1 = plot(1:Scales, MSx, c=RGB(8/255, 63/255, 77/255), lw=3)
        scatter!(1:Scales, MSx, markersize=6, c=RGB(1, 0, 1),
        xlabel = "Scale Factor", ylabel = "Entropy Value", 
        guidefont = font(12, "arial", RGB(7/255, 54/255, 66/255)),
        tickfontsize = 10, tickfontfamily="arial", legend=false,
        title = "Multiscale $(Mobj.Func) ($(titlecase(Methodx))-graining method)",
        plot_titlefontsize=16, plot_titlefontcolor=RGB(7/255, 54/255, 66/255)) #ylim=(0,maximum(MSx)+.2),
        display(p1)
    end

    return MSx, CI
    end


    function coarse(Z,sx)
        Ns = Int(floor(size(Z,1)/sx))
        Y = mean(reshape(Z[1:sx*Ns],sx,Ns),dims=1)
        return Y
    end
    function modified(Z,sx)
       #= Ns = size(Z,1) - sx + 1
        Y = zeros(Ns)
        for k = 1:Ns
            Y[k] = mean(Z[k:k+sx-1])
        end=#
        Y = (conv(Z,ones(Int,sx))/sx)[sx:end-sx+1]
        return Y 
    end
    function imf(Z,sx)
        Y = sum(Z[1:sx,:],dims=1)
        return Y
    end
    function timeshift(Z,sx)
        Y =  reshape(Z[1:Int(sx*floor(length(Z)/sx))],
                (sx,Int(floor(length(Z)/sx))))
        return Y
    end
    function generalized(Z,sx)
        Ns = floor(Int, length(Z)/sx)
        Y = var(reshape(Z[1:sx*Ns],sx,Ns)',corrected=false,dims=2)
        return Y
    end
    

    function PkFind(X)
        Nx = length(X)
        Indx = zeros(Int,Nx);
        for n = 2:Nx-1
            if X[n-1]< X[n] > X[n+1]
                Indx[n] = n
    
            elseif X[n-1] < X[n] == X[n+1]
                k = 1
                Indx[n] = n
                while (n+k)<Nx && X[n] == X[n+k]
                    Indx[n+k] = n+k
                    k+=1
                end
                n+=k
            end
        end
        Indx = Indx[Indx.!==0]
        return Indx
    end
    
    function EMD(X, Scales::Int)
        Xt = copy(X); N = size(Xt,1); n=1; IMFs = zeros(Scales+1,N)
        MaxER = 20;  MinTN = 2;  #Xt .-= mean(Xt)
        r1 = Xt
    
        while n <= Scales 
            r0 = Xt
            x = 0;            
            Upx = PkFind(r0);  Lwx = PkFind(-r0)  
            UpEnv = CubicSpline(r0[Upx], Upx)  #Spline1D(Upx,r0[Upx],k=3,bc="nearest")
            LwEnv = CubicSpline(r0[Lwx], Lwx)  #Spline1D(Lwx,r0[Lwx],k=3,bc="nearest")
            r1 = r0.- (UpEnv(1:N) .+ LwEnv(1:N))./2 #r0.- (UpEnv.(1:N) .+ LwEnv.(1:N))./2   
            RT = (sum(r0.*r0) - sum(r1.*r1))/sum(r0.*r0)
            length(vcat(Upx,Lwx)) <= MinTN ? (LOG = "Decomposition hit minimal extrema criteria."; break) : nothing
    
            while x < 100 && RT > 0.2
                r0 = 1*r1
                Upx = PkFind(r0);  Lwx = PkFind(-r0)  
                UpEnv = CubicSpline(r0[Upx], Upx)  #Spline1D(Upx,r0[Upx],k=3,bc="nearest")
                LwEnv = CubicSpline(r0[Lwx], Lwx)  #Spline1D(Lwx,r0[Lwx],k=3,bc="nearest")
                r1 = r0.- (UpEnv(1:N) .+ LwEnv(1:N))./2 #r0.- (UpEnv.(1:N) .+ LwEnv.(1:N))./2     
                RT = (sum(r0.*r0) - sum(r1.*r1))/sum(r0.*r0)
                x += 1;
                10*log10(sqrt(sum(r0.*r0))/sqrt(sum(r1.*r1))) > MaxER ? 
                (LOG = "Decomposition hit energy ratio criteria."; break) : nothing
            end
    
            IMFs[n,:] = r1
            Xt .-= r1
            IMFs[Scales+1,:] = r0 .+ mean(X)
            n+=1
        end
        LOG = "All went well :) "
        return IMFs, LOG
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