module _XMSEn
export XMSEn
using Statistics: std, mean, median, var
using Dierckx: Spline1D
using Plots
#=function __init__()
 @warn("\n\n Methodx option IMF (Intrinisic Mode Function) is not stable.
 Random or highly aperiodic signals may not decompose fully.
 Access to the IMFs decomposed by the empirical mode decomposition (EMD) function
 can be found by calling _MSEn.EMD(`Sig`,`MaxIMFs`).
 A stable EMD function will be included in future releases.\n\n")
end=#
    
    """
        MSx, CI = XMSEn(Sig, Mobj) 

    Returns a vector of multiscale cross-entropy values `MSx` and the complexity 
    index `CI` between the data sequences contained in `Sig` using the parameters 
    specified by the multiscale object `Mobj` over 3 temporal scales with coarse-
    graining `default`. 

        MSx,CI = MSEn(Sig::AbstractArray{T,2} where T<:Real, Mobj::NamedTuple; 
                         Scales::Int=3, Methodx::String="coarse", RadNew::Int=0, Plotx::Bool=false)

    Returns a vector of multiscale cross-entropy values `MSx` and the complexity 
    index `CI` of the data sequences contained in `Sig` using the parameters 
    specified by the multiscale object `Mobj` and the following 'keyword' arguments:

    # Arguments:
    `Scales`   - Number of temporal scales, an integer > 1   (default: 3) \n
    `Method`   - Graining method, one of the following:\n
                 {`"coarse", "modified", "imf", "timeshift"`}  [default: 'coarse'] 
                 For further info on graining procedures, see the Entropyhub guide. \n
    `RadNew`   - Radius rescaling method, an integer in the range [1 4].
                 When the entropy specified by `Mobj` is `SampEn` or `ApEn`, 
                 `RadNew` allows the radius threshold to be updated at each 
                 time scale (Xt). If a radius value is specified by `Mobj` (`r`),
                 this becomes the rescaling coefficient, otherwise it is set
                 to 0.2 (default). The value of `RadNew` specifies one of the 
                 following methods: \n
                 [1]    Standard Deviation          - r*std(Xt) \n
                 [2]    Variance                    - r*var(Xt) \n
                 [3]    Mean Absolute Deviation     - r*mean_ad(Xt) \n
                 [4]    Median Absolute Deviation   - r*med_ad(Xt) \n
    `Plotx`    - When `Plotx` == true, returns a plot of the entropy value at each
                 time scale (i.e. the multiscale entropy curve)  [default: false]
    
    `For further info on these graining procedures see the EntropyHub guide.`

    # See also `MSobject`, `MSEn`, `cXMSEn`, `rXMSEn`, `hXMSEn`, `XSampEn`, `XApEn`, `XFuzzEn`

    # References:
        [1] Rui Yan, Zhuo Yang, and Tao Zhang,
            "Multiscale cross entropy: a novel algorithm for analyzing two
            time series." 
            5th International Conference on Natural Computation. 
            Vol. 1, pp: 411-413 IEEE, 2009.

        [2] Madalena Costa, Ary Goldberger, and C-K. Peng,
            "Multiscale entropy analysis of complex physiologic time series."
            Physical review letters
            89.6 (2002): 068102.

        [3] Vadim V. Nikulin, and Tom Brismar,
            "Comment on “Multiscale entropy analysis of complex physiologic
            time series”." 
            Physical review letters 
            92.8 (2004): 089803.

        [4] Madalena Costa, Ary L. Goldberger, and C-K. Peng. 
            "Costa, Goldberger, and Peng reply." 
            Physical Review Letters
            92.8 (2004): 089804.

        [5] Antoine Jamin, et al,
            "A novel multiscale cross-entropy method applied to navigation 
            data acquired with a bike simulator." 
            41st annual international conference of the IEEE EMBC
            IEEE, 2019.

        [6] Antoine Jamin and Anne Humeau-Heurtier. 
            "(Multiscale) Cross-Entropy Methods: A Review." 
            Entropy 
            22.1 (2020): 45.

    """
    function XMSEn(Sig::AbstractArray{T,2} where T<:Real, Mobj::NamedTuple; 
        Scales::Int=3, Methodx::String="coarse", RadNew::Int=0, Plotx::Bool=false)

    size(Sig,1) == 2 ? Sig = Sig' : nothing

    (size(Sig,1)>10) ? nothing : error("Sig:   must be a numeric vector" )
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
                with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (lowercase(Methodx) in ["coarse","modified","imf","timeshift"]) ? nothing :
        error("Method:  must be one of the following string names -
                'coarse','modified','imf','timeshift'")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("XSampEn","XApEn"))) ? nothing :
        error("RadNew:     must be 0, or an integer in range [1 4] with 
                entropy function 'XSampEn' or 'XApEn'")
        
    if lowercase(Methodx)=="imf" 
        Imfx, _ = EMD(Sig[:,1],Scales-1) 
        Imfy, _ = EMD(Sig[:,2],Scales-1) 
        (Scales >= size(Imfx,1) || Scales >= size(Imfy,1)) ? nothing : 
        @warn("Max number of IMF's decomposed from EMD is less than number of Scales.
            MSEn evaluated over $(size(Sig,1)) scales instead of $Scales.")
        #=Sig = zeros(size(Sig,1), 2, Scales)
        Sig[:,1,:] = transpose(Imfx[1:Scales,:])
        Sig[:,2,:] = transpose(Imfy[1:Scales,:])
        sum(all(Sig[:,1,:].==0,dims=1))==0 ? nothing : 
        Sig = Sig[:, :, vec(collect(all(Sig[:,1,:] .!= 0, dims=1)))]
        sum(all(Sig[:,2,:].==0,dims=1))==0 ? nothing :
        Sig = Sig[:, :, vec(collect(all(Sig[:,2,:] .!= 0, dims=1)))] =#

        Sig = zeros(Scales, size(Sig,1), 2)
        Sig[:,:,1] = Imfx #[1:Scales,:]
        Sig[:,:,2] = Imfy # [1:Scales,:]

        sum(all(Sig[:,:,1].==0,dims=2))==0 ? nothing : 
        Sig = Sig[:, :, vec(collect(all(Sig[:,1,:] .!= 0, dims=2)))]
        sum(all(Sig[:,:,2].==0,dims=2))==0 ? nothing :
        Sig = Sig[:, :, vec(collect(all(Sig[:,2,:] .!= 0, dims=2)))]
        Scales = size(Sig,1)
    end

    MSx = zeros(Scales)
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
    Func2 = getfield(_XMSEn,Symbol(lowercase(Methodx)))
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
            Cy = ("Standard Deviation","Variance","Mean Abs Deviation", "Median Abs Deviation")
            @warn("No radius value provided in Mobj.
                Default set to 0.2*$(Cy[RadNew]) of each new time-series.")            
            Cx = .2
        end
    end

    for T = 1:Scales
        print(" .")
        Temp = Func2(Sig,T) 

        if lowercase(Methodx) == "timeshift"
            Tempx = zeros(T)
            for k = 1:T
                RadNew > 0 ? Args = (Args..., r=Cx*Rnew(Temp[k,:,:])) : nothing
                Tempy = Mobj.Func(Temp[k,:,:]; Args...)
                typeof(Tempy)<:Tuple ? Tempx[k] = Tempy[1][end] : Tempx[k] = Tempy[end]
                # Tempx[k] = Tempy[1][end]
            end
            Temp2 = mean(Tempx)
        else
            RadNew > 0 ? Args = (Args..., r=Cx*Rnew(Temp)) : nothing   
            Tempx = Mobj.Func(Temp; Args...)
            typeof(Tempx)<:Tuple ? Temp2 = Tempx[1][end] : Temp2 = Tempx[end]
            # Temp2 = Tempx[1][end] 
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
        plot_titlefontsize=16, plot_titlefontcolor=RGB(7/255, 54/255, 66/255))
        display(p1)
    end

    return MSx, CI
    end


    function coarse(Z,sx)
        Ns = Int(floor(size(Z,1)/sx))
        T1 = mean(reshape(Z[1:sx*Ns,1],sx,Ns),dims=1)[:]
        T2 = mean(reshape(Z[1:sx*Ns,2],sx,Ns),dims=1)[:]
        return hcat(T1,T2)
    end
    function modified(Z,sx)
        Ns = size(Z,1) - sx + 1
        Y = zeros(Ns,2)
        for k = 1:Ns
            Y[k,1] = mean(Z[k:k+sx-1,1])
            Y[k,2] = mean(Z[k:k+sx-1,2])
        end
        return Y 
    end
    function imf(Z,sx)
        Y = sum(Z[1:sx,:,:],dims=1)
        return Y[1,:,:]
    end
    function timeshift(Z,sx)
        Y = zeros(sx,Int(floor(size(Z,1)/sx)),2)
        Y[:,:,1] =  reshape(Z[1:Int(sx*floor(size(Z,1)/sx)),1],
                    (sx,Int(floor(size(Z,1)/sx))))
        Y[:,:,2] =  reshape(Z[1:Int(sx*floor(size(Z,1)/sx)),2],
                    (sx,Int(floor(size(Z,1)/sx))))
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
            UpEnv = Spline1D(Upx,r0[Upx],k=3,bc="nearest")
            LwEnv = Spline1D(Lwx,r0[Lwx],k=3,bc="nearest")
            r1 = r0.- (UpEnv.(1:N) .+ LwEnv.(1:N))./2     
            RT = (sum(r0.*r0) - sum(r1.*r1))/sum(r0.*r0)
            length(vcat(Upx,Lwx)) <= MinTN ? (LOG = "Decomposition hit minimal extrema criteria."; break) : nothing
    
            while x < 100 && RT > 0.2
                r0 = 1*r1
                Upx = PkFind(r0);  Lwx = PkFind(-r0)  
                UpEnv = Spline1D(Upx,r0[Upx],k=3,bc="nearest")
                LwEnv = Spline1D(Lwx,r0[Lwx],k=3,bc="nearest")
                r1 = r0.- (UpEnv.(1:N) .+ LwEnv.(1:N))./2  
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