module _XMSEn
export XMSEn, EMD
using Statistics: std, mean, median, var
# using Dierckx: Spline1D
using DataInterpolations: CubicSpline
using Plots
using DSP: conv

#=function __init__()
 @warn("\n\n Methodx option IMF (Intrinisic Mode Function) is not stable.
 Random or highly aperiodic signals may not decompose fully.
 Access to the IMFs decomposed by the empirical mode decomposition (EMD) function
 can be found by calling _MSEn.EMD(`Sig`,`MaxIMFs`).
 A stable EMD function will be included in future releases.\n\n")
end=#
    
    """
        MSx, CI = XMSEn(Sig1, Sig2, Mobj) 

    Returns a vector of multiscale cross-entropy values `MSx` and the complexity 
    index `CI` between the data sequences contained in `Sig1` and `Sig2` using the parameters 
    specified by the multiscale object `Mobj` over 3 temporal scales with coarse-
    graining `default`. 

        MSx,CI = MSEn(Sig1::AbstractVector{T} where T<:Real, Sig2::AbstractVector{T} where T<:Real, Mobj::NamedTuple; 
                         Scales::Int=3, Methodx::String="coarse", RadNew::Int=0, Plotx::Bool=false)

    Returns a vector of multiscale cross-entropy values `MSx` and the complexity 
    index `CI` of the data sequences contained in `Sig1` and `Sig2` using the parameters 
    specified by the multiscale object `Mobj` and the following 'keyword' arguments:

    # Arguments:
    `Scales`   - Number of temporal scales, an integer > 1   (default: 3) \n
    `Method`   - Graining method, one of the following:\n
                 {`"coarse", "modified", "imf", "timeshift","generalized"`}  [default: 'coarse'] 
                 For further info on graining procedures, see the Entropyhub guide. \n
    `RadNew`   - Radius rescaling method, an integer in the range [1 4].
                 When the entropy specified by `Mobj` is `XSampEn` or `XApEn`, 
                 `RadNew` allows the radius threshold to be updated at each 
                 time scale (Xt). If a radius value is specified by `Mobj` (`r`),
                 this becomes the rescaling coefficient, otherwise it is set
                 to 0.2 (default). The value of `RadNew` specifies one of the 
                 following methods: \n
                 [1]    Pooled Standard Deviation          - r*std(Xt) \n
                 [2]    Pooled Variance                    - r*var(Xt) \n
                 [3]    Total Mean Absolute Deviation      - r*mean_ad(Xt) \n
                 [4]    Total Median Absolute Deviation    - r*med_ad(Xt) \n
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
    function XMSEn(Sig1::AbstractVector{T} where T<:Real, Sig2::AbstractVector{T} where T<:Real, Mobj::NamedTuple; 
                     Scales::Int=3, Methodx::String="coarse", RadNew::Int=0, Plotx::Bool=false)

    (size(Sig1,1)>=10)  && (size(Sig2,1)>=10) ? nothing : error("Sig1/Sig2:   sequences must have >= 10 values")
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
                with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (lowercase(Methodx) in ["coarse","modified","imf","timeshift","generalized"]) ? nothing :
        error("Method:  must be one of the following string names -
                'coarse','modified','imf','timeshift','generalized'")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("XSampEn","XApEn"))) ? nothing :
        error("RadNew:     must be 0, or an integer in range [1 4] with 
                entropy function 'XSampEn' or 'XApEn'")
        
    String(Symbol(Mobj.Func))=="XSampEn" ? Mobj = merge(Mobj,(Vcp=false,)) : nothing

    if lowercase(Methodx)=="imf"
        Imfx, _ = EMD(Sig1,Scales-1) 
        Imfy, _ = EMD(Sig2,Scales-1) 
        Sig1 = Imfx # [1:Scales,:]'
        Sig2 = Imfy # [1:Scales,:]'
        # If any of the IMFs are just zeros, then take only the IMFs that aren't
        sum(all(Sig1.==0,dims=2))==0 ? Tp1 = ones(Int,size(Imfx,1)) : Tp1 = vec(collect(all(Sig1 .!= 0, dims=2)));
        sum(all(Sig2.==0,dims=2))==0 ? Tp2 = ones(Int,size(Imfy,1)) : Tp2 = vec(collect(all(Sig2 .!= 0, dims=2)));
        if !all(Bool.(Tp1.+Tp2.-1)) || (size(Imfx,1)<Scales || size(Imfy,1)<Scales) 
            Sig1 = Sig1[Bool.(Tp1.+Tp2.-1),:]
            Sig2 = Sig2[Bool.(Tp1.+Tp2.-1),:]
            @warn("Max number of IMF's decomposed from EMD is less than number of Scales.
                MSEn evaluated over $(size(Sig,1)) scales instead of $Scales.")
            Scales = size(Sig1,1);
        end
    end

    MSx = zeros(Scales)
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
    Func2 = getfield(_XMSEn,Symbol(lowercase(Methodx)))
    if RadNew > 0
        if RadNew == 1
            Rnew = (x,y) -> sqrt((var(x)*(size(x,1)-1) + var(y)*(size(y,1)-1))/(size(x,1)+size(y,1)-1))
        elseif RadNew == 2
            Rnew = (x,y) -> ((var(x)*(size(x,1)-1) + var(y)*(size(y,1)-1))/(size(x,1)+size(y,1)-1))
        elseif RadNew == 3
            Rnew = (x,y) -> mean(abs.(vcat(x,y) .- mean(vcat(x,y))))
        elseif RadNew == 4
            Rnew = (x,y) -> median(abs.(vcat(x,y) .- median(vcat(x,y))))    
        end

        if haskey(Mobj,:r)
            Cx = Mobj.r
        else
            Cy = ("Pooled Standard Deviation","Pooled Variance","Total Mean Abs Deviation", "Total Median Abs Deviation")
            @warn("No radius value provided in Mobj.
                Default set to 0.2*$(Cy[RadNew]) of each new time-series.")            
            Cx = .2
        end
    end

    for T = 1:Scales
        print(" .")
        TempA, TempB = Func2(Sig1, Sig2, T) 

        if lowercase(Methodx) == "timeshift"
            Tempx = zeros(T)
            for k = 1:T
                RadNew > 0 ? Args = (Args..., r=Cx*Rnew(TempA[k,:], TempB[k,:])) : nothing
                Tempy = Mobj.Func(TempA[k,:], TempB[k,:]; Args...)
                typeof(Tempy)<:Tuple ? Tempx[k] = Tempy[1][end] : Tempx[k] = Tempy[end]
                # Tempx[k] = Tempy[1][end]
            end
            Temp2 = mean(Tempx)
        else
            RadNew > 0 ? Args = (Args..., r=Cx*Rnew(TempA,TempB)) : nothing   
            Tempx = Mobj.Func(TempA, TempB; Args...)
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


    function coarse(Za, Zb, sx)
        Na = Int(floor(size(Za,1)/sx))
        Nb = Int(floor(size(Zb,1)/sx))        
        Y1 = mean(reshape(Za[1:sx*Na],sx,Na),dims=1)[:]
        Y2 = mean(reshape(Zb[1:sx*Nb],sx,Nb),dims=1)[:]
        return Y1, Y2
    end
    function modified(Za, Zb, sx)
      #=  Ns = size(Z,1) - sx + 1
        Y = zeros(Ns,2)
        for k = 1:Ns
            Y[k,1] = mean(Z[k:k+sx-1,1])
            Y[k,2] = mean(Z[k:k+sx-1,2])
        end=#

        Y1 = (conv(Za,ones(Int, sx))/sx)[sx:end-sx+1][:]
        Y2 = (conv(Zb,ones(Int, sx))/sx)[sx:end-sx+1][:]
        return Y1, Y2
    end
    function imf(Za, Zb, sx)
        Y1 = sum(Za[1:sx,:],dims=1)[:]
        Y2 = sum(Zb[1:sx,:],dims=1)[:]
        return Y1, Y2
    end
    function timeshift(Za, Zb, sx)
       # Y1 = zeros(sx,Int(floor(size(Za,1)/sx)))
       # Y2 = zeros(sx,Int(floor(size(Zb,1)/sx)))
        Y1 =  reshape(Za[1:Int(sx*floor(size(Za,1)/sx))],
                    (sx,Int(floor(size(Za,1)/sx))))
        Y2 =  reshape(Zb[1:Int(sx*floor(size(Zb,1)/sx))],
                    (sx,Int(floor(size(Zb,1)/sx))))
        return Y1, Y2
    end
    function generalized(Za, Zb, sx)
        Na = floor(Int, size(Za,1)/sx)
        Nb = floor(Int, size(Zb,1)/sx)
        Y1 = var(reshape(Za[1:sx*Na],sx,Na)', corrected=false, dims=2)[:]
        Y2 = var(reshape(Zb[1:sx*Nb],sx,Nb)', corrected=false, dims=2)[:]
        return Y1, Y2
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
            r0 = Xt;
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
                r1 = r0.- (UpEnv(1:N) .+ LwEnv(1:N))./2  # r0.- (UpEnv.(1:N) .+ LwEnv.(1:N))./2  
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