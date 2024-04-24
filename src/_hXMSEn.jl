module _hXMSEn
export hXMSEn
using Statistics: std, mean, median, var
using Plots
    """
        MSx, Sn, CI = hXMSEn(Sig1, Sig2, Mobj)

    Returns a vector of cross-entropy values (`MSx`) calculated at each node 
    in the hierarchical tree, the average cross-entropy value across all 
    nodes at each scale (`Sn`), and the complexity index (`CI`) of the hierarchical 
    tree (i.e. sum(`Sn`)) between the data sequences contained in `Sig1` and `Sig2` using
    the parameters specified by the multiscale object (`Mobj`) over 3 temporal
    scales (default).
    The entropy values in `MSx` are ordered from the root node (S.00) to the
    Nth subnode at scale T (S.TN): i.e. S.00, S.10, S.11, S.20, S.21, S.22,
    S.23, S.30, S.31, S.32, S.33, S.34, S.35, S.36, S.37, S.40, ... , S.TN.
    The average cross-entropy values in Sn are ordered in the same way, with the
    value of the root node given first: i.e. S0, S1, S2, ..., ST
     
        MSx, Sn, CI = hXMSEn(Sig1::AbstractVector{T} where T<:Real, Sig2::AbstractVector{T} where T<:Real, Mobj::NamedTuple; 
                                 Scales::Int=3, RadNew::Int=0, Plotx::Bool=false)
    
    Returns a vector of cross-entropy values (`MSx`) calculated at each node 
    in the hierarchical tree, the average cross-entropy value across all
    nodes at each scale (`Sn`), and the complexity index (`CI`) of the entire
    hierarchical tree between the data sequences contained in `Sig1` and `Sig2` using 
    the following name/value pair arguments:

    # Arguments:
    `Scales`   - Number of temporal scales, an integer > 1   (default: 3)
                 At each scale (T), entropy is estimated for 2^(T-1) nodes. \n
    `RadNew`   - Radius rescaling method, an integer in the range [1 4].
                 When the entropy specified by `Mobj` is `XSampEn` or `XApEn`, 
                 `RadNew` allows the radius threshold to be updated at each 
                 node in the tree. If a radius value is specified by `Mobj` (`r`),
                 this becomes the rescaling coefficient, otherwise it is set
                 to 0.2 (default). The value of `RadNew` specifies one of the 
                 following methods: \n
                 [1]    Pooled Standard Deviation          - r*std(Xt) \n
                 [2]    Pooled Variance                    - r*var(Xt) \n
                 [3]    Total Mean Absolute Deviation      - r*mean_ad(Xt) \n
                 [4]    Total Median Absolute Deviation    - r*med_ad(Xt) \n
    `Plotx`    - When `Plotx` == true, returns a plot of the average cross-entropy 
                 value at each time scale (i.e. the multiscale entropy curve)
                 and a hierarchical graph showing the entropy value of each node
                 in the hierarchical tree decomposition.  (default: false) \n
  
    # See also `MSobject`, `XMSEn`, `rXMSEn`, `cXMSEn`, `XSampEn`, `XApEn`, `hMSEn`
  
    # References:
        [1] Matthew W. Flood (2021), 
            "EntropyHub - An open source toolkit for entropic time series analysis"
            PLoS ONE 16(11):e0295448, 
            DOI:  10.1371/journal.pone.0259448
            https://www.EntropyHub.xyz
    
        [2]   Rui Yan, Zhuo Yang, and Tao Zhang,
            "Multiscale cross entropy: a novel algorithm for analyzing two
            time series." 
            5th International Conference on Natural Computation. 
            Vol. 1, pp: 411-413 IEEE, 2009.
  
        [3] Ying Jiang, C-K. Peng and Yuesheng Xu,
            "Hierarchical entropy analysis for biological signals."
            Journal of Computational and Applied Mathematics
            236.5 (2011): 728-742.
  
    """
    function hXMSEn(Sig1::AbstractVector{T} where T<:Real, Sig2::AbstractVector{T} where T<:Real, Mobj::NamedTuple; 
        Scales::Int=3, RadNew::Int=0, Plotx::Bool=false)

    (size(Sig1,1)>=10)  && (size(Sig2,1)>=10) ? nothing : error("Sig1/Sig2:   sequences must have >= 10 values")
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
            with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("XSampEn","XApEn"))) ? nothing :
    error("RadNew:  must be 0, or an integer in range [1 4] with entropy function `XSampEn` or `XApEn`")
    
    String(Symbol(Mobj.Func))=="XSampEn" ? Mobj = merge(Mobj,(Vcp=false,)) : nothing            

    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
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

    XX, YY, Na, Nb = Hierarchy(Sig1, Sig2, Scales)
    MSx = zeros(size(XX,1))
    for T in eachindex(XX[:,1]) # = 1:size(XX,1)
        print(" .")
        TempA = XX[T,1:Int(Na/(2^(floor(log2(T)))))]
        TempB = YY[T,1:Int(Nb/(2^(floor(log2(T)))))]
        RadNew > 0 ? Args = (Args..., r=Cx*Rnew(TempA, TempB)) : nothing      
        Temp2 = Mobj.Func(TempA, TempB; Args...)
        typeof(Temp2)<:Tuple ? MSx[T] = Temp2[1][end] : MSx[T] = Temp2[end]
    end

    Sn = zeros(Scales)
    for t = 1:Scales
        Sn[t] = mean(MSx[2^(t-1):(2^t)-1])
    end

    CI = sum(Sn)
    print("\n")
    if any(isnan.(MSx))
        println("Some entropy values may be undefined.")
    end

    if Plotx
        p1 = plot(1:Scales, Sn, c=RGB(8/255, 63/255, 77/255), lw=3)
        scatter!(1:Scales, Sn, markersize=6, c=RGB(1, 0, 1),
        xlabel = "Scale Factor", ylabel = "Entropy Value", 
        guidefont = font(12, "arial", RGB(7/255, 54/255, 66/255)),
        tickfontsize = 10, tickfontfamily="arial", legend=false,
        title = "Hierarchical Multiscale $(Mobj.Func) Entropy", 
        grid=false, xlim=(0.5,Scales+.5), ylim=(minimum(Sn)-.25,maximum(Sn)+.25),
        plot_titlefontsize=16, plot_titlefontcolor=RGB(7/255, 54/255, 66/255))

        N = 2^(Scales-1)
        x = zeros(2*N - 1)  
        x[1] = N        
        y = Scales.*(Scales .- floor.(Int, log2.(1:(2*N)-1)))
        for k = 1:2*N-1
            Q = floor(Int, log2(k))
            P = floor(Int,k/2)           
            if k>1              
                Bool(k%2) ? x[k] = x[P] + N/(2^Q) : x[k] = x[P] - N/(2^Q)
            end
        end
                    
        Edges = hcat(repeat(1:N-1,inner=2),2:2*N-1) 
        labx = [k for k in string.(round.(MSx,digits=3))]
        p2 = scatter(x,y,markersize=10*(MSx.-(minimum(MSx)).+1)./(abs(minimum(MSx))+1),
        c=RGB(1,0,1), legend=false, xticks=false, yticks=false)  #, txt=labx, markerfontsize=8)
        annotate!(x, y.+1, labx, 12-Scales)
        for k = 1:length(x)-1        
            plot!(x[Edges[k,:]],y[Edges[k,:]],c=RGB(8/255,63/255,77/255),
            lw=2.5, ylim=(Scales-1, maximum(y)+1), xlim=(0,2*N), grid=false)            
        end
        px = plot()
        pt = plot(p1, px, p2, layout = grid(3,1, heights=[0.2,0.1, 0.7]))
        display(pt)  

    end

    return MSx, Sn, CI
    end

    function Hierarchy(Za, Zb,sx)
        Na = Int(2^floor(log2(size(Za,1))))      
        Nb = Int(2^floor(log2(size(Zb,1))))
        if mod(log2(size(Za,1)),1) != 0
            @warn("Only first $(Na) samples were used in hierarchical decomposition.
                The last $(size(Za,1)-Na) samples of the data sequence were ignored.")
        end
        if mod(log2(size(Zb,1)),1) != 0
            @warn("Only first $(Nb) samples were used in hierarchical decomposition.
                The last $(size(Zb,1)-Nb) samples of the data sequence were ignored.")
        end

        if Na/(2^(sx-1)) < 8
           error("Data length ($(Na)) of Sig1 is too short to estimate entropy at the lowest subtree.
           Consider reducing the number of scales.") 
        elseif  Nb/(2^(sx-1)) < 8
            error("Data length ($(Nb)) of Sig2 is too short to estimate entropy at the lowest subtree.
            Consider reducing the number of scales.") 
        end
        
        Za = Za[1:Na,:]
        Zb = Zb[1:Nb,:]
        
        U1 = zeros((2^sx)-1,Na);  
        U2 = zeros((2^sx)-1,Nb);
        U1[1,:] = Za;  
        U2[1,:] = Zb;
        
        p=2
        for k = 1:sx-1
            for n = 1:2^(k-1)
                Temp = U1[2^(k-1)+n-1,:]       
                U1[p,1:Int(Na/2)]  = (Temp[1:2:end] + Temp[2:2:end])/2
                U1[p+1,1:Int(Na/2)]= (Temp[1:2:end] - Temp[2:2:end])/2
                p +=2
            end
        end

        p=2
        for k = 1:sx-1
            for n = 1:2^(k-1)
                Temp = U2[2^(k-1)+n-1,:]       
                U2[p,1:Int(Nb/2)]  = (Temp[1:2:end] + Temp[2:2:end])/2
                U2[p+1,1:Int(Nb/2)]= (Temp[1:2:end] - Temp[2:2:end])/2
                p +=2
            end
        end

        return U1, U2, Na, Nb
    end
end

"""
MM = zeros(2*Scales - 1,T)
mid = Int(ceil(T/2))
for t = 1:Scales
    MM[2*t -1, mid-(2^(t-1))+1:2:mid-1+(2^(t-1))] = MSx[2^(t-1):(2^t)-1]
end  

b = bar3(MM)

zlim([min(MSx(MSx>0)) max(MSx)]); 
cmap = colormap('cool');  colormap([1 1 1; cmap])
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
axis square, view([30 55])        
yticks(1:2:(2*T -1));  yticklabels(0:T-1)
xticks(''); %     xticks(1:2:(2*p.Results.Scales - 1))    
xlabel('Subtree Nodes','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
ylabel('Scale Factor','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
zlabel('Entropy','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
title(sprintf('Hierarchical Multiscale (%s) Entropy',func2str(Y{1})),...
    'FontSize',16,'FontWeight','bold','Color',[7 54 66]/255)  
"""

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