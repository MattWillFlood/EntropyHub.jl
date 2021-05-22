module _hMSEn
export hMSEn
using Statistics: std, mean, median, var
using Plots        
      
    """
        MSx, Sn, CI = hMSEn(Sig, Mobj) 

    Returns a vector of entropy values (`MSx`) calculated at each node in the
    hierarchical tree, the average entropy value across all nodes at each 
    scale (Sn), and the complexity index (`CI`) of the hierarchical tree
    (i.e. sum(Sn)) for the data sequence (`Sig`) using the parameters specified
    by the multiscale object (`Mobj`) over 3 temporal scales (default).
    The entropy values in MSx are ordered from the root node (S_00) to the
    Nth subnode at scale T (S_TN): i.e. S_00, S_10, S_11, S_20, S_21, S_22,
    S_23, S_30, S_31, S_32, S_33, S_34, S_35, S_36, S_37, S_40, ... , S_TN.
    The average entropy values in Sn are ordered in the same way, with the
    value of the root node given first: i.e. S0, S1, S2, ..., ST

        MSx, Sn, CI = hMSEn(Sig::AbstractArray{T,1} where T<:Real, Mobj::NamedTuple; 
                                Scales::Int=3, RadNew::Int=0, Plotx::Bool=false)

    Returns a vector of entropy values (`MSx`) calculated at each node in the
    hierarchical tree, the average entropy value across all nodes at each 
    scale (Sn), and the complexity index (CI) of the entire hierarchical
    tree for the data sequence (`Sig`) using the following 'keyword' arguments:

    # Arguments:
    `Scales`   - Number of temporal scales, an integer > 1   (default: 3)
                 At each scale (T), entropy is estimated for 2^(T-1) nodes.
    `RadNew`   - Radius rescaling method, an integer in the range [1 4].
                 When the entropy specified by `Mobj` is `SampEn` or `ApEn`, 
                 `RadNew` allows the radius threshold to be updated at each 
                 node in the tree. If a radius value is specified by `Mobj` (`r`),
                 this becomes the rescaling coefficient, otherwise it is set
                 to 0.2 (default). The value of `RadNew` specifies one of the 
                 following methods:
                 [1]    Standard Deviation          - r*std(Xt)
                 [2]    Variance                    - r*var(Xt)
                 [3]    Mean Absolute Deviation     - r*mean_ad(Xt)
                 [4]    Median Absolute Deviation   - r*med_ad(Xt,1)
    `Plotx`    - When `Plotx` == true, returns a plot of the average entropy 
                 value at each time scale (i.e. the multiscale entropy curve)
                 and a hierarchical graph showing the entropy value of each node
                 in the hierarchical tree decomposition.  (default: false)
  
    # See also `MSobject`, `MSEn`, `cMSEn`, `rMSEn`, `SampEn`, `ApEn`, `XMSEn`

    # References:
        [1] Ying Jiang, C-K. Peng and Yuesheng Xu,
          "Hierarchical entropy analysis for biological signals."
          Journal of Computational and Applied Mathematics
          236.5 (2011): 728-742.
    """
    function hMSEn(Sig::AbstractArray{T,1} where T<:Real, Mobj::NamedTuple; 
        Scales::Int=3, RadNew::Int=0, Plotx::Bool=false)

    (size(Sig,1)>10) ? nothing : error("Sig:   must be a numeric vector" )
    (length(Mobj) >= 1) ? nothing :  error("Mobj:    must be a multiscale entropy object created 
            with the function EntropyHub.MSobject")
    (Scales>1) ? nothing : error("Scales:     must be an integer > 1")
    (RadNew==0 || (RadNew in 1:4 && String(Symbol(Mobj.Func)) in ("SampEn","ApEn"))) ? nothing :
    error("RadNew:  must be 0, or an integer in range [1 4] with entropy function `SampEn` or `ApEn`")
        
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

    XX, N = Hierarchy(Sig, Scales)
    MSx = zeros(size(XX,1))
    Args = NamedTuple{keys(Mobj)[2:end]}(Mobj)
    for T = 1:size(XX,1)
        print(". ")
        Temp = XX[T,1:Int(N/(2^(floor(log2(T)))))]
        RadNew > 0 ? Args = (Args..., r=Cx*Rnew(Temp[:])) : nothing      
        Temp2 = Mobj.Func(Temp[:]; Args...)
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
        pt = plot(p1, px, p2, layout = grid(3,1, heights=[0.25,.1, 0.8]))

        display(pt)

    end

    return MSx, Sn, CI
    end

    function Hierarchy(Z,sx)
        N = Int(2^floor(log2(size(Z,1))))
        if mod(log2(size(Z,1)),1) != 0
            @warn("Only first $(N) samples were used in hierarchical decomposition.
                The last $(length(Z)-N) samples of the data sequence were ignored.")
        end
        if N/(2^(sx-1)) < 8
           error("Data length ($(N)) is too short to estimate entropy at the lowest subtree.
           Consider reducing the number of scales.") 
        end
        
        Z = Z[1:N]
        U = zeros((2^sx)-1,N)
        U[1,:] = Z;  p=2
        for k = 1:sx-1
            for n = 1:2^(k-1)
                Temp = U[2^(k-1)+n-1,:]       
                U[p,1:Int(N/2)]  = (Temp[1:2:end] + Temp[2:2:end])/2
                U[p+1,1:Int(N/2)]= (Temp[1:2:end] - Temp[2:2:end])/2
                p=p+2
            end
        end
        return U, N
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