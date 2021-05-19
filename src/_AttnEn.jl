module _AttnEn
export AttnEn
using StatsBase: Histogram, fit
    """
    # Av4, (Hxx,Hnn,Hxn,Hnx) = AttnEn(`Sig`) 

    Returns the attention entropy (`Av4`) calculated as the average of the
    sub-entropies (`Hxx`,`Hxn`,`Hnn`,`Hnx`) estimated from the data sequence
    (`Sig`) using a base-2 logarithm.

    # Av4, (Hxx, Hnn, Hxn, Hnx) = AttnEn(`Sig`, `Logx`, ___ )

    Returns the attention entropy (`Av4`) and the sub-entropies (`Hxx`,`Hnn`,`Hxn`,`Hnx`)
    from the data sequence (`Sig`) where,
    Hxx:    entropy of local-maxima intervals
    Hnn:    entropy of local minima intervals
    Hxn:    entropy of intervals between local maxima and subsequent minima
    Hnx:    entropy of intervals between local minima and subsequent maxima

    # Arguments:
    `Logx`  - Logarithm base, a positive scalar  
              (Enter 0 for natural logarithm)

    See also `EnofEn`, `SpecEn`, `XSpecEn`, `PermEn`, `MSEn`
  
    # References:
      [1] Jiawei Yang, et al.,
            "Classification of Interbeat Interval Time-series Using 
            Attention Entropy." 
            IEEE Transactions on Affective Computing 
            (2020)
  
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
    function AttnEn(Sig::AbstractArray{T,1} where T<:Real; Logx::Real=2)
        
    (Logx == 0) ? Logx = exp(1) : nothing
    N = size(Sig,1)
    (N > 10) ? nothing : error("Sig:   must be a numeric vector")
    (Logx>0) ? nothing : error("Logx:  must be a positive number > 0")

    Xmax = PkFind(Sig)
    Xmin = PkFind(-Sig)
    Txx = diff(Xmax)
    Tnn = diff(Xmin)
    Temp = diff(sort(vcat(Xmax, Xmin)))

    if isempty(Xmax) 
        error("No local maxima found!") 
    elseif isempty(Xmin)
        error("No local minima found!") 
    end
    
    (Xmax[1]<Xmin[1]) ? (Txn = Temp[1:2:end]; Tnx = Temp[2:2:end]) :
            (Txn = Temp[2:2:end]; Tnx = Temp[1:2:end])

    Edges = -0.5:N
    Pnx = fit(Histogram,Tnx,Edges).weights
    Pnn = fit(Histogram,Tnn,Edges).weights
    Pxx = fit(Histogram,Txx,Edges).weights
    Pxn = fit(Histogram,Txn,Edges).weights

    Pnx = Pnx[Pnx.!=0]/size(Tnx,1)
    Pxn = Pxn[Pxn.!=0]/size(Txn,1)
    Pnn = Pnn[Pnn.!=0]/size(Tnn,1)
    Pxx = Pxx[Pxx.!=0]/size(Txx,1)

    Hxx = -sum(Pxx.*(log.(Logx,Pxx)))
    Hxn = -sum(Pxn.*(log.(Logx,Pxn)))
    Hnx = -sum(Pnx.*(log.(Logx,Pnx)))
    Hnn = -sum(Pnn.*(log.(Logx,Pnn)))
    Av4 = (Hnn + Hxx + Hxn + Hnx)/4

    return Av4, (Hxx,Hnn,Hxn,Hnx)
    end

    function PkFind(X)
        Nx = size(X,1)
        Indx = zeros(Int,Nx);
        for n = 2:Nx-1
            if X[n-1]< X[n] > X[n+1]
                Indx[n] = n

            elseif X[n-1] < X[n] == X[n+1]
                k = 1
                while (n+k)<Nx && X[n] == X[n+k]
                    k +=1
                end
                if X[n] > X[n+k]
                    Indx[n] = n + floor((k-1)/2)
                end
            end
        end
        Indx = Indx[Indx.!==0]
    return Indx
    end

end