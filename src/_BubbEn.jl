module _BubbEn
export BubbEn
using GroupSlices
    """
        Bubb, H = BubbEn(Sig) 

    Returns the bubble entropy (`Bubb`) and the conditional Rényi entropy (`H`)
    estimates of dimension m = 2 from the data sequence (`Sig`) using 
    the default parameters: 
    embedding dimension = 2, time delay = 1, logarithm = natural 

        Bubb, H = BubbEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, Logx::Real=exp(1))

    Returns the bubble entropy (`Bubb`) estimate of the data sequence (`Sig`)  
    using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, an integer > 1   \n
              BubbEn returns estimates for each dimension [2,...,m]
    `tau`   - Time Delay, a positive integer    \n
    `Logx`  - Logarithm base, a positive scalar \n

    # See also `PhasEn`, `MSEn`

    # References:
        [1] George Manis, M.D. Aktaruzzaman and Roberto Sassi,
            "Bubble entropy: An entropy almost free of parameters."
            IEEE Transactions on Biomedical Engineering
            64.11 (2017): 2711-2718.

    """
    function BubbEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, Logx::Real=exp(1))

    N = size(Sig,1)
    (N > 10) ? nothing :  error("Sig:   must be a numeric vector")
    (m > 1) ? nothing :   error("m:     must be an integer > 1")
    (tau >0) ? nothing :  error("tau:   must be an integer > 0")
    (Logx>0) ? nothing :  error("Logx:     must be a positive scalar > 0")

    Sx = zeros(N,m+1)
    H = zeros(m+1)
    Sx[:,1] = Sig
    for k = 2:m+1    
        Sx[1:N-(k-1)*tau,k] = Sig[1+(k-1)*tau:N]
        Swapx  = BubbSort(Sx[1:N-(k-1)*tau,1:k])

        Locs = getindex.(indexin(Swapx, unique(Swapx)))
        Temp = unique(Locs)
        
        p = zeros(size(Temp,1))
        for n in Temp 
            p[n] = sum(Locs.==n);
        end  
        p ./= (N-(k-1)*tau)
        H[k] = -log(Logx,  sum(p.^2))
        
        if round(sum(p),digits=6) != 1
            @warn("Potential error in detected swap number")
        end    
    end
        
    Bubb = diff(H)./log.(Logx,  (2:m+1)./(0:m-1))
    Bubb = Bubb[2:end]

    return Bubb, H
    end


    function BubbSort(Data)  
        x,N2 = size(Data)
        swaps = zeros(Int, x)
        for y = 1:x
            t = 1
            while t <= N2-1
                for kk = 1:N2-t
                    if Data[y,kk] > Data[y,kk+1]
                        temp = Data[y,kk]
                        Data[y,kk] = Data[y,kk+1]
                        Data[y,kk+1] = temp
                        swaps[y] = swaps[y] + 1
                    end
                end
                t = t + 1;
            end
        end
        bsorted = Data;
        return swaps # bsorted
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