module _RangEn
export RangEn
using Statistics: mean
    """
        Rangx, A, B = RangEn(Sig) 
        
    Returns the range entropy estimate (`Rangx`) and the number of matched state
    vectors (`m: B`, `m+1: A`) estimated from the data sequence (`Sig`)
    using the sample entropy algorithm and the following default parameters: 
    embedding dimension = 2, time delay = 1, radius threshold = 0.2, logarithm = natural.    

        Rangx, A, B = RangEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, r::Real=0.2, Methodx::String="SampEn", Logx::Real=exp(1))
        
    Returns the range entropy estimates (`Rangx`) for dimensions = `m`
    estimated for the data sequence (`Sig`) using the specified keyword arguments:
       
    # Arguments:    
    `m`         - Embedding Dimension, a positive integer\n
    `tau`       - Time Delay, a positive integer\n
    `r`         - Radius Distance Threshold, a positive value between 0 and 1\n
    `Methodx`   - Base entropy method, either 'SampEn' [default] or 'ApEn'\n
    `Logx`      - Logarithm base, a positive scalar  \n

    # See also `ApEn`, `SampEn`, `FuzzEn`,  `MSEn`

    # References:
        [1] Omidvarnia, Amir, et al.
            "Range entropy: A bridge between signal complexity and self-similarity"
            Entropy 
            20.12 (2018): 962.
            
        [2] Joshua S Richman and J. Randall Moorman. 
            "Physiological time-series analysis using approximate entropy
            and sample entropy." 
            American Journal of Physiology-Heart and Circulatory Physiology 
            2000

    """
    function RangEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1, 
                 r::Real=0.2, Methodx::String="SampEn", Logx::Real=exp(1))
    
    N = length(Sig)
    (N>10) ? nothing : error("Sig:   must be a numeric vector")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (r>=0) && (r<=1) ? nothing : error("r:     must be a scalar must be a value between 0 and 1")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")
    (lowercase(Methodx) in ["sampen", "apen"]) ? nothing : error("Methodx must be either 'ApEn' or 'SampEn'")  
                 
    if lowercase(Methodx) == "sampen"         
        Nx = N - m*tau
        Sx = zeros(Nx,m+1)
        for k = 1:m+1        
            Sx[:,k] = Sig[1 + (k-1)*tau:Nx + (k-1)*tau]
        end
       
        A = zeros(Int, Nx)
        B = zeros(Int, Nx)
        for k = 1:(Nx - 1)          
            Dxy = abs.(repeat(Sx[k,1:end-1],1,Nx-k) .- Sx[k+1:end,1:end-1]')'            
            Mx = maximum(Dxy,dims=2)
            Mn = minimum(Dxy,dims=2)
            RR = (Mx .- Mn)./(Mx .+ Mn) .<= r
            B[k] = sum(RR)

            if B[k]>0  
                Dxy2 = abs.(repeat(Sx[k,:],1,B[k]) .- Sx[k+1:end,:][RR[:],:]')'#Sx[k+1:end,:][all.(eachrow(RR)),:]')'              
                Mx = maximum(Dxy2,dims=2)   
                Mn = minimum(Dxy2,dims=2)
                RR2 = (Mx .- Mn)./(Mx .+ Mn) .<= r                        
                A[k] = sum(RR2)
            end
        end
                    
        Rangx = -log(sum(A)/sum(B))/log(Logx) 
       
        return Rangx, A, B
        
    elseif lowercase(Methodx) == "apen"
        Nx = N - (m-1)*tau
        Sx = zeros(Nx,m)
        for k = 1:m    
            Sx[:,k] = Sig[1 + (k-1)*tau:Nx + (k-1)*tau]
        end               
        Sx = hcat(Sx, vcat(Sig[m*tau + 1:end],zeros(tau)))    

        B = zeros(Int, Nx)
        A = zeros(Int, Nx-tau)
        for k in 1:Nx                    
            Dxy = abs.(repeat(Sx[k,1:end-1],1,Nx)' .- Sx[:,1:end-1])             
            Mx = maximum(Dxy,dims=2)
            Mn = minimum(Dxy,dims=2)
            RR = (Mx .- Mn)./(Mx .+ Mn) .<= r             
            B[k] = sum(RR)
                     
            if k <= (Nx - tau)
                RR[end-tau+1:end] .= false
                Dxy2 = abs.(repeat(Sx[k,:],1,sum(RR))' .-  Sx[RR[:],:]) #Sx[all.(eachrow(RR)),:])            
                Mx2 = maximum(Dxy2,dims=2)
                Mn2 = minimum(Dxy2,dims=2)
                RR2 = (Mx2 .- Mn2)./(Mx2 .+ Mn2) .<= r                     
                A[k] = sum(RR2)
            end
        end

        Ax = mean(log.(A./(Nx-tau))/log(Logx))
        Bx = mean(log.(B./Nx)/log(Logx))
        Ap = Bx - Ax
       
        return Ap, A, B  
    
    else
        error("Methodx must be either 'ApEn' or 'SampEn'")
    end

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