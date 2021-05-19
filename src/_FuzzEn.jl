module _FuzzEn
export FuzzEn
using Statistics: mean, std
using LinearAlgebra: UpperTriangular, I, inv
    """
    #    Fuzz, Ps1, Ps2 = FuzzEn(`Sig`) 

    Returns the fuzzy entropy estimates `Fuzz` and the average fuzzy distances
    (`m`:Ps1, `m+1`:Ps2) for `m` = [1,2] estimated from the data sequence `Sig`
    using the default parameters: embedding dimension = 2, time delay = 1, 
    fuzzy function (`Fx`) = "default", fuzzy function parameters (`r`) = [0.2, 2],
    logarithm = natural

    #    Fuzz, Ps1, Ps2 = FuzzEn(`Sig`, 'keyword' = value, ...)

    Returns the fuzzy entropy estimates `Fuzz` for dimensions = [1,...,`m`]
    estimated for the data sequence `Sig` using the specified keyword arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer  [default: 2]\n
    `tau`   - Time Delay, a positive integer        [default: 1]\n
    `Fx`    - Fuzzy function name, one of the following:
                {"sigmoid","modsampen","default","gudermannian","linear"}\n
    `r`     - Fuzzy function parameters, a 1 element scalar or a 2 element
              tuple of positive values. The `r` parameters for each fuzzy
              function are defined as follows:      [default: [.2 2]]\n
                sigmoid:      r(1) = divisor of the exponential argument
                              r(2) = value subtracted from argument (pre-division)
                modsampen:    r(1) = divisor of the exponential argument
                              r(2) = value subtracted from argument (pre-division)
                default:      r(1) = divisor of the exponential argument
                              r(2) = argument exponent (pre-division)
                gudermannian: r  = a scalar whose value is the numerator of
                                    argument to gudermannian function:
                                    GD(x) = atan(tanh(r/x))
                linear:       r  = an integer value. When r = 0, the
                                    argument of the exponential function is 
                                    normalised between [0 1]. When r = 1,
                                    the minimuum value of the exponential 
                                    argument is set to 0.   \n                          
    `Logx`  - Logarithm base, a positive scalar  [default: natural]

    ## For further information on keyword arguments, see the EntropyHub guide.

    # See also `SampEn`, `ApEn`, `PermEn`, `DispEn`, `XFuzzEn`, `FuzzEn2D`, `MSEn`
  
    # References:
      [1] Weiting Chen, et al.
            "Characterization of surface EMG signal based on fuzzy entropy."
            IEEE Transactions on neural systems and rehabilitation engineering
            15.2 (2007): 266-272.
  
      [2] Hong-Bo Xie, Wei-Xing He, and Hui Liu
            "Measuring time series regularity using nonlinear
            similarity-based sample entropy."
            Physics Letters A
            372.48 (2008): 7140-7146.
            
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
    function FuzzEn(Sig::AbstractArray{T,1} where T<:Real; m::Int=2, tau::Int=1,
         r::Union{Real,Tuple{Real,Real}}=(.2,2), Fx::String="default", Logx::Real=exp(1))

    N = size(Sig)[1]
    (N>10) ? nothing : error("Sig:   must be a numeric vector")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (minimum(r)>=0 && length(r)<=2) ? 
        nothing : error("r:  must be a scalar or 2 element tuple of positive values")
    (lowercase(Fx) in ["default","sigmoid","modsampen","gudermannian","linear"]) ?
        nothing : error("Fx:    must be one of the following strings -
        'default', 'sigmoid', 'modsampen', 'gudermannian', 'linear'")
    (Logx>0) ? nothing : error("Logx:   must be a positive number > 0")
            
    if length(r) == 2 && lowercase(Fx)=="linear"
        r = 0;
        print("Multiple values for r entered. Default value (0) used.\n") 
    elseif length(r) == 2 && lowercase(Fx)=="gudermannian"
        r = r[1]
        print("Multiple values for r entered. First value used.\n") 
    end

    m += 1      
    Fun = getfield(_FuzzEn,Symbol(lowercase(Fx)))
    Sx = zeros(N,m)    
    for k = 1:m
        Sx[1:N-(k-1)*tau,k] = Sig[1+(k-1)*tau:end]
    end
    Ps1 = zeros(m)
    Ps2 = zeros(m-1)
    Ps1[1] = .5
    
    for k = 2:m
        N1 = N - k*tau;    N2 = N - (k-1)*tau;
        T2 = Sx[1:N2,1:k] .- mean(Sx[1:N2,1:k],dims=2)
        d2 = zeros(N2-1,N2-1)
        
        for p = 1:N2-1
            Mu2 = maximum(abs.(repeat(transpose(T2[p,:]),N2-p) - T2[p+1:end,:]),dims=2)
            d2[p,p:end] = Fun(Mu2,r)
        end
        d1 = d2[1:N1-1,1:N1-1]
        Ps1[k] = sum(d1)/(N1*(N1-1))
        Ps2[k-1] = sum(d2)/(N2*(N2-1))
    end

    Fuzz = (log.(Logx, Ps1[1:end-1]) - log.(Logx,Ps2))

    return Fuzz, Ps1, Ps2
    end


    function sigmoid(x,r)
        if length(r) == 1
            error("When Fx = 'Sigmoid', r must be a two-element tuple.")
        end
        y = inv.(1 .+ exp.((x.-r[2])/r[1]))
        return y    
    end

    function modsampen(x,r)
        if length(r) == 1
            error("When Fx = 'Modsampen', r must be a two-element tuple.")
        end
        y = inv.(1 .+ exp.((x.-r[2])/r[1]))
        return y    
    end

    function default(x,r)   
        if length(r) == 1
            error("When Fx = 'Default', r must be a two-element tuple.")
        end
        y = exp.(-(x.^r[2])/r[1])
        return y         
    end

    function gudermannian(x,r)
        if r <= 0
            error("When Fx = 'Gudermannian', r must be a scalar > 0.")
        end
        y = atan.(tanh.(r[1]./x))    
        y ./= maximum(y)    
        return y    
    end

    function linear(x,r)
        if r == 0 && length(x)>1
            y = exp.(-(x .- minimum(x))/(maximum(x)-minimum(x)))
        elseif r == 1
            y = exp.(-(x .- minimum(x)))
        elseif r == 0 && length(x)==1
            y = [0]
        else
            error("When Fx = 'Linear', r must be 0 or 1.")
        end
        return y
    end

end