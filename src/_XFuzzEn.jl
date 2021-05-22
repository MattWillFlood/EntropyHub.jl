module _XFuzzEn
export XFuzzEn
using Statistics: mean, std
    """
        XFuzz, Ps1, Ps2 = XFuzzEn(Sig) 

    Returns the cross-fuzzy entropy estimates (`XFuzz`) and the average fuzzy
    distances (m:Ps1, m+1:Ps2) for m = [1,2] estimated for the data sequences
    contained in 'Sig', using the default parameters: embedding dimension = 2,
    time delay = 1, fuzzy function (Fx) = 'default', 
    fuzzy function parameters (r) = [0.2, 2], logarithm = natural

        XFuzz, Ps1, Ps2 = XFuzzEn(Sig::AbstractArray{T,2} where T<:Real; m::Int=2, tau::Int=1, r::Union{Real,Tuple{Real,Real}}=(.2,2), Fx::String="default", Logx::Real=exp(1))

    Returns the cross-fuzzy entropy estimates (`XFuzz`) for dimensions = [1,...,m]
    estimated for the data sequences in `Sig` using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Embedding Dimension, a positive integer   [default: 2]    \n
    `tau`   - Time Delay, a positive integer            [default: 1]    \n
    `Fx`    - Fuzzy funtion name, one of the following: 
              {`default`,`sigmoid`,`modsampen`,`gudermannian`,`linear`}\n
    `r`     - Fuzzy function parameters, a scalar or a 2 element tuple 
              of positive values. The `r` parameters for each fuzzy
              function are defined as follows:\n
              sigmoid:       r(1) = divisor of the exponential argument
                             r(2) = value subtracted from argument (pre-division)
              modsampen:     r(1) = divisor of the exponential argument
                             r(2) = value subtracted from argument (pre-division)
              default:       r(1) = divisor of the exponential argument
                             r(2) = argument exponent (pre-division)
              gudermannian:  r    = a scalar whose value is the numerator of
                                    argument to gudermannian function:
                                    GD(x) = atan(tanh(`r`/x)).
              linear:        r    = an integer value. When `r` = 0, the
                                    argument of the exponential function is 
                                    normalised between [0 1]. When `r` = 1,
                                    the minimuum value of the exponential 
                                    argument is set to 0.     \n                     
    `Logx`  - Logarithm base, a positive scalar  \n

    For further information on the 'keyword' arguments, see the EntropyHub guide.

    # See also `FuzzEn`, `XSampEn`, `XApEn`, `FuzzEn2D`, `XMSEn`, `MSEn`

    # References:
        [1] Hong-Bo Xie, et al.,
            "Cross-fuzzy entropy: A new method to test pattern synchrony of
            bivariate time series." 
            Information Sciences 
            180.9 (2010): 1715-1724.
  
    """
    function XFuzzEn(Sig::AbstractArray{T,2} where T<:Real; m::Int=2, tau::Int=1,
         r::Union{Real,Tuple{Real,Real}}=(.2,2), Fx::String="default", Logx::Real=exp(1))

    (size(Sig,2) > size(Sig,1)) ? Sig = transpose(Sig) : nothing
    N = size(Sig,1)
    (N>10 && size(Sig,2)==2) ? nothing :  error("Sig:   must be a 2-columns matrix")
    (m > 0) ? nothing : error("m:     must be an integer > 0")
    (tau>0) ? nothing : error("tau:   must be an integer > 0")
    (minimum(r)>=0 && length(r)<=2) ? nothing : 
        error("r:  must be a scalar or 2 element vector of positive values")
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

    S1 = Sig[:,1]; S2 = Sig[:,2]
    m += 1      
    Fun = getfield(_XFuzzEn,Symbol(lowercase(Fx)))
    Sx1 = zeros(N,m)
    Sx2 = zeros(N,m)  
    for k = 1:m
        Sx1[1:N-(k-1)*tau,k] = S1[1 + (k-1)*tau:N]
        Sx2[1:N-(k-1)*tau,k] = S2[1 + (k-1)*tau:N]
    end

    Ps1 = zeros(m)
    Ps2 = zeros(m-1)
    Ps1[1] = 1
    for k = 2:m
        N1 = N - k*tau
        N2 = N - (k-1)*tau
        A = Sx1[1:N2,1:k] .- mean(Sx1[1:N2,1:k],dims=2)
        B = Sx2[1:N2,1:k] .- mean(Sx2[1:N2,1:k],dims=2)
        d2 = zeros(N2,N2)
        for p = 1:N2
            Mu2 = maximum(abs.(transpose(A[p,:]) .- B),dims=2)
            d2[p,:] = Fun(Mu2,r)
        end    
        Ps1[k] = mean(d2[1:N1,1:N1])
        Ps2[k-1] = mean(d2)
    end

    XFuzz = log.(Logx, Ps1[1:end-1]) .- log.(Logx, Ps2)    
    return XFuzz, Ps1, Ps2
    end

    function sigmoid(x,r)
        if length(r) == 1
            error("When Fx = 'Sigmoid', r must be a two-element vector.")
        end
        y = inv.(1 .+ exp.((x.-r[2])/r[1]))
        return y    
    end

    function modsampen(x,r)
        if length(r) == 1
            error("When Fx = 'Modsampen', r must be a two-element vector.")
        end
        y = inv.(1 .+ exp.((x.-r[2])/r[1]))
        return y    
    end

    function default(x,r)   
        if length(r) == 1
            error("When Fx = 'Default', r must be a two-element vector.")
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