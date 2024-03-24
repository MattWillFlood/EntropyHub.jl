module _FuzzEn2D
export FuzzEn2D
using Statistics: mean, std
    """
        Fuzz2D = FuzzEn2D(Mat) 

    Returns the bidimensional fuzzy entropy estimate (`Fuzz2D`) estimated for 
    the data matrix (`Mat`) using the default parameters: time delay = 1,
    fuzzy function (Fx) = 'default', fuzzy function parameters (r) = [0.2, 2],
    logarithm = natural, template matrix size = [floor(H/10) floor(W/10)], 
    (where H and W represent the height and width of the data matrix 'Mat') \n
    ** The minimum dimension size of Mat must be > 10.**

        Fuzz2D = FuzzEn2D(Mat::AbstractArray{T,2} where T<:Real; m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10), 
                            tau::Int=1, r::Union{Real,Tuple{Real,Real}}=(.2*std(Mat, corrected=false),2), 
                                Fx::String="default", Logx::Real=exp(1), Lock::Bool=true)

    Returns the bidimensional fuzzy entropy (`Fuzz2D`) estimates for the data
    matrix (`Mat`) using the specified 'keyword' arguments:

    # Arguments:
    `m`     - Template submatrix dimensions, an integer scalar (i.e. the same
              height and width) or a two-element vector of integers 
              [height, width] with a minimum value > 1.
              (default: [floor(H/10) floor(W/10)])  \n
    `tau`   - Time Delay, a positive integer   (default: 1)  \n
    `Fx`    - Fuzzy function name, one of the following: 
              {`"sigmoid", "modsampen", "default", "gudermannian",`
              `"bell", "triangular", "trapezoidal1", "trapezoidal2",`
              `"z_shaped", "gaussian", "constgaussian"`}\n 
    `r`     - Fuzzy function parameters, a 1 element scalar or a 2 element
              vector of positive values. The 'r' parameters for each fuzzy
              function are defined as follows:\n
              sigmoid:        r(1) = divisor of the exponential argument
                              r(2) = value subtracted from argument (pre-division)
              modsampen:      r(1) = divisor of the exponential argument
                              r(2) = value subtracted from argument (pre-division)
              default:        r (1) = divisor of the exponential argument
                              r(2) = argument exponent (pre-division)
              gudermannian:   r  = a scalar whose value is the numerator of
                                   argument to gudermannian function:
                                   GD(x) = atan(tanh(r/x))
              triangular:     r = a scalar whose value is the threshold (corner point) of the triangular function.
              trapezoidal1:   r = a scalar whose value corresponds to the upper (2r) and lower (r) corner points of the trapezoid.
              trapezoidal2:   r(1) = a value corresponding to the upper corner point of the trapezoid.
                              r(2) = a value corresponding to the lower corner point of the trapezoid.
              z_shaped:       r = a scalar whose value corresponds to the upper (2r) and lower (r) corner points of the z-shape.
              bell:           r(1) = divisor of the distance value
                              r(2) = exponent of generalized bell-shaped function
              gaussian:       r = a scalar whose value scales the slope of the Gaussian curve.
              constgaussian:  r = a scalar whose value defines the lower threshod and shape of the Gaussian curve.                    
              [DEPRICATED] linear:       r  = an integer value. When r = 0, the
                                argument of the exponential function is 
                                normalised between [0 1]. When r = 1,
                                the minimuum value of the exponential 
                                argument is set to 0.   \n    
    `Logx`  - Logarithm base, a positive scalar    (default: natural)\n
    `Lock`  - By default, FuzzEn2D only permits matrices with a maximum
              size of 128 x 128 to prevent memory errors when storing data on RAM. 
              e.g. For Mat = [200 x 200], m = 3, and tau = 1, FuzzEn2D 
              creates a vector of 753049836 elements. To enable matrices
              greater than [128 x 128] elements, set `Lock` to false.
              (default: true)\n
              ` WARNING: unlocking the permitted matrix size may cause
              your Julia IDE to crash.` \n

    # See also `SampEn2D`, `FuzzEn`, `XFuzzEn`

    # References:
        [1] Luiz Fernando Segato Dos Santos, et al.,
            "Multidimensional and fuzzy sample entropy (SampEnMF) for
            quantifying H&E histological images of colorectal cancer."
            Computers in biology and medicine 
            103 (2018): 148-160.

        [2] Mirvana Hilal and Anne Humeau-Heurtier,
            "Bidimensional fuzzy entropy: Principle analysis and biomedical
            applications."
            41st Annual International Conference of the IEEE (EMBC) Society
            2019.

        [3] Hamed Azami, et al.
            "Fuzzy Entropy Metrics for the Analysis of Biomedical Signals: 
            Assessment and Comparison"
            IEEE Access
            7 (2019): 104833-104847

    """
    function FuzzEn2D(Mat::AbstractArray{T,2} where T<:Real; 
        m::Union{Int,Tuple{Int,Int}}=floor.(Int, size(Mat)./10), 
        tau::Int=1, r::Union{Real,Tuple{Real,Real}}=(.2*std(Mat,corrected=false),2), 
        Fx::String="default", Logx::Real=exp(1), Lock::Bool=true)
        
    NL, NW = size(Mat)
    ((NL > 128 || NW > 128) && Lock) ? 
        error("To prevent memory errors, matrix width & length must have <= 128 elements.
            To estimate FuzzEn2D  for the current matrix ($NL,$NW) change Lock to 'false'.
            Caution: unlocking the safe matrix size may cause the Julia IDE to crash.") :
            nothing

    length(m)==1 ?  (mL = m; mW = m) : (mL = m[1];  mW = m[2])
    (NL > 10 && NW > 10) ? nothing : 
    error("Number of rows and columns in Mat must be > 10")
    (minimum(m)>1) ? nothing : error("m:  must be an integer > 1, or a 2 element integer tuple > 1")
    (tau > 0) ? nothing : error("tau:   must be an integer > 0")
    (minimum(r)>=0) ? nothing : error("r:     must be a positive value")
        #error("r:     must be 2 element tuple of positive values")
    (lowercase(Fx) in ["default","sigmoid","modsampen","gudermannian","bell", "z_shaped",
        "triangular", "trapezoidal1","trapezoidal2","gaussian","constgaussian"]) ?
        nothing : error("Fx:    must be one of the following strings -
        'default', 'sigmoid', 'modsampen', 'gudermannian', 'bell', 'z_shaped',
        'triangular', 'trapezoidal1','trapezoidal2','gaussian','constgaussian'")
    (Logx>0) ? nothing :  error("Logx:     must be a positive number > 0")

    if length(r) == 2 && lowercase(Fx)=="linear"
        r = 0;
        print("Multiple values for r entered. Default value (0) used.\n") 
    elseif length(r) == 2 && lowercase(Fx)=="gudermannian"
        r = r[1]
        print("Multiple values for r entered. First value used.\n") 
    end

    Fun = getfield(_FuzzEn2D,Symbol(lowercase(Fx)))
    NL = NL - mL*tau
    NW = NW - mW*tau
    X1 = zeros(mL,mW,NL*NW)
    X2 = zeros(mL+1,mW+1,NL*NW)
    p = 0;
    for k = 1:NL
        for n = 1:NW
            p += 1
            Temp2 = Mat[k:tau:(mL)*tau+k,n:tau:(mW)*tau+n]
            Temp1 = Temp2[1:end-1,1:end-1]
            X1[:,:,p] = Temp1 .- mean(Temp1)
            X2[:,:,p] = Temp2 .- mean(Temp2)
        end
    end
    p = size(X1,3)
    p != NL*NW ? @warn("Potential error with submatrix division.") : nothing
    Ny = p*(p-1)/2;
    Ny > 300000000 ? @warn("Number of pairwise distance calculations is $Ny") : nothing
    Y1 = zeros(p-1)
    Y2 = zeros(p-1)
    for k = 1:p-1
        Temp1 = maximum(abs.(X1[:,:,k+1:end] .- X1[:,:,k]),dims=(1,2))
        Y1[k] = sum(Fun(Temp1[:], r))
        Temp2 = maximum(abs.(X2[:,:,k+1:end] .- X2[:,:,k]),dims=(1,2))
        Y2[k] = sum(Fun(Temp2[:], r))
    end

    Fuzz2D = -log(Logx, sum(Y2)/sum(Y1))
    return Fuzz2D
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

    """
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
    """

    function triangular(x,r)
        length(r)==1 ? nothing : error("When Fx = 'Triangular', r must be a scalar > 0.")
        y = 1 .- (x./r)
        y[x .> r] .= 0
        return y
    end

    function trapezoidal1(x, r)
        length(r)==1 ? nothing : ("When Fx = 'Trapezoidal1', r must be a scalar > 0.")
        y = zeros(length(x))
        y[x .<= r*2] = 2 .- (x[x .<= r*2]./r)
        y[x .<= r] .= 1
        return y
    end

    function trapezoidal2(x, r)
        (r isa Tuple) && (length(r)==2) ? nothing : error("When Fx = 'Trapezoidal2', r must be a two-element tuple.")
        y = zeros(length(x))
        y[x .<= maximum(r)] = 1 .- (x[x .<= maximum(r)] .- minimum(r))./(maximum(r)-minimum(r))
        y[x .<= minimum(r)] .= 1
        return y
    end

    function z_shaped(x, r)
        length(r)==1 ? nothing : error("When Fx = 'Z_shaped', r must be a scalar > 0.")
        y = zeros(length(x))
        y[x .<= 2*r] .= 2*(((x[x .<= 2*r] .- 2*r)./r).^2)
        y[x .<= 1.5*r] .= 1 .- (2*(((x[x .<= 1.5*r] .- r)/r).^2))
        y[x .<= r] .= 1
        return y
    end

    function bell(x, r)
        (r isa Tuple) && length(r)==2 ? nothing : error("When Fx = 'Bell', r must be a two-element tuple.")
        y = inv.(1 .+ abs.(x./r[1]).^(2*r[2]))
        return y
    end

    function gaussian(x, r)
        length(r)==1 ? nothing : error("When Fx = 'Gaussian', r must be a scalar > 0.")
        y = exp.(-((x.^2)./(2*(r.^2))))
        return y
    end

    function constgaussian(x, r)
        length(r)==1 ? nothing : error("When Fx = 'ConstGaussian', r must be a scalar > 0.")
        y = ones(length(x))
        y[x .> r] = exp.(-log(2)*((x[x .> r] .- r)./r).^2)
        return y
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