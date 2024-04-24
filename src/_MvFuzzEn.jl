module _MvFuzzEn
export MvFuzzEn
using Statistics: std, mean
"""
    MFuzz, B0, Bt, B1 = MvFuzzEn(Data) 

 Returns the multivariate fuzzy entropy estimate (`MFuzz`) and the 
 average vector distances (`m`: `B0`; joint total `m+1` subspace: `Bt`; 
 all possible `m+1` subspaces: `B1`), from the M multivariate sequences
 in `Data` using the default parameters: 
 embedding dimension = 2*ones(M,1), time delay = ones(M,1), 
 fuzzy membership function = "default", fuzzy function parameters= [0.2, 2],
 logarithm = natural, data normalization = false,

 
!!! note 

    The entropy value returned as `MFuzz` is estimated using the "full" 
    method [i.e.  -log(Bt/B0)] which compares delay vectors across all possible `m+1` 
    expansions of the embedding space as applied in [1][3]. Contrary to
    conventional definitions of sample entropy, this method does not provide a
    lower bound of 0!!
    Thus, it is possible to obtain negative entropy values for multivariate 
    fuzzy entropy, even for stochastic processes...

    Alternatively, one can calculate `MFuzz` via the "naive" method, 
    which ensures a lower bound of 0, by using the average vector distances
    for an individual `m+1` subspace (B1) [e.g. -log(B1(1)/B0)],
    or the average for all `m+1` subspaces [i.e. -log(mean(B1)/B0)].

    To maximize the number of points in the embedding process, this algorithm 
    uses N - max(m * tau) delay vectors and _*not*_ N - max(m) * max(tau) as employed 
    in [1] and [3].

 
-------------------------------------------------------------
      
    MFuzz, B0, Bt, B1 = MvFuzzEn(Data::AbstractArray{T} where T<:Real; m::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, tau::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, r::Union{Real,Tuple{Real,Real}}=(.2,2.0), Fx::String="default", Logx::Real=exp(1), Norm::Bool=false)
        
 Returns the multivariate sample entropy estimates (`MSamp`) estimated
 from the M multivariate data sequences in `Data` using the specified 
 keyword arguments:

 # Arguments:
 `Data`  - Multivariate dataset, NxM matrix of N (>10) observations (rows) and M (cols) univariate data sequences\n
 `m`     - Embedding Dimension, a vector of M positive integers\n
 `tau`   - Time Delay, a vector of M positive integers\n
 `Fx`    - Fuzzy function name, one of the following: 
            {`"sigmoid", "modsampen", "default", "gudermannian",`
            `"bell", "triangular", "trapezoidal1", "trapezoidal2",`
            `"z_shaped", "gaussian", "constgaussian"`}\n 
 `r`      - Fuzzy function parameters, a 1 element scalar or a 2 element
                tuple of positive values. The `r` parameters for each fuzzy
                function are defined as follows:      [default: [.2 2]]\n
                default:        r(1) = divisor of the exponential argument
                                r(2) = argument exponent (pre-division)
                sigmoid:        r(1) = divisor of the exponential argument
                                r(2) = value subtracted from argument (pre-division)
                modsampen:      r(1) = divisor of the exponential argument
                                r(2) = value subtracted from argument (pre-division)
                gudermannian:   r  = a scalar whose value is the numerator of
                                    argument to gudermannian function:
                                    GD(x) = atan(tanh(`r`/x))
                triangular:     r = a scalar whose value is the threshold (corner point) of the triangular function.
                trapezoidal1:   r = a scalar whose value corresponds to the upper (2r) and lower (r) corner points of the trapezoid.
                trapezoidal2:   r(1) = a value corresponding to the upper corner point of the trapezoid.
                                r(2) = a value corresponding to the lower corner point of the trapezoid.
                z_shaped:       r = a scalar whose value corresponds to the upper (2r) and lower (r) corner points of the z-shape.
                bell:           r(1) = divisor of the distance value
                                r(2) = exponent of generalized bell-shaped function
                gaussian:       r = a scalar whose value scales the slope of the Gaussian curve.
                constgaussian:  r = a scalar whose value defines the lower threshod and shape of the Gaussian curve.\n
 `Logx`   - Logarithm base, a positive scalar \n 
 `Norm`   - Normalisation of all M sequences to unit variance, a boolean\n
        
 # See also `MvSampEn`, `FuzzEn`, `XFuzzEn`, `FuzzEn2D`, `MSEn`, `MvPermEn`

 # References:
    [1] Ahmed, Mosabber U., et al. 
        "A multivariate multiscale fuzzy entropy algorithm with application
        to uterine EMG complexity analysis." 
        Entropy 19.1 (2016): 2.

    [2] Azami, Alberto Fernández, Javier Escudero. 
        "Refined multiscale fuzzy entropy based on standard deviation for 
        biomedical signal analysis." 
        Medical & biological engineering & computing 55 (2017): 2037-2052.

    [3] Ahmed Mosabber Uddin, Danilo P. Mandic
        "Multivariate multiscale entropy analysis."
        IEEE signal processing letters 19.2 (2011): 91-94.

"""
    function MvFuzzEn(Data::AbstractArray{T,2} where T<:Real; m::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, tau::Union{AbstractArray{T} where T<:Int, Nothing}=nothing, 
                r::Union{Real,Tuple{Real,Real}}=(.2,2.0), Fx::String="default", Logx::Real=exp(1), Norm::Bool=false)

    N, Dn = size(Data)
    isnothing(m) ? m = 2*ones(Int, Dn) : nothing
    isnothing(tau) ? tau = ones(Int, Dn) : nothing
    Norm ? Data = Data./std(Data, dims=1, corrected=false) : nothing

    (N>10) && (Dn>1) ? nothing : error("Data:   must be an NxM matrix where N>10 and M>1")
    ndims(m)==1  && length(m)==Dn  && all(m.>0) && eltype(m)<:Int ? nothing : error("m:  vector of M positive integers")
    ndims(tau)==1  && length(tau)==Dn  && all(tau.>0) && eltype(tau)<:Int ? nothing : error("tau:  vector of M positive integers")
    (minimum(r)>=0 && length(r)<=2) ? 
    nothing : error("r:  must be a scalar or 2 element tuple of positive values")
    (lowercase(Fx) in ["default","sigmoid","modsampen","gudermannian","bell", "z_shaped",
    "triangular", "trapezoidal1","trapezoidal2","gaussian","constgaussian"]) ?
    nothing : error("Fx:    must be one of the following strings -
    'default', 'sigmoid', 'modsampen', 'gudermannian', 'bell', 'z_shaped',
    'triangular', 'trapezoidal1','trapezoidal2','gaussian','constgaussian'")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")

    Fun = getfield(_MvFuzzEn,Symbol(lowercase(Fx)))

    Nx = N - maximum((m.-1).*tau)
    Ny = N - maximum(m.*tau)
    Vex = zeros(Nx,sum(m))
    q = 1;
    for k = 1:Dn
        for p = 1:m[k]
            Vex[:,q] = Data[1+(p-1)*tau[k]:Nx+(p-1)*tau[k],  k]
            q += 1
        end    
    end
    Count0 = Distx(Vex .- mean(Vex,dims=2), r, Fun)
    B0 = sum(Count0)/(Nx*(Nx-1)/2)

    B1 = zeros(Dn)
    Temp = cumsum(m)
    Vez = Inf.*ones(1,sum(m)+1)
    for k = 1:Dn
        Sig = Data[1+m[k]*tau[k]:Ny+m[k]*tau[k], k]
        Vey = hcat(Vex[1:Ny, 1:Temp[k]], Sig,  Vex[1:Ny, Temp[k]+1:end])
        Vez = vcat(Vez, Vey)
        Count1 = Distx(Vey .- mean(Vey,dims=2), r, Fun);
        B1[k] = sum(Count1)/(Ny*(Ny-1)/2)
    end
    Vez = Vez[2:end,:]
    Count1 = Distx(Vez .- mean(Vez,dims=2), r, Fun)
    Bt = sum(Count1)/(Dn*Ny*((Dn*Ny)-1)/2)
    
    MFuzz = -log(Bt/B0)/log(Logx)
    return MFuzz, B0, Bt, B1
    end


    function Distx(Vex, r, Fun)
        Nt = size(Vex)[1]
        Counter = zeros(Nt-1,Nt-1)
        for x=1:Nt-1
            Counter[x,x:end] = Fun(maximum(abs.(Vex[x+1:end,:] .- Vex[x:x,:]), dims=2)[:], r)
        end
        return Counter
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

    function triangular(x,r)
        length(r)==1 ? nothing : error("When Fx = 'Triangular', r must be a scalar > 0.")
        y = 1 .- (x./r)
        y[x .> r] .= 0
        return y
    end

    function trapezoidal1(x, r)
        length(r)==1 ? nothing : error("When Fx = 'Trapezoidal1', r must be a scalar > 0.")
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
        y[x .<= 2*r] = 2*(((x[x .<= 2*r] .- 2*r)./r).^2)
        y[x .<= 1.5*r] = 1 .- (2*(((x[x .<= 1.5*r] .- r)/r).^2))
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

"""Copyright 2024 Matthew W. Flood, EntropyHub
  
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

For Terms of Use see https://github.com/MattWillFlood/EntropyHub"""