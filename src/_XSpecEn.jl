module _XSpecEn
export XSpecEn
using FFTW: fft
using DSP: conv
    """
        XSpec, BandEn = XSpecEn(Sig) 

    Returns the cross-spectral entropy estimate (`XSpec`) of the full cross-
    spectrum and the within-band entropy (`BandEn`) estimated between the data 
    sequences contained in `Sig` using the default  parameters: 
    N-point FFT = 2 * max(length(`Sig1`/`Sig2`)) + 1, normalised band edge frequencies = [0 1],
    logarithm = base 2, normalisation = w.r.t # of spectrum/band frequency 
    values.

        XSpec, BandEn = XSpecEn(Sig1::Union{AbstractMatrix{T}, AbstractVector{T}} where T<:Real, Sig2::Union{AbstractVector{T} where T<:Real, Nothing} = nothing;  N::Union{Nothing,Int}=nothing, Freqs::Tuple{Real,Real}=(0,1), Logx::Real=exp(1), Norm::Bool=true)

    Returns the cross-spectral entropy (`XSpec`) and the within-band entropy 
    (`BandEn`) estimate between the data sequences contained in `Sig1` and `Sig2` using the
    following specified 'keyword' arguments:

    # Arguments:
    `N`     - Resolution of spectrum (N-point FFT), an integer > 1    \n
    `Freqs` - Normalised band edge frequencies, a scalar in range [0 1]
              where 1 corresponds to the Nyquist frequency (Fs/2).
              **Note: When no band frequencies are entered, BandEn == SpecEn** \n
    `Logx`  - Logarithm base, a positive scalar     [default: base 2]
              ** enter 0 for natural log**  \n
    `Norm`  - Normalisation of `XSpec` value:
              [false]  no normalisation.
              [true]   normalises w.r.t # of spectrum/band frequency values [default] \n

    For more info, see the EntropyHub guide

    # See also `SpecEn`, `fft`, `XDistEn`, `periodogram`, `XSampEn`, `XApEn`
    
    # References:
        [1]  Matthew W. Flood,
            "XSpecEn - EntropyHub Project"
            (2021) https://github.com/MattWillFlood/EntropyHub

 
    """
    function XSpecEn(Sig1::Union{AbstractMatrix{T}, AbstractVector{T}} where T<:Real, Sig2::Union{AbstractVector{T} where T<:Real, Nothing} = nothing; 
        N::Union{Nothing,Int}=nothing, Freqs::Tuple{Real,Real}=(0,1), Logx::Real=exp(1), Norm::Bool=true)
        
    if all(isa.((Sig1,Sig2), AbstractVector))
        N1 = size(Sig1,1);  N2 = size(Sig2,1)  
        S1 = copy(Sig1); S2 = copy(Sig2)
    elseif (minimum(size(Sig1))==2 && (Sig2 isa Nothing)) 
        argmin(size(Sig1)) == 2 ? nothing : Sig1 = Sig1'
        S1 = Sig1[:,1]; S2 = Sig1[:,2];
        N1 = maximum(size(Sig1)); N2 = maximum(size(Sig1));
    else   error("""Sig1 and Sig2 must be 2 separate vectors 
                \t\t\t - OR - 
                Sig1 must be 2-column matrix and Sig2 nothing""")
    end
    
    N isa Nothing ? N = 2*max(N1,N2) + 1 : N = 2*N1 + 1;

    (N1>=10 && N2>=10) ? nothing :  error("Sig1/Sig2:   sequences must have >= 10 values")
    (N > 1) ? nothing :  error("N:     must be an integer > 1")
    (0<=Freqs[1]<1 && 0<Freqs[2]<=1 && Freqs[1]<Freqs[2]) ? nothing :
        error("Freq:    must be a two element tuple with values in range [0 1].
                The values must be in increasing order.")
    (Logx>0) ? nothing : error("Logx:     must be a positive number > 0")
        
    Freqs = collect(Freqs)
    Fx = Int(ceil(N/2))
    Freqs = Int.(round.(Freqs.*Fx))
    Freqs[Freqs.==0] .= 1

    if Freqs[1] > Freqs[2]
        error("Lower band frequency must come first.")
    elseif Freqs[2]-Freqs[1]<1
        error("Spectrum resoution too low to determine bandwidth.") 
    elseif minimum(Freqs)<0 || maximum(Freqs)>Fx
        error("Freqs must be normalized w.r.t sampling frequency [0 1].")
    end

    Temp = conv(S1,S2)
    N <= size(Temp,1) ? Temp = Temp[1:N] : Temp = vcat(Temp,zeros(N-size(Temp,1)))
    
    Pt = abs.(fft(Temp));
    Pxx = Pt[1:Fx]/sum(Pt[1:Fx]);
    XSpec = -transpose(Pxx)*log.(Logx, Pxx)
    Pband = (Pxx[Freqs[1]:Freqs[2]])/sum(Pxx[Freqs[1]:Freqs[2]]);
    BandEn = -transpose(Pband)*log.(Logx, Pband)

    if Norm
        XSpec = XSpec/log(Logx, Fx)
        BandEn = BandEn/log(Logx, Freqs[2]-Freqs[1]+1)
    end

    return XSpec, BandEn
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