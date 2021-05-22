module _SpecEn
export SpecEn
using FFTW: fft
using DSP: conv
    """
        Spec, BandEn = SpecEn(Sig) 

    Returns the spectral entropy estimate of the full spectrum (`Spec`)
    and the within-band entropy (`BandEn`) estimated from the data sequence (`Sig`)
    using the default  parameters: 
    N-point FFT = 2*len(`Sig`) + 1, normalised band edge frequencies = [0 1],
    logarithm = base 2, normalisation = w.r.t # of spectrum/band frequency values.

        Spec, BandEn = SpecEn(Sig::AbstractArray{T,1} where T<:Real; N::Int=1 + (2*size(Sig,1)), Freqs::Tuple{Real,Real}=(0,1), Logx::Real=exp(1), Norm::Bool=true)

    Returns the spectral entropy (`Spec`) and the within-band entropy (`BandEn`)
    estimate for the data sequence (`Sig`) using the specified 'keyword' arguments:

    # Arguments:
    `N`     - Resolution of spectrum (N-point FFT), an integer > 1  \n
    `Freqs` - Normalised band edge frequencies, a 2 element tuple with values \n
              in range [0 1] where 1 corresponds to the Nyquist frequency (Fs/2).
              Note: When no band frequencies are entered, BandEn == SpecEn
    `Logx`  - Logarithm base, a positive scalar (enter 0 for natural log) \n
    `Norm`  - Normalisation of `Spec` value:\n
              [false]  no normalisation.
              [true]   normalises w.r.t # of spectrum/band frequency values - default.

    For more info, see the EntropyHub guide.

    # See also `XSpecEn`, `fft`, `MSEn`,  `XMSEn`
  
    # References:
        [1] G.E. Powell and I.C. Percival,
            "A spectral entropy method for distinguishing regular and 
            irregular motion of Hamiltonian systems." 
            Journal of Physics A: Mathematical and General 
            12.11 (1979): 2053.
  
        [2] Tsuyoshi Inouye, et al.,
            "Quantification of EEG irregularity by use of the entropy of 
            the power spectrum." 
            Electroencephalography and clinical neurophysiology 
            79.3 (1991): 204-210.
  

    """
    function SpecEn(Sig::AbstractArray{T,1} where T<:Real; N::Int=1 + (2*size(Sig,1)), 
        Freqs::Tuple{Real,Real}=(0,1), Logx::Real=exp(1), Norm::Bool=true)
        
    (size(Sig)[1] > 4) ? nothing : error("Sig:   must be a numeric vector")
    (N > 1) ? nothing :  error("N:     must be an integer > 1")
    (0<=Freqs[1]<1 && 0<Freqs[2]<=1 && Freqs[1]<Freqs[2]) ? nothing :
        error("Freq:    must be a two element tuple with values in range [0 1].
                The values must be in increasing order.")
    (Logx>0) ? nothing :   error("Logx:     must be a positive number > 0")
    
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

    Temp = conv(Sig,Sig)
    N <= size(Temp,1) ? Temp = Temp[1:N]  : Temp = vcat(Temp,zeros(N-size(Temp)[1]))
    
    Pt = abs.(fft(Temp))
    Pxx = Pt[1:Fx]/sum(Pt[1:Fx])
    Spec = -(transpose(Pxx)*log.(Logx, Pxx))
    Pband = (Pxx[Freqs[1]:Freqs[2]])/sum(Pxx[Freqs[1]:Freqs[2]])
    BandEn = -(transpose(Pband)*log.(Logx, Pband))

    if Norm
        Spec = Spec/(log(Logx, Fx));
        BandEn = BandEn/(log(Logx, Freqs[2]-Freqs[1]+1));
    end

    return Spec, BandEn
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