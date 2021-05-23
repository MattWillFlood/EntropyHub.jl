# Cross Entropies

__*Functions for estimating the cross-entropy between two univariate time series.*__

`The following functions also form the cross-entropy method used by Multiscale Cross-Entropy functions.`


!!! tip "Signals for mutliscale cross-entropy functions"

    For cross-entropy and multiscale cross-entropy functions, the two time series signals are passed as a two-column or two-row matrix. 
    At present, it is not possible  to pass signals of different lengths separately. 


These functions are directly available when EntropyHub is imported:

```
julia> using EntropyHub
julia> names(EntropyHub)
```
```
 :ApEn
 :AttnEn
 :BubbEn
   ⋮
 :hXMSEn
 :rMSEn
 :rXMSEn
```

```@docs
EntropyHub.XApEn
EntropyHub.XSampEn
EntropyHub.XFuzzEn
EntropyHub.XK2En
EntropyHub.XPermEn
EntropyHub.XCondEn
EntropyHub.XDistEn
EntropyHub.XSpecEn
```