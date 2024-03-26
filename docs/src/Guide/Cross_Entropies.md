# Cross Entropies

__*Functions for estimating the cross-entropy between two univariate time series.*__

`The following functions also form the cross-entropy method used by Multiscale Cross-Entropy functions.`

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