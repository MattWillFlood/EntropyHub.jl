```@meta
CollapsedDocStrings = true
Description = "Multiscale entropy between two time series signals"
```
# Multiscale Cross-Entropies

__*Functions for estimating the multiscale entropy between two univariate time series.*__

Just as one can calculate multiscale entropy using any Base entropy, the same functionality is possible with multiscale cross-entropy using any Cross-entropy function:
    [`XApEn`](@ref), [`XSampEn`](@ref), [`XK2En`](@ref), [`XCondEn`](@ref), [`XPermEn`](@ref), [`XSpecEn`](@ref), [`XDistEn`](@ref), [`XFuzzEn`](@ref).

To do so, we again use the [`MSobject`](@ref) function to pass a multiscale object (`Mobj`) to the multiscale cross-entropy functions.

!!! info "NOTE:"

    Multiscale cross-entropy functions have three positional arguments:

    1. the first data seuqence, `Sig1` (a vector of >10 elements),
    2. the second data seuqence, `Sig2` (a vector of > 10 elements),
    3. the multiscale entropy object, `Mobj` -> see [`MSobject`](@ref)


```@docs
EntropyHub.XMSEn
EntropyHub.cXMSEn
EntropyHub.rXMSEn
EntropyHub.hXMSEn
```