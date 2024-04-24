```@meta
CollapsedDocStrings = true
Description = "Multivariate multiscale entropy of multivariate time series "
```
# Multivariate Multiscale Entropies

__*Functions for estimating the Multivariate multiscale entropy of a multivariate time series dataset.*__

Just as one can calculate multiscale entropy using any Base or Cross- entropy, the same functionality is possible with multivariate multiscale entropy using any Multivariate function:
    [`MvSampEn`](@ref), [`MvFuzzEn`](@ref), [`MvDispEn`](@ref), [`MvPermEn`](@ref), [`MvCoSiEn`](@ref).

To do so, we again use the [`MSobject`](@ref) function to pass a multiscale object (`Mobj`) to the multivariate multiscale entropy functions.

!!! info "NOTE:"

    Multivariate multiscale entropy functions have two positional arguments:

    1. the multivariate dataset, `Data` (an NxM matrix of N observations (>10 elements), and M time series (>1)),
    2. the multiscale entropy object, `Mobj` -> see [`MSobject`](@ref)

```@docs
EntropyHub.MvMSEn
EntropyHub.cMvMSEn
```