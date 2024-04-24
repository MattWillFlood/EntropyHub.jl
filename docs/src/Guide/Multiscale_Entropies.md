```@meta
CollapsedDocStrings = true
Description = "Multiscale entropy of single time series signals"
```
# Multiscale Entropies

__*Functions for estimating the multiscale entropy  of a univariate time series.*__

Multiscale entropy can be calculated using any of the Base entropies: 
    ([`ApEn`](@ref), [`AttnEn`](@ref), [`BubbEn`](@ref), [`CondEn`](@ref), [`CoSiEn`](@ref), [`DistEn`](@ref), [`DivEn`](@ref), 
    [`DispEn`](@ref), [`EnofEn`](@ref), [`FuzzEn`](@ref), [`GridEn`](@ref), [`IncrEn`](@ref), [`K2En`](@ref),
    [`PermEn`](@ref), [`PhasEn`](@ref), [`RangEn`](@ref), [`SampEn`](@ref), [`SlopEn`](@ref), [`SpecEn`](@ref), [`SyDyEn`](@ref)).

!!! info "NOTE:"

    Multiscale cross-entropy functions have two positional arguments:

    1. the data sequence, `Sig` (a vector > 10 elements),
    2. the multiscale entropy object, `Mobj` -> see [`MSobject`](@ref)



```@docs
EntropyHub.MSobject
EntropyHub.MSEn
EntropyHub.cMSEn
EntropyHub.rMSEn
EntropyHub.hMSEn
```

