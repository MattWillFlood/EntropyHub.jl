# Multiscale Entropies

__*Functions for estimating the multiscale entropy between of a univariate time series.*__

Multiscale entropy can be calculated using any of the Base entropies: 
    (`ApEn`, `AttnEn`, `BubbEn`, `CondEn`, `CoSiEn`, `DistEn`, 
    `DispEn`, `EnofEn`, `FuzzEn`, `GridEn`, `IncrEn`, `K2En`,
    `PermEn`, `PhasEn`, `SampEn`, `SlopEn`, `SpecEn`, `SyDyEn`).

!!! info "NOTE:"

    Multiscale cross-entropy functions have two positional arguments:

    1. the time series signal, `Sig` (a vector > 10 elements),
    2. the multiscale entropy object, `Mobj`.


```@docs
EntropyHub.MSobject
```

The following functions use the multiscale entropy object shown above.

```@docs
EntropyHub.MSEn
EntropyHub.cMSEn
EntropyHub.rMSEn
EntropyHub.hMSEn
```

