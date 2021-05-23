# Multiscale Cross-Entropies

__*Functions for estimating the multiscale entropy between two univariate time series.*__

Just as one can calculate multiscale entropy using any Base entropy, the same functionality is possible with multiscale cross-entropy using any Cross-entropy function:
    (`XApEn`, `XSampEn`, `XK2En`, `XCondEn`, `XPermEn`, `XSpecEn`, `XDistEn`, `XFuzzEn`).

To do so, we again use the `MSobject` function to pass a multiscale object (`Mobj`) to the multiscale cross-entropy functions.

!!! info "NOTE:"

    Multiscale cross-entropy functions have two positional arguments:

    1. the time series signals, `Sig` (an Nx2 matrix),
    2. the multiscale entropy object, `Mobj`.

!!! tip "Signals for mutliscale cross-entropy functions"

    For cross-entropy and multiscale cross-entropy functions, the two time series signals are passed as a two-column or two-row matrix. 
    At present, it is not possible  to pass signals of different lengths separately. 


[`EntropyHub.MSobject`](@ref)

The following functions use the multiscale entropy object shown above.


```@docs
EntropyHub.XMSEn
EntropyHub.cXMSEn
EntropyHub.rXMSEn
EntropyHub.hXMSEn
```