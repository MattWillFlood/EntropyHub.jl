```@meta
CollapsedDocStrings = true
Description = "Entropy of two-dimensional (2D) data"
```
# Bidimensional Entropies

__*Functions for estimating the entropy of a two-dimensional univariate matrix.*__

While EntropyHub functions primarily apply to time series data, with the following
bidimensional entropy functions one can estimate the entropy of two-dimensional (2D)
matrices. Hence, bidimensional entropy functions are useful for applications such as image analysis.

!!! danger "IMPORTANT: Locked Matrix Size"

    Each bidimensional entropy function ([`SampEn2D`](@ref), [`FuzzEn2D`](@ref), [`DistEn2D`](@ref), [`DispEn2D`](@ref), 
    [`EspEn2D`](@ref), [`PermEn2D`](@ref)) has an important keyword argument - `Lock`. Bidimensional entropy functions are
    "locked" by default (`Lock == true`) to only permit matrices with a maximum size of 128 x 128.

    The reason for this is because there are hundreds of millions of pairwise calculations
    performed in the estimation of bidimensional entropy, so memory errors often
    occur when storing data on RAM.

    e.g. For a matrix of size [200 x 200], an embedding dimension (`m`) = 3, and a time
    delay (`tau`) = 1, there are 753,049,836 pairwise matrix comparisons (6,777,448,524
    elemental subtractions).
    To pass matrices with sizes greater than [128 x 128], set `Lock = false`.

    `CAUTION: unlocking the permitted matrix size may cause your Julia IDE to crash.`

    

```@docs
EntropyHub.SampEn2D
EntropyHub.FuzzEn2D
EntropyHub.DistEn2D
EntropyHub.DispEn2D
EntropyHub.PermEn2D
EntropyHub.EspEn2D
```