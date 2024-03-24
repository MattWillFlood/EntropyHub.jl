# Base Entropies

__*Functions for estimating the entropy of a single univariate time series.*__

`The following functions also form the base entropy method used by Multiscale functions.`

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
ApEn
SampEn
FuzzEn
K2En
PermEn
CondEn
DistEn
SpecEn
DispEn
SyDyEn
IncrEn
CoSiEn
PhasEn
SlopEn
BubbEn
GridEn
EnofEn
AttnEn
RangEn
DivEn
```