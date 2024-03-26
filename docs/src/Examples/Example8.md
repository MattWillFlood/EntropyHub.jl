# Example 8: Composite Multiscale Cross-Approximate Entropy

Import two signals of uniformly distributed pseudorandom integers in the range [1 8] and
create a multiscale entropy object with the following parameters:
`EnType` = XApEn(), embedding dimension = 2, time delay = 2, radius distance threshold = 0.5
```@example
using EntropyHub # hide
X = ExampleData("randintegers2");
Mobj = MSobject(XApEn, m = 2, tau = 2, r = 0.5)
Mobj # hide
``` 

Calculate the comsposite multiscale cross-approximate entropy over 3 temporal scales
where the radius distance threshold value (`r`) specified by `Mobj` becomes scaled by the
variance of the signal at each scale.
```@example
using EntropyHub # hide
X = ExampleData("randintegers2"); # hide 
Mobj = MSobject(XApEn, m = 2, tau = 2, r = 0.5) # hide
MSx, _ = cXMSEn(X[:,1], X[:,2], Mobj, Scales = 3, RadNew = 1)
MSx # hide
```