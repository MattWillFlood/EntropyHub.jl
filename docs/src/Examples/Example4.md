# Example 4: Cross-Distribution Entropy w/ Different Binning Methods

Import a signal of pseudorandom integers in the range [1, 8] and calculate the cross-
distribution entropy with an embedding dimension of 5, a time delay (`tau`) of 3, and 'Sturges' bin selection method.

```@example
using EntropyHub # hide
X = ExampleData("randintegers2");
XDist, _ = XDistEn(X[:,1], X[:,2], m = 5, tau = 3);
println(" ") # hide
println("XDist = ", XDist) # hide
``` 

Use Rice's method to determine the number of histogram bins and return the probability of each bin (`Ppi`).

```@example 
using EntropyHub # hide
X = ExampleData("randintegers2"); # hide
XDist, Ppi = XDistEn(X[:,1], X[:,2], m = 5, tau = 3, Bins = "rice")
println("XDist = ", XDist) # hide
println("Ppi = ", Ppi) #hide
```
