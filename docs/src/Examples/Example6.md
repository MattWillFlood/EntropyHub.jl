# Example 6: Multiscale Increment Entropy

Import a signal of uniformly distributed pseudorandom integers in the range [1 8] and
create a multiscale entropy object with the following parameters:
`EnType` = IncrEn(), embedding dimension = 3, a quantifying resolution = 6, normalization = true.

```@example
using EntropyHub # hide
X = ExampleData("randintegers");
Mobj = MSobject(IncrEn, m = 3, R = 6, Norm = true)
```

Calculate the multiscale increment entropy over 5 temporal scales using the modified
graining procedure where:

``y_j^{(\tau)} =\frac{1}{\tau } \sum_{i=\left(j-1\right)\tau +1}^{j\tau } x{_i},    1<= j <= \frac{N}{\tau }``



```@example
using EntropyHub # hide
X = ExampleData("randintegers"); # hide
Mobj = MSobject(IncrEn, m = 3, R = 6, Norm = true) # hide
MSx, _ = MSEn(X, Mobj, Scales = 5, Methodx = "modified");
MSx # hide
```