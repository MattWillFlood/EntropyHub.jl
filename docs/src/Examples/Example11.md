# Example 11: Multivariate Dispersion Entropy

Import a vector of 4096 uniformly distributed random integers in range [1 8] and convert it to a multivariate set of 4 sequences with 1024 samples each. 

```@example
using EntropyHub # hide
X = ExampleData("randintegers");
Data = reshape(X, 1024, 4);
Data # hide
```

Calculate the multivariate dispersion entropy and reverse dispersion entropy for embedding dimensions (m) = [1,1,2,3], using a 7-symbol transform.

```@example
using EntropyHub # hide
X = ExampleData("randintegers"); # hide
Data = reshape(X, 1024, 4); # hide
MDisp, RDE = MvDispEn(Data, m = [1,1,2,3], c = 7);
MDisp  # hide
```   
```@example
using EntropyHub # hide
X = ExampleData("randintegers"); # hide
Data = reshape(X, 1024, 4); # hide
MDisp, RDE = MvDispEn(Data, m = [1,1,2,3], c = 7); # hide
RDE  # hide
``` 

Perform the same calculation but normalize the output entropy estimate w.r.t the number of unique dispersion patterns

```@example
using EntropyHub # hide
X = ExampleData("randintegers"); # hide
Data = reshape(X, 1024, 4); # hide
MDisp, RDE = MvDispEn(Data, m = [1,1,2,3], c = 7, Norm = true);
MDisp  # hide
``` 
```@example
using EntropyHub # hide
X = ExampleData("randintegers"); # hide
Data = reshape(X, 1024, 4); # hide
MDisp, RDE = MvDispEn(Data, m = [1,1,2,3], c = 7, Norm = true); # hide
RDE  # hide
``` 
        
Compare the results above (``Methodx == 'v1'``) with those obtained using the *mvDE* method (``Methodx=='v2'``), returning estimates for each value from 1, ..., max(m)

```@example
using EntropyHub # hide
X = ExampleData("randintegers"); # hide
Data = reshape(X, 1024, 4); # hide
MDisp, RDE = MvDispEn(Data, m = [1,1,2,3], c = 7, Norm = true, Methodx = "v2")
MDisp  # hide
``` 
```@example
using EntropyHub # hide
X = ExampleData("randintegers"); # hide
Data = reshape(X, 1024, 4); # hide
MDisp, RDE = MvDispEn(Data, m = [1,1,2,3], c = 7, Norm = true, Methodx = "v2") # hide
RDE  # hide
``` 

