# Example 1: Sample Entropy

Import a signal of normally distributed random numbers [mean = 0; SD = 1], and calculate the
sample entropy for each embedding dimension (m) from 0 to 4.

```@example
using EntropyHub # hide
X = ExampleData("gaussian");
Samp, _ = SampEn(X, m = 4);
Samp # hide
```

Select the last value to get the sample entropy for m = 4.
```@example
using EntropyHub # hide
X = ExampleData("gaussian"); # hide
Samp, _ = SampEn(X, m = 4); # hide
Samp[end]
```

Calculate the sample entropy for each embedding dimension (m) from 0 to 4 with a time delay (tau) of 2 samples.
```@example
using EntropyHub # hide
X = ExampleData("gaussian"); #hide
Samp, Phi1, Phi2 = SampEn(X, m = 4, tau = 2)
```
