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
X = ExampleData("gaussian"); # hide
Samp, Phi1, Phi2 = SampEn(X, m = 4, tau = 2)
```

Import a signal of uniformly distributed random numbers in the range [-1, 1] and calculate the sample entropy for an embedding dimension (m) of 5, a time delay of 2, and a threshold radius of 0.075.      
Return the conditional probability (Vcp) and the number of overlapping matching vector pairs of lengths m+1 (Ka) and m (Kb), respectively.
```@example
using EntropyHub # hide
X = ExampleData("uniform"); # hide
Samp, _, _, Vcp_Ka_Kb = SampEn(X, m = 5, tau = 2, r = 0.075, Vcp = true)
Vcp, Ka, Kb = Vcp_Ka_Kb
println("Vcp = ", Vcp) # hide
println("Ka = ", Ka) # hide
println("Kb = ", Kb) # hide
```



