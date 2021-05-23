# Example 3: Phase Entropy w/ Second Order Difference Plot

Import the x and y components of the Henon system of equations.

```
Data = ExampleData("henon");

using Plots
scatter(Data[:,1], Data[:,2],
markercolor = "green", markerstrokecolor = "black",
markersize = 3, background_color = "black",grid = false)
```

![Henon](../assets/henonjl.png)

Calculate the phase entropy of the y-component in bits (logarithm base 2) without normalization using 7 angular partitions and return the second-order difference plot.
```@example
using EntropyHub #  hide
Data = ExampleData("henon"); # hide
Y = Data[:,2];
Phas = PhasEn(Y, K = 7, Norm = false, Logx = 2, Plotx = true)
```

![Phas1](../assets/phasx1jl.png)

Calculate the phase entropy of the x-component using 11 angular partitions, a time delay of 2, and return the second-order difference plot.
```@example
using EntropyHub # hide
Data = ExampleData("henon"); # hide
X = Data[:,1];
Phas = PhasEn(X, K = 11, tau = 2, Plotx = true)
```

![Phas2](../assets/phasx2jl.png)
