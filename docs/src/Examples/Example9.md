# Example 9: Hierarchical Multiscale corrected Cross-Conditional Entropy

Import the x and y components of the Henon system of equations and create a multiscale
entropy object with the following parameters:
`EnType` = XCondEn(), embedding dimension = 2, time delay = 2, number of symbols = 12,
logarithm base = 2, normalization = true
```
Data = ExampleData("henon");
Mobj = MSobject(XCondEn, m = 2, tau = 2, c = 12,  Logx = 2, Norm = true)

using Plots
scatter(Data[:,1], Data[:,2], markercolor = "green", markerstrokecolor = "black",
markersize = 3, background_color = "black", grid = false)
```

![Henon](../assets/henonjl.png)

Calculate the hierarchical multiscale corrected cross-conditional entropy over 4 temporal
scales and return the average cross-entropy at each scale (`Sn`), the complexity index (`Ci`),
and a plot of the multiscale entropy curve and the hierarchical tree with the cross-entropy
value at each node.
```@example
using EntropyHub # hide
Data = ExampleData("henon"); # hide
Mobj = MSobject(XCondEn, m = 2, tau = 2, c = 12,  Logx = 2, Norm = true) # hide
MSx, Sn, Ci = hXMSEn(Data[:,1], Data[:,2], Mobj, Scales = 4, Plotx = true)
```

![hXMSEn](../assets/hXMSEnjl.png)