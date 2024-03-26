# Example 10: Bidimensional Fuzzy Entropy

Import an image of a Mandelbrot fractal as a matrix.
```
X = ExampleData("mandelbrot_Mat");

using Plots
heatmap(X, background_color="black")
```
![AlmondBread](../assets/mandelbrotjl.png)

Calculate the bidimensional fuzzy entropy in trits (logarithm base 3) with a template
matrix of size [8 x 5], and a time delay (`tau`) of 2 using a `'constgaussian'` fuzzy membership function (r=24).

```@example
using EntropyHub # hide
X = ExampleData("mandelbrot_Mat"); # hide
FE2D = FuzzEn2D(X, m = (8, 5), tau = 2, Fx = "constgaussian", r = 24, Logx = 3)
```