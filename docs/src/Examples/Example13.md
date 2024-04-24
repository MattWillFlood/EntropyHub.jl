# Example 13: Window Data Tool

This example shows how to use the [`WindowData()`](@ref) function to divide univariate/multivariate data into
subsequences.

_______________________________________________________________________________

### Ex. 1
Import a sequence of uniformly distributed random numbers.
```
X = ExampleData("uniform");
```

Extract a series of subsequences by dividing the sequence into windows of length 1024 with no overlap.
```@example
using EntropyHub # hide
X = ExampleData("uniform"); # hide
WinData, Log = WindowData(X, WinLen = 1024);
WinData # hide
```
```@example
using EntropyHub # hide
X = ExampleData("uniform"); # hide
WinData, Log = WindowData(X, WinLen = 1024); # hide
Log # hide
```

Repeat the previous step but change the number of overlapping window samples to 256.
```@example
using EntropyHub # hide
X = ExampleData("uniform"); # hide
WinData, Log = WindowData(X, WinLen = 1024, Overlap = 256);
WinData # hide
```
```@example
using EntropyHub # hide
X = ExampleData("uniform"); # hide
WinData, Log = WindowData(X, WinLen = 1024, Overlap = 256); # hide
Log # hide
```
_______________________________________________________________________________

### Ex. 2
Window a range of numbers (1:1234) into 5 windows.
```@example
using EntropyHub # hide
WinData, Log = WindowData(1:1234);
WinData # hide
```
```@example
using EntropyHub # hide
WinData, Log = WindowData(1:1234); # hide
Log # hide
```

Repeat the previous step, but this time retain any remaining samples that do not fill the last window.
```@example
using EntropyHub # hide
WinData, Log = WindowData(1:1234, Mode="include");
WinData # hide
```
```@example
using EntropyHub # hide
_, Log = WindowData(1:1234); # hide
Log # hide
```
Note that the last vector in `WinData` contains only 4 values.
```@example
using EntropyHub # hide
WinData, _ = WindowData(1:1234, Mode="include"); # hide
WinData[end] # hide
```
_______________________________________________________________________________

### Ex. 3
Generate a multivariate dataset of 5 uniformly-distributed random number sequences (N=3333) and 
divide the dataset into subsets of 777 samples.
```@example
using EntropyHub # hide
using Random
X = rand(MersenneTwister(0), 3333, 5)
WinData, Log = WindowData(X, WinLen = 777)
WinData # hide
```
```@example
using EntropyHub # hide
using Random # hide
X = rand(MersenneTwister(0), 3333, 5) # hide
WinData, Log = WindowData(X, WinLen = 777) # hide
Log # hide
```

Repeat the previous step including 55 samples of overlap and retain any remaining samples that do not fill a window.
```@example
using EntropyHub # hide
using Random # hide
X = rand(MersenneTwister(0), 3333, 5) # hide
WinData, Log = WindowData(X, WinLen = 777, Overlap = 55, Mode = "include")
WinData # hide
```
```@example
using EntropyHub # hide
using Random # hide
X = rand(MersenneTwister(0), 3333, 5) # hide
WinData, Log = WindowData(X, WinLen = 777, Overlap = 55, Mode = "include") # hide
Log # hide
```
