```@contents
Pages = ["Example1.md", "Example2.md", "Example3.md", "Example4.md", "Example5.md",
        "Example6.md","Example7.md","Example8.md","Example9.md","Example10.md", "Example11.md"]
```

# Examples:


The following sections provide some basic examples of EntropyHub functions. 
These examples are merely a snippet of the full range of EntropyHub functionality. 

In the following examples, signals / data are imported into Julia using the ExampleData() function. 
To use this function as shown in the examples below, __*an internet connection is required*__.


```@docs
EntropyHub.ExampleData
```


!!! tip "IMPORTANT TO NOTE"
    
    Parameters of the base or cross- entropy methods are passed to multiscale and multiscale cross- functions using the multiscale entropy object using MSobject.
    Base and cross- entropy methods are declared with MSobject() using any Base or Cross- entropy function.
    See the MSobject example in the following sections for more info.

!!! warning "Hierarchical Multiscale Entropy (+ Multiscale Cross-Entropy)"
    
    In hierarchical multiscale entropy (hMSEn) and hierarchical multiscale cross-entropy (hXMSEn) functions, the length of the time series signal(s) is halved at each scale. 
    Thus, hMSEn and hXMSEn only use the first 2^N data points where 2^N <= the length of the original time series signal.
    i.e. For a signal of 5000 points, only the first 4096 are used. For a signal of 1500 points, only the first 1024 are used.

!!! danger "BIDIMENSIONAL ENTROPIES"
    
    Each bidimensional entropy function (SampEn2D, FuzzEn2D, DistEn2D) has an important keyword argument - `Lock`. 
    Bidimensional entropy functions are "locked" by default (`Lock == true`) to only permit matrices with a maximum  size of 128 x 128.
