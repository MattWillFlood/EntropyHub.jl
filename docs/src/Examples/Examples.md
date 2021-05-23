```@contents
Pages = ["Example1.md", "Example2.md", "Example3.md", "Example4.md", "Example5.md",
        "Example6.md","Example7.md","Example8.md","Example9.md","Example10.md"]
```

# Examples:


The following sections provide some basic examples of EntropyHub functions. 
These examples are merely a snippet of the full range of EntropyHub functionality. 

In the following examples, signals / data are imported into Julia using the ExampleData() function. 
To use this function as shown in the examples below, __*an internet connection is required*__.


ExampleData() accepts any of the following strings: 

* `"uniform"`         -    vector of uniformly distributed random numbers in range [0 1]
* `"gaussian"`        -    vector of normally distributed random numbers with mean = 0, SD = 1
* `"randintegers"`    -    vector of uniformly distributed pseudorandom integers in range [1 8]
* `"chirp"`           -    vector of chirp signal with the following parameters:   f0 = .01; t1 = 4000; f1 = .025
* `"lorenz"`          -    3-column matrix: X, Y, Z components of the Lorenz system  (alpha: 10; beta: 8/3 ; rho: 28); [Xo = 10; Yo = 20; Zo = 10]
* `"henon"`           -    2-column matrix: X, Y components of the Henon attractor (alpha: 1.4; beta: 0.3); [Xo = 0; Yo = 0]
* `"uniform2"`        -    2-column matrix: uniformly distributed random numbers in range [0 1]
* `"gaussian2"`       -    2-column matrix: normally distributed random numbers with mean = 0, SD = 1
* `"randintegers2"`   -    2-column matrix: uniformly distributed pseudorandom integers in range [1 8]
* `"uniform_Mat"`     -    Matrix of uniformly distributed random numbers in range [0 1]
* `"gaussian_Mat"`    -    Matrix of normally distributed random numbers with mean = 0; SD = 1
* `"randintegers_Mat"`-    Matrix of uniformly distributed pseudorandom integers in range [1 8]
* `"mandelbrot_Mat"`  -    Matrix of image of fractal generated from the mandelbrot set
* `"entropyhub_Mat"`  -    Matrix of image of the entropyhub logo


!!! tip "IMPORTANT TO NOTE"
    
    For cross-entropy and multiscale cross-entropy functions, the two time series signals are passed as a two-column or two-row matrix. 
    At present, it is not possible  to pass signals of different lengths separately. 

    Parameters of the base or cross- entropy methods are passed to multiscale and multiscale cross- functions using the multiscale entropy object using MSobject.
    Base and cross- entropy methods are declared with MSobject() using any [Base](../Guide/Base_Entropies.html) or [Cross-](../Guide/Cross_Entropies.html) entropy function.
    See the MSobject example in the following sections for more info.

!!! warning "Hierarchical Multiscale Entropy (+ Multiscale Cross-Entropy)"
    
    In hierarchical multiscale entropy (hMSEn) and hierarchical multiscale cross-entropy (hXMSEn) functions, the length of the time series signal(s) is halved at each scale. 
    Thus, hMSEn and hXMSEn only use the first 2^N data points where 2^N <= the length of the original time series signal.
    i.e. For a signal of 5000 points, only the first 4096 are used. For a signal of 1500 points, only the first 1024 are used.

!!! danger "BIDIMENSIONAL ENTROPIES"
    
    Each bidimensional entropy function (SampEn2D, FuzzEn2D, DistEn2D) has an important keyword argument - `Lock`. 
    Bidimensional entropy functions are "locked" by default (`Lock == true`) to only permit matrices with a maximum  size of 128 x 128.
