# EntropyHub: An open-source toolkit for entropic time series analysis
__*Julia Edition*__

<p align="center">
<img src="https://github.com/MattWillFlood/EntropyHub/blob/main/Graphics/EntropyHub_JuliaLogo.png" alt = "EntropyHub Git" width="250" height="340" />
</p>

## About

Information and uncertainty can be regarded as two sides of the same coin: 
the more uncertainty there is, the more information we gain by removing that 
uncertainty. In the context of information and probability theory, ***Entropy*** 
quantifies that uncertainty. 

The concept of entropy has its origins in 
[classical physics](http://www.scholarpedia.org/article/Entropy "Scholarpedia")
under the second law of thermodynamics, a law 
[considered to underpin our fundamental understanding](https://www.penguin.co.uk/books/301539/the-order-of-time/9780141984964.html "Rovelli") 
of [time in physics](https://en.wikipedia.org/wiki/Time_in_physics "Wiki Time"). 
Attempting to analyse the analog world around
us requires that we measure time in discrete steps, but doing so compromises 
our ability to measure entropy accurately. Various measures have been derived 
to estimate entropy (uncertainty) from discrete time series, each seeking to 
best capture the uncertainty of the system under examination. This has resulted 
in many entropy statistics from approximate entropy and sample entropy, to
multiscale sample entropy and refined-composite multiscale cross-sample entropy.


As the number of statisitcal entropy measures grows, it becomes more difficult
to identify, contrast and compare the performance of each measure. To overcome
this, we have developed EntropyHub - an open-source toolkit designed to 
integrate the many established entropy methods into one package. The goal of 
EntropyHub is to provide a comprehensive set of functions with a simple and 
consistent syntax that allows the user to augment parameters at the command 
line, enabling a range from basic to advanced entropy methods to be implemented
with ease.

***It is important to clarify that the entropy functions herein described 
estimate entropy in the context of probability theory and information theory as
defined by Shannon, and not thermodynamic or other entropies from classical physics.***

## Installation

There are two ways to install EntropyHub for Julia.

#### Method 1:
   1. In Julia, open the package REPL by typing `]`. The command line should appear as: 
   
      `@vX.Y.  pkg>        `
      
      Where X and Y refer to your version of Julia.
  
   2. Type:
   
      `add EntropyHub`
	
      (Note: this is __case sensitive__)
      
   Alternatively, one can use the Pkg module to perform the same procedure:
   
   `using Pkg`
	
   `Pkg.add("EntropyHub")`
    
#### Method 2:
   1. In Julia, open the package REPL by typing `]`. The command line should appear as: 
   
      `@vX.Y.  pkg>        `
      
      Where X and Y refer to your version of Julia.
  
   2. Type:
   
      `add https://github.com/MattWillFlood/EntropyHub.jl`
	
      (Note: this is __case sensitive__)
      


### System Requirements 
There are several package dependencies which will be installed alongside EntropyHub (if not already installed):

DSP, FFTW, HTTP, Random, Plots, StatsBase, StatsFuns, GroupSlices, 
Statistics, DelimitedFiles, Combinatorics, LinearAlgebra, Dierckx, Clustering

EntropyHub was designed using Julia 1.5 and is intended for use with Julia versions >= 1.2.


## Documentation & Help 

A key advantage of EntropyHub is the comprehensive documentation available to help users to make the most of the toolkit.
 
To learn more about a specific function, one can do so easily from the command line by typing: `?`, which will open the julia help system, and then typing the function name.

For example:

	julia> ?  
	help?> SampEn	  # Documentation on sample entropy function

	julia> ?  
	help?> XSpecEn    # Documentation on cross-spectral entropy function

	julia> ?
	help?> hXMSEn     # Documentation on hierarchical multiscale cross-entropy function

All information on the EntropyHub package is detailed in the *EntropyHub Guide*, a .pdf document available [here](https://github.com/MattWillFlood/EntropyHub/blob/main/EntropyHub%20Guide.pdf).


## Functions

EntropyHub functions fall into 5 categories: 

    * Base                functions for estimating the entropy of a single univariate time series.
    * Cross               functions for estimating the entropy between two univariate time series.
    * Bidimensional       functions for estimating the entropy of a two-dimensional univariate matrix.
    * Multiscale          functions for estimating the multiscale entropy of a single univariate time series using any of the Base entropy functions.
    * Multiscale Cross    functions for estimating the multiscale entropy between two univariate time series using any of the Cross-entropy functions.

#### The following tables outline the functions available in the EntropyHub package.

*When new entropies are published in the scientific literature, efforts will
be made to incorporate them in future releases.*

### Base Entropies:

Entropy Type   |  Function Name 
---|---
Approximate Entropy                               	  |	ApEn
Sample Entropy                                		  |	SampEn
Fuzzy Entropy                                 		  |	FuzzEn
Kolmogorov Entropy                            		  |	K2En
Permutation Entropy                           		  |	PermEn
Conditional Entropy                           		  |	CondEn
Distribution Entropy                          		  |	DistEn
Spectral Entropy                              		  |	SpecEn
Dispersion Entropy                            		  |	DispEn
Symbolic Dynamic Entropy                          	  |	SyDyEn
Increment Entropy                                 	  |	IncrEn
Cosine Similarity Entropy                         	  |	CoSiEn
Phase Entropy                                         |	PhasEn
Slope Entropy                                      	  |	SlopEn
Bubble Entropy                                		  |	BubbEn
Gridded Distribution Entropy                          |	GridEn
Entropy of Entropy                            	      |	EnofEn
Attention Entropy                                     |	AttnEn

_______________________________________________________________________

### Cross Entropies:

Entropy Type   |  Function Name 
--|--
Cross Sample Entropy                                  |	XSampEn
Cross Approximate Entropy                             |	XApEn
Cross Fuzzy Entropy                                   |	XFuzzEn
Cross Permutation Entropy                             |	XPermEn
Cross Conditional Entropy                             |	XCondEn
Cross Distribution Entropy                            |	XDistEn
Cross Spectral Entropy                          	  |	XSpecEn
Cross Kolmogorov Entropy                              |	XK2En
	
_______________________________________________________________________

### Bidimensional Entropies

Entropy Type   |  Function Name 
--|--
Bi-Dimensional Sample Entropy                         |	SampEn2D
Bi-Dimensional Fuzzy Entropy                          |	FuzzEn2D
Bi-Dimensional Distribution Entropy                   |	DistEn2D
	
_________________________________________________________________________

### Multiscale Entropy Functions

Entropy Type   |  Function Name 
--|--
Multiscale Entropy                                    | MSEn
Composite/Refined-Composite Multiscale Entropy        | cMSEn
Refined Multiscale Entropy                            | rMSEn
Hierarchical Multiscale Entropy                 | hMSEn
	
_________________________________________________________________________

### Multiscale Cross-Entropy Functions
Entropy Type   |  Function Name 
--|--
Multiscale Cross-Entropy                              |   XMSEn
Composite/Refined-Composite Multiscale Cross-Entropy  |   cXMSEn
Refined Multiscale Cross-Entropy                            |   rXMSEn
Hierarchical Multiscale Cross-Entropy                 |   hXMSEn





## License and Terms of Use
EntropyHub is licensed under the Apache License (Version 2.0) and is free to
use by all on condition that the following reference be included on any outputs
realized using the software:
 
        Matthew W. Flood and Bernd Grimm, 
        EntropyHub: An Open-Source Toolkit for Entropic Time Series Analysis,
        2021 www.EntropyHub.xyz

__________________________________________________________________


        © Copyright 2021 Matthew W. Flood, EntropyHub
        Licensed under the Apache License, Version 2.0 (the "License");
        you may not use this file except in compliance with the License.
        You may obtain a copy of the License at
        
                 http://www.apache.org/licenses/LICENSE-2.0
        
        Unless required by applicable law or agreed to in writing, software
        distributed under the License is distributed on an "AS IS" BASIS,
        WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
        See the License for the specific language governing permissions and
        limitations under the License.
        
        For Terms of Use see https://github.com/MattWillFlood/EntropyHub



## Contact

If you find this package useful, please consider starring it on GitHub, 
MatLab File Exchange, PyPI or Julia Packages as this helps us to gauge user 
satisfaction.

If you have any questions about the package or identify any issues, 
please do not hesitate to contact us at:    info@entropyhub.xyz


***Thank you*** for using EntropyHub.

Yours in research,

Matt

        
        
        
        
         ___  _   _  _____  _____  ____  ____  _     _          
        |  _|| \ | ||_   _||     \|    ||    || \   / |   ___________ 
        | \_ |  \| |  | |  |   __/|    ||  __| \ \_/ /   /  _______  \
        |  _|| \ \ |  | |  |   \  |    || |     \   /   |  /  ___  \  |
        | \_ | |\  |  | |  | |\ \ |    || |      | |    | |  /   \  | | 
        |___||_| \_|  |_|  |_| \_||____||_|      |_|   _|_|__\___/  | | 
         _   _  _   _  ____                           / |__\______\/  | 
        | | | || | | ||    \     An open-source      |  /\______\__|_/ 
        | |_| || | | ||    |     toolkit for         | |  /   \  | | 
        |  _  || | | ||    \     entropic time-      | |  \___/  | |          
        | | | || |_| ||     \    series analysis     |  \_______/  |
        |_| |_|\_____/|_____/                         \___________/ 
        
        
<p  align="center">
	<img src="https://github.com/MattWillFlood/EntropyHub/blob/main/Graphics/EntropyHubLogo3.png" width="250" height="350"/>
</p>

