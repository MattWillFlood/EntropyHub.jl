# EntropyHub: An open-source toolkit for entropic data analysis
__*Julia Edition*__

<p align="center">
<img src="https://github.com/MattWillFlood/EntropyHub/blob/main/Graphics/EntropyHub_JuliaLogo.png" alt = "EntropyHub Git" width="250" height="340" />
</p>


## Latest Update
### v2.0
__*----- New multivariate methods -----*__       
Five new multivariate entropy functions incorporating several method-specific variations        
        > [Multivariate Sample Entropy](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.84.061918)      
        > [Multivariate Fuzzy Entropy](https://www.mdpi.com/1099-4300/19/1/2) [++ many fuzzy functions]        
        > [Multivariate Dispersion Entropy](https://www.mdpi.com/1099-4300/21/9/913) [++ many symbolic sequence transforms]          
        > [Multivariate Cosine Similarity Entropy](https://www.mdpi.com/1099-4300/24/9/1287)        
        > Multivariate Permutation Entropy  [++ *amplitude-aware*, *edge*, *phase*, *weighted* and *modified* variants]       

__*----- New multivariate multiscale methods -----*__       
Two new multivariate multiscale entropy functions        
        > [Multivariate Multiscale Entropy](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.84.061918) [++ coarse, modified and generalized graining procedures]          
        > [Composite and Refined-composite Multivariate Multiscale Entropy](https://link.springer.com/article/10.1007/s11517-017-1647-5)      

__*----- Extra signal processing tools -----*__     
**WindowData()** is a new function that allows users to segment data (univariate or multivariate time series) into windows with/without overlapping samples! This allows users to calculate entropy on subsequences of their data to perform analyses with greater time resolution.        

*Other little fixes...*
   
__*----- Docs edits -----*__     
        - Examples in the www.EntropyHub.xyz documentation were updated to match the latest package syntax.     
        

## Welcome
This toolkit provides a wide range of functions to calculate different entropy statistics.
There is an ever-growing range of information-theoretic and dynamical systems entropy measures presented in the scientific literature. 
The goal of EntropyHub is to integrate the many established entropy methods in one open-source package.


## About

Information and uncertainty can be regarded as two sides of the same coin: 
the more uncertainty there is, the more information we gain by removing that 
uncertainty. In the context of information and probability theory, ***Entropy*** 
quantifies that uncertainty. 
Attempting to analyse the analog world around
us requires that we measure time in discrete steps, but doing so compromises 
our ability to measure entropy accurately. Various measures have been derived 
to estimate entropy (uncertainty) from discrete time series, each seeking to 
best capture the uncertainty of the system under examination. This has resulted 
in many entropy statistics from approximate entropy and sample entropy, to
multiscale sample entropy and refined-composite multiscale cross-sample entropy.

The goal of EntropyHub is to provide a comprehensive set of functions with a simple and 
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
Statistics, DelimitedFiles, Combinatorics, LinearAlgebra, DataInterpolations, Clustering

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

EntropyHub functions fall into 8 categories: 

    * Base                       functions for estimating the entropy of a single univariate time series.
    * Cross                      functions for estimating the entropy between two univariate time series.
    * Multivariate               functions for estimating the entropy of a multivariate dataset.
    * Bidimensional              functions for estimating the entropy of a two-dimensional univariate matrix.
    * Multiscale                 functions for estimating the multiscale entropy of a single univariate time series using any of the Base entropy functions.
    * Multiscale Cross           functions for estimating the multiscale entropy between two univariate time series using any of the Cross-entropy functions.
    * Multivariate Multiscale    functions for estimating the multivariate multiscale entropy of multivariate dataset using any of the Multivariate-entropy functions.
    * Other                      Supplementary functions for various tasks related to EntropyHub and signal processing.

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
Phase Entropy                                        |	PhasEn
Slope Entropy                                    	  |	SlopEn
Bubble Entropy                                		  |	BubbEn
Gridded Distribution Entropy                         |	GridEn
Entropy of Entropy                            	     |	EnofEn
Attention Entropy                                    |	AttnEn
Range Entropy                                        |   RangEn
Diversity Entropy                                    |   DivEn

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
Cross Spectral Entropy                            	   |	XSpecEn
Cross Kolmogorov Entropy                              |	XK2En
	
_______________________________________________________________________

### Multivariate Entropies:

Entropy Type   |  Function Name 
--|--
Multivariate Sample Entropy                                  |	MvSampEn
Multivariate Fuzzy Entropy                                   |	MvFuzzEn
Multivariate Permutation Entropy                             |	MvPermEn
Multivariate Dispersion Entropy                              |	MvDispEn
Multivariate Cosine Similarity Entropy                       |	MvCoSiEn

_______________________________________________________________________

### Bidimensional Entropies

Entropy Type   |  Function Name 
--|--
Bidimensional Sample Entropy                         |	SampEn2D
Bidimensional Fuzzy Entropy                          |	FuzzEn2D
Bidimensional Distribution Entropy                   |	DistEn2D
Bidimensional Dispersion Entropy                     |	DispEn2D
Bidimensional Permutation Entropy                    |	PermEn2D
Bidimensional Espinosa Entropy                       |	EspEn2D
	
_________________________________________________________________________

### Multiscale Entropy Functions

Entropy Type   |  Function Name 
--|--
Multiscale Entropy                                    | MSEn
Composite/Refined-Composite Multiscale Entropy        | cMSEn
Refined Multiscale Entropy                            | rMSEn
Hierarchical Multiscale Entropy                       | hMSEn
	
_________________________________________________________________________

### Multiscale Cross-Entropy Functions
Entropy Type   |  Function Name 
--|--
Multiscale Cross-Entropy                              |   XMSEn
Composite/Refined-Composite Multiscale Cross-Entropy  |   cXMSEn
Refined Multiscale Cross-Entropy                      |   rXMSEn
Hierarchical Multiscale Cross-Entropy                 |   hXMSEn

_________________________________________________________________________

### Multivariate Multiscale Entropy Functions

Entropy Type   |  Function Name 
--|--
Multivariate Multiscale Entropy                                    | MvMSEn
Composite/Refined-Composite Multivariate Multiscale Entropy        | cMvMSEn

_________________________________________________________________________

### Other Functions

Entropy Type   |  Function Name 
--|--
Example Data Import Tool            |  ExampleData
Window Data Tool                    |  WindowData
Multiscale Entropy Object           |  MSobject


## License and Terms of Use
EntropyHub is licensed under the Apache License (Version 2.0) and is free to
use by all on condition that the following reference be included on any outputs
realized using the software:
 
        Matthew W. Flood (2021), 
        EntropyHub: An Open-Source Toolkit for Entropic Time Series Analysis,
        PLoS ONE 16(11):e0259448
        DOI:   10.1371/journal.pone.0259448
        www.EntropyHub.xyz

__________________________________________________________________


        © Copyright 2024 Matthew W. Flood, EntropyHub
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

If you have any questions about the package, please do not hesitate to contact us at: info@entropyhub.xyz
If you identify any bugs, please contact us at: fix@entropyhub.xyz
If you need any help installing or using the toolkit, please contact us at: help@entropyhub.xyz


***Thank you*** for using EntropyHub.

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