# EntropyHub: An open-source toolkit for entropic data analysis
__*Julia Edition*__

<p align="center">
<img src="https://github.com/MattWillFlood/EntropyHub/blob/main/Graphics/EntropyHub_JuliaLogo.png" alt = "EntropyHub Git" width="250" height="340" />
</p>


## Latest Update
### v1.0
__*----- New entropy methods -----*__             
Two new base entropy functions (and their multiscale versions) have been added:         
        > [Diversity Entropy](https://ieeexplore.ieee.org/document/9194995)             
        > [Range Entropy](https://www.mdpi.com/1099-4300/20/12/962)             

__*----- New fuzzy membership functions -----*__          
Several new fuzzy membership functions have been added to FuzzEn, XFuzzEn and FuzzEn2D to provide more options for mapping the degree of similarity between embedding vectors.          
These include *trapezoidal*, *triangular* and *gaussian*, among others.                 
Further info on these membership functions can be found [here.](https://hal.science/hal-02267711/document)              

__*----- Phase Permutation Entropy -----*__               
A new variant - '*phase*' permutation entropy - has been added to PermEn.                  
This method employs a hilbert transformation of the data sequence, based on the methods outlined [here.](https://doi.org/10.1016/j.physa.2020.125686)           

__*----- Cross-Entropy with different length sequences -----*__           
EntropyHub v1.0 now allows for cross-entropy (and multiscale cross-entropy) estimation with different length signals (_except XCondEn and XPermEn_).            
As a result, the new cross-entropy functions require a separate input for each sequence (Sig1, Sig2).                   

__*----- Refined-Composite Multiscale Fuzzy Entropy -----*__              
In addition to the refined-composite multiscale sample entropy that was available in earlier versions, now one can estimate the refined-composite multiscale fuzzy entropy based on the method outlined [here.](https://link.springer.com/article/10.1007/s11517-017-1647-5)    
What's more, refined-composite multicale cross-fuzzy entropy is also available, and both can be estimated using any of the fuzzy membership functions in FuzzEn or XFuzzEn.             

__*----- Generalized Multiscale Entropy -----*__          
Generaized multiscale entropy and generalized multiscale cross-entropy can now be estimated. Just choose the '*generalized*' as the graining procedure in MSEn or XMSEn.        

__*----- Variance of sample entropy estimate -----*__             
Based on the [method outlined by Lake et al.,](https://journals.physiology.org/doi/epdf/10.1152/ajpregu.00069.2002) it is now possible to obtain a measure of the variance in the sample entropy estimate.              
This is achieved by approximating the number of overlapping embedding vectors.          
To do so, just set the parameter '*Vcp*'==true in SampEn and XSampEn, but note that doing so requires *a lot* of computer memory.              

*Several little bugs and inconsistencies have also been fixed in this release. We want to thank all of you who have identified and alerted us to these bugs.*       
*Most of these bugs have been noted via the [GitHub issues portal](https://github.com/MattWillFlood/EntropyHub/issues).*        

__*----- Bug fixes -----*__             
        - The DispEn2D function in python has now fixed [this issue](https://github.com/MattWillFlood/EntropyHub/issues/8).     
        - The type hint for FuzzEn in python has been updated [as requested](https://github.com/MattWillFlood/EntropyHub/issues/1).             
        - [Compatbility issues with EntropyHub.jl](https://github.com/MattWillFlood/EntropyHub.jl/issues/3) are now resolved.           
        - A bug in the K2En python function led to incorrect entropy estimates for data sequences with many equal values. This has been corrected.              
    
__*----- Other Changes -----*__             
        - The *'equal'* method for discretizing data in DispEn and DispEn2D has been updated to be consistent across Python, MatLab and Julia. This is unlikely to have impacted any users previously.          
        - The zeroth dimension (m=0) estimate of ApEn and XApEn has been changed to -phi(1).                    
        - The default radius threshold distance for XApEn, XSampEn and XK2En has been changed to use the *pooled* standard deviation [i.e. 0.2*SDpooled(X,Y)].   
        

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
Phase Entropy                                        |	PhasEn
Slope Entropy                                        |	SlopEn
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
Cross Spectral Entropy                          	   |	XSpecEn
Cross Kolmogorov Entropy                              |	XK2En
	
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

