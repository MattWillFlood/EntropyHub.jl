```@meta
CurrentModule = EntropyHub
```
```@raw html
    <div style="display: flex; justify-content: flex-end">    
        <div class="__dimensions_badge_embed__" data-doi="10.1371/journal.pone.0259448"  data-style="small_circle"></div><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
        <div  class="altmetric-embed" data-badge-type='donut' data-badge-popover='right' data-altmetric-id="116252437"></div><script type="text/javascript" src="https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js"></script>
        <script data-name="BMC-Widget" data-cfasync="false" src="https://cdnjs.buymeacoffee.com/1.0.0/widget.prod.min.js" data-id="EntropyHub" data-description="Support me on Buy me a coffee!" data-message="Want to support EntropyHub? You can support the project by buying us a coffee!" data-color="#FF813F" data-position="Right" data-x_margin="18" data-y_margin="18"></script>
    </div>
```

![EH4J](./assets/logo.png)  


# EntropyHub
__*An Open-Source Toolkit For Entropic Time Series Analysis*__

---

#### `EntropyHub.jl` is part of the EntropyHub project.
#### For more info visit: [www.EntropyHub.xyz](https://www.EntropyHub.xyz)

!!! tip ""

    Also available with: [Matlab](https://www.mathworks.com/matlabcentral/fileexchange/94185-entropyhub) // [Python](https://pypi.org/project/EntropyHub/) 

---

## Latest Updates
### v2.0
*----- New multivariate methods -----*

Five new multivariate entropy functions incorporating several method-specific variations
   + [Multivariate Sample Entropy](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.84.061918)       
   + [Multivariate Fuzzy Entropy](https://www.mdpi.com/1099-4300/19/1/2) [++ many fuzzy functions]              
   + [Multivariate Dispersion Entropy](https://www.mdpi.com/1099-4300/21/9/913) [++ many symbolic sequence transforms]         
   + [Multivariate Cosine Similarity Entropy](https://www.mdpi.com/1099-4300/24/9/1287)               
   + Multivariate Permutation Entropy  [++ *amplitude-aware*, *edge*, *phase*, *weighted* and *modified* variants]              

*----- New multivariate multiscale methods -----*

Two new multivariate multiscale entropy functions

   + [Multivariate Multiscale Entropy](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.84.061918) [++ coarse, modified and generalized graining procedures]                
   + [Composite and Refined-composite Multivariate Multiscale Entropy](https://link.springer.com/article/10.1007/s11517-017-1647-5)             

*----- Extra signal processing tools -----*  

[`WindowData()`](@ref) is a new function that allows users to segment data (univariate or multivariate time series) into windows with/without overlapping samples! This allows users to calculate entropy on subsequences of their data to perform analyses with greater time resolution.        

__**Other little fixes...**__
   
*----- Docs edits -----*

   - Examples in the www.EntropyHub.xyz documentation were updated to match the latest package syntax.    
        
---

## Introduction

This toolkit provides a wide range of functions to calculate different entropy statistics. 
There is an ever-growing range of information-theoretic and dynamical systems entropy measures presented in the scientific literature. The goal of **EntropyHub.jl** is to integrate the many established entropy methods in one open-source package with an extensive documentation and consistent syntax [that is also accessible in multiple programming languages ([Matlab](https://www.entropyhub.xyz/matlab/EHmatlab.html), [Python](https://www.entropyhub.xyz/python/EHpython.html))].


### About 

Information and uncertainty can be regarded as two sides of the same coin: 
the more uncertainty there is, the more information we gain by removing that 
uncertainty. In the context of information and probability theory, **Entropy** 
quantifies that uncertainty. 

Various measures have been derived 
to estimate entropy (uncertainty) from discrete time series, each seeking to 
best capture the uncertainty of the system under examination. This has resulted 
in many entropy statistics from approximate entropy and sample entropy, to
multiscale sample entropy and refined-composite multiscale cross-sample entropy.

The goal of EntropyHub is to provide a comprehensive set of functions with a simple and 
consistent syntax that allows the user to augment parameters at the command 
line, enabling a range from basic to advanced entropy methods to be implemented
with ease.

!!! warning "NOTE:"
    It is important to clarify that the entropy functions herein described estimate entropy 
    in the context of probability theory and information theory as defined by Shannon, 
    and not thermodynamic or other entropies from classical physics.

## Installation
Using the Julia REPL:

```
julia> using Pkg; Pkg.add("EntropyHub")
```

or

```
julia> ] 
pkg> add EntropyHub
```

To get the latest version of EntropyHub directly from GitHub:
```
julia> ] 
pkg> add https://github.com/MattWillFlood/EntropyHub.jl
```

## Citing
EntropyHub is licensed under the Apache License (Version 2.0) and is free to use by
all on condition that the following reference be included on any outputs realized using the
software:
```
Matthew W. Flood (2021),
EntropyHub: An Open-Source Toolkit for Entropic Time Series Analysis,
PLoS ONE 16(11):e0259448
DOI:  10.1371/journal.pone.0259448
www.EntropyHub.xyz 
```
---

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

If you find this package useful, please consider starring it on [GitHub](https://github.com/MattWillFlood/EntropyHub.jl) 
and Julia Packages (or MatLab File Exchange and PyPI). This helps us to gauge user satisfaction.

## Contact

For general queries and information about EntropyHub, contact: `info@entropyhub.xyz`   

If you have any questions or need help using the package, please contact us at: `help@entropyhub.xyz`     

If you notice or identify any issues, please do not hesitate to contact us at: `fix@entropyhub.xyz`     


We will do our best to help you with any relevant issues that you may have.

If you come across any errors or technical issues, you can raise these under the issues tab on the EntropyHub.jl GitHub page. 
Similarly, if you have any suggestions or recommendations on how this package can be improved, please let us know.

__Thank you for using EntropyHub,__

Matt

## Function List

EntropyHub functions fall into 8 categories: 

* `Base`                           functions for estimating the entropy of a single univariate time series.
* `Cross`                          functions for estimating the entropy between two univariate time series.
* `Bidimensional`                  functions for estimating the entropy of a two-dimensional univariate matrix.
* `Multiscale`                     functions for estimating the multiscale entropy of a single univariate time series using any of the Base entropy functions.
* `Multiscale Cross`               functions for estimating the multiscale entropy between two univariate time series using any of the Cross-entropy functions.
* `Multivariate Multiscale`        functions for estimating the multivariate multiscale entropy of multivariate dataset using any of the Multivariate-entropy functions.
* `Other`                          Supplementary functions for various tasks related to EntropyHub and signal processing.


```@index
```


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




Documentation for [EntropyHub](https://github.com/MattWillFlood/EntropyHub.jl).