```@meta
CurrentModule = EntropyHub
```

![EH4J](./assets/logo.png)


# EntropyHub
__*An Open-Source Toolkit For Entropic Time Series Analysis*__

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

If you find this package useful, please consider starring it on [GitHub](https://github.com/MattWillFlood/EntropyHub.jl) 
and Julia Packages (or MatLab File Exchange and PyPI). This helps us to gauge user satisfaction.

## Function List

EntropyHub functions fall into 5 categories: 

* `Base`                functions for estimating the entropy of a single univariate time series.
* `Cross`               functions for estimating the entropy between two univariate time series.
* `Bidimensional`       functions for estimating the entropy of a two-dimensional univariate matrix.
* `Multiscale`          functions for estimating the multiscale entropy of a single univariate time series using any of the Base entropy functions.
* `Multiscale Cross`    functions for estimating the multiscale entropy between two univariate time series using any of the Cross-entropy functions.

## Contact

For general queries and information about EntropyHub, contact: `info@entropyhub.xyz`   

If you have any questions or need help using the package, please contact us at: `help@entropyhub.xyz`     

If you notice or identify any issues, please do not hesitate to contact us at: `fix@entropyhub.xyz`     


We will do our best to help you with any relevant issues that you may have.

If you come across any errors or technical issues, you can raise these under the issues tab on the EntropyHub.jl GitHub page. 
Similarly, if you have any suggestions or recommendations on how this package can be improved, please let us know.

__Thank you for using EntropyHub,__

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




Documentation for [EntropyHub](https://github.com/MattWillFlood/EntropyHub.jl).

```@index
```
