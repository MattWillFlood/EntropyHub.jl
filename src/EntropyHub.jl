module EntropyHub

export 
    ApEn,
    AttnEn, 
    BubbEn, 
    CoSiEn, 
    CondEn, 
    DispEn, 
    DistEn,
    DivEn,
    EnofEn,
    FuzzEn,
    GridEn,
    IncrEn,
    K2En,
    PermEn,
    PhasEn,
    RangEn,
    SampEn,
    SlopEn,
    SyDyEn,
    SpecEn,

    XApEn,
    XCondEn,
    XDistEn,
    XFuzzEn,
    XK2En,
    XPermEn,
    XSpecEn,
    XSampEn,

    SampEn2D,
    FuzzEn2D,
    DistEn2D,
    DispEn2D,
    PermEn2D,
    EspEn2D,

    MSobject,
    MSEn,
    cMSEn,
    rMSEn,
    hMSEn,
    XMSEn,
    cXMSEn,
    rXMSEn,
    hXMSEn,

    EMD,
    ExampleData

# using Reexport

include("./_ApEn.jl"),
include("./_AttnEn.jl"),
include("./_BubbEn.jl"),
include("./_CondEn.jl"),
include("./_CoSiEn.jl"),
include("./_DispEn.jl"),
include("./_DistEn.jl"),
include("./_DivEn.jl"),
include("./_FuzzEn.jl"),
include("./_EnofEn.jl"),
include("./_GridEn.jl"),
include("./_IncrEn.jl"),
include("./_K2En.jl"),
include("./_PermEn.jl"),
include("./_PhasEn.jl"),
include("./_RangEn.jl"),
include("./_SampEn.jl"),
include("./_SpecEn.jl"),
include("./_SyDyEn.jl"),
include("./_SlopEn.jl"),

include("./_XApEn.jl"),
include("./_XFuzzEn.jl"),
include("./_XDistEn.jl"),
include("./_XCondEn.jl"),
include("./_XSpecEn.jl"),
include("./_XSampEn.jl"),
include("./_XK2En.jl"),
include("./_XPermEn.jl"),

include("./_MSobject.jl"),
include("./_ExampleData.jl"),

include("./_MSEn.jl"),
include("./_XMSEn.jl"),
include("./_cMSEn.jl"),
include("./_cXMSEn.jl"),
include("./_rMSEn.jl"),
include("./_rXMSEn.jl"),
include("./_hMSEn.jl"),
include("./_hXMSEn.jl"),

include("./_SampEn2D.jl"),
include("./_FuzzEn2D.jl"),
include("./_DistEn2D.jl"),
include("./_DispEn2D.jl"),
include("./_PermEn2D.jl"),
include("./_EspEn2D.jl"),

#= 
@reexport using ._ApEn
@reexport using ._MSobject

@reexport using ._AttnEn
@reexport using ._BubbEn
@reexport using ._CondEn
@reexport using ._CoSiEn
@reexport using ._DistEn
@reexport using ._DispEn
@reexport using ._EnofEn
@reexport using ._FuzzEn =#


# Base Entropies:
using ._ApEn: ApEn
using ._AttnEn: AttnEn
using ._BubbEn: BubbEn
using ._CoSiEn: CoSiEn
using ._CondEn: CondEn
using ._DispEn: DispEn
using ._DistEn: DistEn
using ._DivEn: DivEn
using ._EnofEn: EnofEn
using ._FuzzEn: FuzzEn
using ._GridEn: GridEn
using ._IncrEn: IncrEn
using ._K2En: K2En
using ._PermEn: PermEn
using ._PhasEn: PhasEn
using ._RangEn: RangEn
using ._SampEn: SampEn
using ._SlopEn: SlopEn
using ._SpecEn: SpecEn
using ._SyDyEn: SyDyEn

# Cross Entropies:
using ._XApEn: XApEn
using ._XCondEn: XCondEn
using ._XDistEn: XDistEn
using ._XFuzzEn: XFuzzEn
using ._XK2En: XK2En
using ._XPermEn: XPermEn
using ._XSampEn: XSampEn
using ._XSpecEn: XSpecEn

# Bidimensional Entropies
using ._SampEn2D: SampEn2D
using ._DistEn2D: DistEn2D
using ._FuzzEn2D: FuzzEn2D
using ._DispEn2D: DispEn2D
using ._PermEn2D: PermEn2D
using ._EspEn2D: EspEn2D

# (cross) Multiscale Entropies
using ._MSobject: MSobject

using ._MSEn: MSEn, EMD
using ._cMSEn: cMSEn
using ._rMSEn: rMSEn
using ._hMSEn: hMSEn

using ._XMSEn: XMSEn
using ._cXMSEn: cXMSEn
using ._rXMSEn: rXMSEn
using ._hXMSEn: hXMSEn

using ._ExampleData: ExampleData


greet() = print(raw"""  
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
    
    
    Please use the following citation on any scientific outputs achieved with the help of EntropyHub:

    Matthew W. Flood,
    EntropyHub: An Open-Source Toolkit for Entropic Time Series Analysis,
    PLoS One 16(11):e0259448 (2021),
    DOI: 10.1371/journal.pone.0259448

    www.EntropyHub.xyz
    
    """)



raw"""
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

EntropyHub functions belong to one of five main classes/categories:
    Base Entropies             >>  e.g. Approximate Entropy (ApEn),
                                        Sample Entropy (SampEn)
    Cross Entropies            >>  e.g. Cross-Approximate Entropy (XApEn)
                                        Cross-Sample Entropy (XSampEn)
    Bidimensional Entropies    >>  e.g. Bidimensional Sample Entropy (SampEn2D)
                                        Bidimensional Fuzzy Entropy (FuzzEn2D)
    Multiscale Entropies       >>  e.g. Multiscale Sample Entropy (MSEn)
                                        Refined Multiscale Sample Entropy (rMSEn)
                                        Composite Multiscale Sample Entropy (cMSEn)
    Multiscale Cross Entropies >>  e.g. Multiscale Cross-Sample Entropy (XMSEn)
                                        Refined Multiscale Cross-Sample Entropy (rXMSEn)

_________________________________________________________________________
Base Entropies                                        |	Function Name	
______________________________________________________|__________________
Approximate Entropy                               	  |	ApEn
Sample Entropy                                		    |	SampEn
Fuzzy Entropy                                 		    |	FuzzEn
Kolmogorov Entropy                            		    | K2En
Permutation Entropy                           		    |	PermEn
Conditional Entropy                           		    |	CondEn
Distribution Entropy                          		    |	DistEn
Spectral Entropy                              		    |	SpecEn
Dispersion Entropy                            		    |	DispEn
Symbolic Dynamic Entropy                          	  |	SyDyEn
Increment Entropy                                 	  |	IncrEn
Cosine Similarity Entropy                         	  |	CoSiEn
Phase Entropy                                         |	PhasEn
Slope Entropy                                      	  |	SlopEn
Bubble Entropy                                		    |	BubbEn
Gridded Distribution Entropy                          |	GridEn
Entropy of Entropy                            		    |	EnofEn
Attention Entropy                                     |	AttnEn
Diversity Entropy                                     | DivEn
Range Entropy                                         | RangEn

_________________________________________________________________________
Cross Entropies                                       |	Function Name
______________________________________________________|__________________
Cross Sample Entropy                                  |	XSampEn
Cross Approximate Entropy                             |	XApEn
Cross Fuzzy Entropy                                   |	XFuzzEn
Cross Permutation Entropy                             |	XPermEn
Cross Conditional Entropy                             |	XCondEn
Cross Distribution Entropy                            |	XDistEn
Cross Spectral Entropy                                |	XSpecEn
Cross Kolmogorov Entropy                              |	XK2En
	
_________________________________________________________________________
Bidimensional Entropies                              |	Function Name
_____________________________________________________|__________________
Bidimensional Sample Entropy                         |	SampEn2D
Bidimensional Fuzzy Entropy                          |	FuzzEn2D
Bidimensional Distribution Entropy                   |	DistEn2D
Bidimensional Dispersion Entropy                     |	DispEn2D
Bidimensional Permutation Entropy                    |	PermEn2D
Bidimensional Espinosa Entropy                       |	EspEn2D
	
_________________________________________________________________________
Multiscale Entropy Functions                          | Function Name
______________________________________________________|__________________ 
Multiscale Entropy Object                             |   MSobject
                                                      |
Multiscale Entropy                                    |   MSEn
Composite/Refined-Composite Multiscale Entropy        |   cMSEn
Refined Multiscale Entropy                            |   rMSEn
Hierarchical Multiscale Entropy Object                |   hMSEn

_________________________________________________________________________
Multiscale Entropies	MSEn                            |	Function Name
_________________________________________________________________________
Multiscale Sample Entropy                             |	
Multiscale Approximate Entropy                        |	
Multiscale Fuzzy Entropy                              |	
Multiscale Permutation Entropy                        |	
Multiscale Dispersion Entropy                         |	
Multiscale Cosine Similarity Entropy                  |	
Multiscale Symblic Dynamic Entropy                    |	  MSobject
Multiscale Conditional Entropy                        |	     +
Multiscale Entropy of Entropy                         |  MSEn / cMSEn
Multiscale Gridded Distribution Entropy               |	rMSEn / hMSEn
Multiscale Slope Entropy                              |
Multiscale Phase Entropy                              |		
Multiscale Kolmogorov Entropy                         |	
Multiscale Distribution Entropy                       |	
Multiscale Bubble Entropy                             |	
Multiscale Increment Entropy			                    |	
Multiscale Attention Entropy                          |	
Multiscale Diversity Entropy                          | 
Multiscale Range Entropy                              | 
	
_________________________________________________________________________
Multiscale Cross-Entropy Functions                    |   Function Name
______________________________________________________|__________________
Multiscale Cross-Entropy Object                       |   MSobject
                                                      |
Multiscale Cross-Entropy                              |   XMSEn
Composite/Refined-Composite Multiscale Cross-Entropy  |   cXMSEn
Refined Multiscale Entropy                            |   rXMSEn
Hierarchical Multiscale Entropy Object                |   hXMSEn

_________________________________________________________________________
Multiscale Cross-Entropies                            |	Function Name
_________________________________________________________________________
Multiscale Cross-Sample Entropy                       |	
Multiscale Cross-Approximate Entropy                  |	
Multiscale Cross-Fuzzy Entropy                        |	MSobject
Multiscale Cross-Permutation Entropy                  |	    +
Multiscale Cross-Distribution Entropy                 |	XMSEn / cXMSEn
Multiscale Cross-Kolmogorov Entropy                   |	rXMSEn / hXMSEn
Multiscale Cross-Conditional Entropy                  |	


  We kindly ask that if you use EntropyHub in your research, to please 
  include the following citation with the appropriate version number,
  as well as original articles upon which functions are derived:

  Matthew W. Flood (2021), 
  "EntropyHub - An open source toolkit for entropic time series analysis"
  PLoS ONE 16(11):e0295448, 
  DOI:  10.1371/journal.pone.0259448
  https://www.EntropyHub.xyz

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
"""


  end