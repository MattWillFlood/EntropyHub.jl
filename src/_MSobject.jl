module _MSobject

export MSobject
    """
        Mobj = MSobject() 

    Returns a multiscale entropy object (`Mobj`) based on that originally 
    proposed by Costa et. al. (2002) using the following default parameters:
    EnType = SampEn(), embedding dimension = 2, time delay = 1, 
    radius = 0.2*SD(`Sig`), logarithm = natural
    
        Mobj = MSobject(EnType::Function)

    Returns a multiscale entropy object by passing the entropy function
    (`EnType`) and the specifying default parameters for that entropy function.
    To see the default parameters for a particular entropy method,          
    type:  `? EntropyHub.EnType`  \n
    (e.g.  `? EntropyHub.SampEn`)

        Mobj = MSobject(EnType::Function; kwargs...)

    Returns a multiscale entropy object using the specified entropy method
    (`EnType`) and the 'keyword' parameters for that particular method.
    To see the default parameters for a particular entropy method,          
    type:   ? EntropyHub.EnType   (e.g.  ? EntropyHub.SampEn)

    EnType can be any of the following entropy functions:

    # Base Entropies:
        -----------------
        `ApEn`      - Approximate Entropy  \n
        `SampEn`    - Sample Entropy   \n
        `FuzzEn`    - Fuzzy Entropy    \n
        `K2En`      - Kolmogorov Entropy   \n
        `PermEn`    - Permutation Entropy	 \n   
        `CondEn`    - Conditional Entropy	   \n 
        `DistEn`    - Distribution Entropy	    \n
        `DispEn`    - Dispersion Entropy	    \n
        `SpecEn`    - Spectral Entropy   \n
        `SyDyEn`    - Symbolic Dynamic Entropy	  \n  
        `IncrEn`    - Increment Entropy	    \n
        `CoSiEn`    - Cosine Similarity Entropy	\n    
        `PhasEn`    - Phase Entropy	    \n
        `SlopEn`    - Slope Entropy       \n 
        `BubbEn`    - Bubble Entropy   \n
        `GridEn`    - Gridded Distribution Entropy  \n  	
        `EnofEn`    - Entropy of Entropy	\n    
        `AttnEn`    - Attention Entropy    \n
        `DivEn`     - Diversity Entropy \n
        `RangEn`    - Range Entropy \n

        
    # Cross Entropies:
        ------------------
        `XApEn`     - Cross-Approximate Entropy    \n 
        `XSampEn`   - Cross-Sample Entropy  \n
        `XFuzzEn`   - Cross-Fuzzy Entropy   \n
        `XK2En`     - Cross-Kolmogorov Entropy  \n
        `XPermEn`   - Cross-Permutation Entropy   \n  
        `XCondEn`   - Cross-Conditional Entropy    \n
        `XDistEn`   - Cross-Distribution Entropy    \n
        `XSpecEn`   - Cross-Spectral Entropy   \n


    # Multivariate Entropies:
        ------------------
        `MvSampEn`   - Multivariate Sample Entropy  \n
        `MvFuzzEn`   - Multivariate Fuzzy Entropy   \n
        `MvPermEn`   - Multivariate Permutation Entropy   \n  
        `MvCoSiEn`   - Multivariate Cosine Similarity Entropy    \n
        `MvDispEn`   - Multivariate Dispersion Entropy    \n

    # See also `MSEn`, `XMSEn`, `MvMSEn`, `rMSEn`, `cMSEn`, `hMSEn`, `rXMSEn`, `cXMSEn`, `hXMSEn`, `cMvMSEn`

    """
    function MSobject(EnType::Function=SampEn; kwargs...)

    Chk = ["ApEn";"SampEn";"FuzzEn";"K2En";"PermEn";"CondEn";"DistEn";"DivEn";
        "DispEn";"SyDyEn";"IncrEn";"CoSiEn";"PhasEn";"SpecEn";"SlopEn";"RangEn";
        "GridEn";"BubbEn";"EnofEn";"AttnEn";
        "XApEn";"XSampEn";"XFuzzEn";"XPermEn";
        "XCondEn";"XDistEn";"XSpecEn";"XK2En";
        "MvSampEn";"MvFuzzEn";"MvPermEn";
        "MvCoSiEn";"MvDispEn"]
    (String(Symbol(EnType)) in Chk) ? nothing :
        error("EnType:      must be a valid entropy function name.
        For more info, type:   julia>  ? EntropyHub.MSobject")

    #Mobj = (Func=getfield(Main.EntropyHub, Symbol(EnType)), kwargs...) 

    Mobj = (Func= EnType, kwargs...) 

    return Mobj
    end

end

"""
Copyright 2024 Matthew W. Flood, EntropyHub

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