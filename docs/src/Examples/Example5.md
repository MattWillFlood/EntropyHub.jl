# Example 5: Multiscale Entropy Object - MSobject()

!!! warning "Note:"

    The base and cross- entropy functions used in the multiscale entropy calculation 
    are declared by passing EntropyHub functions to MSobject(), not string names.

Create a multiscale entropy object (`Mobj`) for multiscale fuzzy entropy, calculated with
an embedding dimension of 5, a time delay (`tau`) of 2, using a sigmoidal fuzzy function with
the `r` scaling parameters (3, 1.2).

```@example
using EntropyHub # hide
Mobj = MSobject(FuzzEn, m = 5, tau = 2, Fx = "sigmoid", r = (3, 1.2))
```

Create a multiscale entropy object (`Mobj`) for multiscale corrected-cross-conditional entropy, 
calculated with an embedding dimension of 6 and using a 11-symbolic data transform.

```@example
using EntropyHub # hide
Mobj = MSobject(XCondEn, m = 6, c = 11)
```

