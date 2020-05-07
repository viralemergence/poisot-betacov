# poisot-betacov

Method 1 - k-NN (files named `Tanimoto`)

Method 2 - LF (files named `LinearFilter`)

This requires *Julia 1.3* to run - there is a `Project` and a `Manifest`
file already, so starting from a fresh *Julia* install all that is required
to get the correct packages is:

~~~ julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
~~~