# interval-methods
Numerical interval methods of root localization for nonlinear systems of equations
Main usage of python library:

import interval_arithmetics as ia
import numpy as np
n = 2
m = 2
X = ia.IntervalVector(np.full((n,), -10), np.full((n,), 9)) #start guess of there root can be (in this case any coordinate between -10 and 9)
f = ia.Multinomial(n, ia.random_multinomial(n, m, seed=1234)) #some function which you can call with numbers or intervals
X_new, has_root = ia.dissect_search_common(f, f.grad, X, steps=5000)

#X_new - small interval box, in which, it is supposed to be root
#has_root - bool, tell you if it is mathematically prooved that X_new has root inside.
#if has_root == False, it means only that it isn't mathematically prooved but steel you can manually count:
Magnitude = (ia.IntervalExtension(f))(X_new).M() # magnitude (upper bound) for module of f on X_new, and if for example
if Magnitude < 1e-3:
    has_root_numerically = True
else:
    has_root_numerically = False
#numerically means that it is guarateed that f has root in box X_new with tolerance = 1e-3
#in this case X_new is approximately [Interval(-1.004482, -1.004483], Interval[0.815660, 0.815661]]
