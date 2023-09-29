# interval-methods
Numerical interval methods of root localization for nonlinear systems of equations<br />
Main usage of python library:<br />
<br />
import interval_arithmetics as ia<br />
import numpy as np<br />
n = 2<br />
m = 2<br />
X = ia.IntervalVector(np.full((n,), -10), np.full((n,), 9)) <br />#start guess of there root can be (in this case any coordinate between -10 and 9)<br />
<br />f = ia.Multinomial(n, ia.random_multinomial(n, m, seed=1234)) <br />#some function which you can call with numbers or intervals<br />
<br />X_new, has_root = ia.dissect_search_common(f, f.grad, X, steps=5000)<br />
<br />
#X_new - small interval box, in which, it is supposed to be root<br />
#has_root - bool, tell you if it is mathematically prooved that X_new has root inside.<br />
#if has_root == False, it means only that it isn't mathematically prooved but steel you can manually count:<br />
<br />

Magnitude = (ia.IntervalExtension(f))(X_new).M() <br />
has_root_numerically = (Magnitude < 1e-3)<br />
<br />
#magnitude (upper bound) for module of f on X_new, and if for example<br />
#numerically means that it is guarateed that f has root in box X_new with tolerance = 1e-3<br />
#in this case X_new is approximately [Interval(-1.004482, -1.004483], Interval[0.815660, 0.815661]]<br />
'''
