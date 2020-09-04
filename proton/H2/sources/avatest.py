#!/usr/bin/env python
from __future__ import division
from pylab import *
d1 = loadtxt("ionizationAvakyan")
d2 = loadtxt("prodHpH2pAvakyan")

xscale('log')
yscale('log')



plot(d1[:,0], d1[:,1], label="e-")
plot(d2[:,0], d2[:,1], label="H2+")
plot(d2[:,0], d2[:,2], label="H+")

legend()
show()

