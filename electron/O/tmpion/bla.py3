#!/usr/bin/env python3
from pylab import *

da = loadtxt("tmp")

plot(da[:,0], da[:,1] / da[:, 5], label="O+(4S)")
plot(da[:,0], da[:,2] / da[:, 5], label="O+(2D)")
plot(da[:,0], da[:,3] / da[:, 5], label="O+(2P)")
plot(da[:,0], da[:,4] / da[:, 5], label="O+(4P)")
legend()
show()
