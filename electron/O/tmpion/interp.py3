#!/usr/bin/env python3
from pylab import *
from scipy.interpolate import interp1d
da = loadtxt("tmp")

f1 = interp1d(da[:,0], da[:,1] / da[:, 5])
f2 = interp1d(da[:,0], da[:,2] / da[:, 5])
f3 = interp1d(da[:,0], da[:,3] / da[:, 5])
f4 = interp1d(da[:,0], da[:,4] / da[:, 5])

db = loadtxt("tmp2")
yscale("log")
plot(db[:,0], db[:,1], label="Total")
plot(db[:,0], db[:,1] * f1(db[:,0]), label="O+(4S)")
plot(db[:,0], db[:,1] * f2(db[:,0]), label="O+(2D)")
plot(db[:,0], db[:,1] * f3(db[:,0]), label="O+(2P)")
plot(db[:,0], db[:,1] * f4(db[:,0]), label="O+(4P)")



def print_array(arr,colnumber=5,precision=4,number=14,type='e'):
	""" Print the array in the formatted way:
	arr : the array
	colnumber : the number of printed columns
	precision : 
	number
	type: e, f, g, i
	"""
	prefix='\t\t\t'
	islogger=False
	formatstr= "%i.%i"%(number,precision)
	formatstr='%'+formatstr+type
	size=len(arr)
	tbl=[ formatstr for i in range(colnumber)]
	colformat='\t'.join(tbl)
	for i in range(int(size/colnumber)):
		print(prefix,colformat % tuple(arr[i*colnumber:(i+1)*colnumber]))
	reste = size%colnumber
	if reste!=0:
		colfort='\t'.join([formatstr for i in range(reste)])
		print(prefix,colfort%tuple(arr[-reste:]))

print_array( db[:,1] * f1(db[:,0]))
print_array( db[:,1] * f2(db[:,0]))
print_array( db[:,1] * f3(db[:,0]))
print_array( db[:,1] * f4(db[:,0]))












print 

#plot(da[:,0], da[:,1] / da[:, 5], label="O+(4S)")
#plot(da[:,0], da[:,2] / da[:, 5], label="O+(2D)")
#plot(da[:,0], da[:,3] / da[:, 5], label="O+(2P)")
#plot(da[:,0], da[:,4] / da[:, 5], label="O+(4P)")
legend()
show()
