#!/usr/bin/env python
from pylab import *
import pylab
import sys


print sys.argv

def print_array(arr,colnumber=5,precision=6,number=14,type='e'):
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
		print prefix,colformat % tuple(arr[i*colnumber:(i+1)*colnumber])
	reste = size%colnumber
	if reste!=0:
		colfort='\t'.join([formatstr for i in range(reste)])
		print prefix,colfort%tuple(arr[-reste:])






if len(sys.argv)!=3:
	print "you should give a arguments : filename and position"
	sys.exit()

position=int(sys.argv[2])
if position==0:
	print "error : position should not be 0"
	sys.exit()



seff=pylab.loadtxt(sys.argv[1])



wavelength=seff[:,0]
data=seff[:,position]

wavelength=12398.42/(wavelength + 0.5)


wavelength=wavelength[::-1]
data=data[::-1]

#print "wavelength",wavelength
#print "data",data


print "\t\t<Process name=\"\" electrons=\"\" threshold=\"\">"
print "\t\t\t<Species>"
print "\t\t\t\t<Specie name=\"\" state=\"\"/>"
print "\t\t\t</Species>"
print "\t\t\t<Egrid unit=\"eV\">"
print_array(wavelength,type='f')
print "\t\t\t</Egrid>"
print "\t\t\t<Cross unit=\"cm2\">"
print_array(data)
print "\t\t\t</Cross>"
print "\t\t</Process>"




