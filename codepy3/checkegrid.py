#!/usr/bin/env python
# -*- coding:Utf-8 -*-


""" Checks to see if cross-section data exist for all processes in the xml file

Examples
----------

To use function, type in the following on the command line:

    python checkegrid.py filename

Arguments
----------
filename : str
    The file location of the xml file that contains the cross-sections

Returns
-------

None

Notes
--------

If there is no issue loading the cross-sections from the xml file, then "No problem" will be printed. If there is a problem, then "Problem" will be printed, along with the name of the species. 
    
Notes created by B Hegyi, 10/19/20

"""

from powerlaw import powerlaw
try:
    import xml.etree.ElementTree as ET # in python >=2.5
except ImportError:
    try:
        import cElementTree as ET # effbot's C module
    except ImportError:
        try:
            import elementtree.ElementTree as ET # effbot's pure Python module
        except ImportError:
            try:
                import lxml.etree as ET # ElementTree API using libxml2
            except ImportError:
                import warnings
                warnings.warn("could not import ElementTree "
                              "(http://effbot.org/zone/element-index.htm)")
                # Or you might just want to raise an ImportError here.
from pylab import *
from io import StringIO
from string import *

def PlotStdNode(vNode):
	leg=""
	if(0==len(vNode.findall("legend"))):
		try:
			leg=vNode.attrib["name"]
		except:
			leg="Elastic"
	else:
		leg=vNode.find("legend").text
	
	dataenergy=loadtxt(StringIO(vNode.find("Egrid").text.replace("\n"," ")))
	
	
	if dataenergy[0] > dataenergy[-1]:
		for i in range(len(dataenergy)-1):
			if not(dataenergy[i] > dataenergy[i+1]):
				print("Problem:", leg)
				print(dataenergy[i], dataenergy[i+1])
				return
	else:
		for i in range(len(dataenergy)-1):
			if not(dataenergy[i] < dataenergy[i+1]):
				print("Problem:", leg)
				print(dataenergy[i], dataenergy[i+1])
				return
	print("No issues!")


if __name__=="__main__":
	if(len(sys.argv)<2):
		print("Please include a file name on the command line")
		sys.exit()
	filename=sys.argv[1]
	print("File name : ",filename)
	
	xscale("log")
	yscale("log")
	
	root=ET.parse(filename).getroot()
	processlist=root.findall(".//Process")
	print("We have found",len(processlist),"processes")
	
	for proc in processlist:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc)
		else:
			print("No issues!")
	#	print proc.attrib["name"]
	#	print len(proc.findall("Ionization"))
	#	print proc.find("Cross").text

	processlist2=root.findall(".//ElasticCrs")
	for proc in processlist2:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc)
		else:
			print("No issues!")
	
	processlist3=root.findall(".//TotalCrs")
	for proc in processlist3:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc)
		else:
			print("No issues!")
	
	print("Finished")
	


