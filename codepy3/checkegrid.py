#!/usr/bin/env python
# -*- coding:Utf-8 -*-

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
				print("PROBLEME", leg)
				print(dataenergy[i], dataenergy[i+1])
				return
	else:
		for i in range(len(dataenergy)-1):
			if not(dataenergy[i] < dataenergy[i+1]):
				print("PROBLEME", leg)
				print(dataenergy[i], dataenergy[i+1])
				return
	print("no pb")


if __name__=="__main__":
	if(len(sys.argv)<2):
		print("veuillez donner un nom de fichier")
		sys.exit()
	filename=sys.argv[1]
	print("Votre nom de fichier : ",filename)
	
	xscale("log")
	yscale("log")
	
	root=ET.parse(filename).getroot()
	processlist=root.findall(".//Process")
	print("Nous avons trouve ",len(processlist),"processus")
	
	for proc in processlist:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc)
		else:
			print("no pb")
	#	print proc.attrib["name"]
	#	print len(proc.findall("Ionization"))
	#	print proc.find("Cross").text

	processlist2=root.findall(".//ElasticCrs")
	for proc in processlist2:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc)
		else:
			print("no pb")
	
	processlist3=root.findall(".//TotalCrs")
	for proc in processlist3:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc)
		else:
			print("no pb")
	
	print("fini")
	


