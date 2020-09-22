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
from glob import glob
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
				return False, leg,dataenergy[i], dataenergy[i+1]
	else:
		for i in range(len(dataenergy)-1):
			if not(dataenergy[i] < dataenergy[i+1]):
				return False, leg,dataenergy[i], dataenergy[i+1]
	return True,"","",""


if __name__=="__main__":
	fichiers=glob("*/*/*xml")
	print("You have",len(fichiers),'Cross section files')
	
	for filename in fichiers:
		try:
			root=ET.parse(filename).getroot()
		except:
			print("OOOPS, bad XML file:",filename)
			continue
		processlist=root.findall(".//Process")
		for proc in processlist:
			if(0==len(proc.findall("Shirai"))):
				rez, leg, E1, E2 = PlotStdNode(proc)
			if(not rez):
				print(filename,"\t =>\t", leg, E1, E2)
		
		processlist2=root.findall(".//ElasticCrs")
		for proc in processlist2:
			if(0==len(proc.findall("Shirai"))):
				rez, leg, E1, E2 = PlotStdNode(proc)
			if(not rez):
				print(filename,"\t =>\t", leg, E1, E2)
		processlist3=root.findall(".//TotalCrs")
		for proc in processlist3:
			if(0==len(proc.findall("Shirai"))):
				rez, leg, E1, E2 = PlotStdNode(proc)
			if(not rez):
				print(filename,"\t =>\t", leg, E1, E2)
		


