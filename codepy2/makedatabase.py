#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
#from powerlaw import powerlaw
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
from StringIO import StringIO
from string import *
from glob import glob
import sys

def GetUnique(root,str):
	processlist=root.findall(str)
	if(len(processlist)==1):
		return  processlist[0].text
	print "Error, impossible to find the node "+str+" in your file"
	sys.exit()


def ReadNode2(node):
	mynode={}
	notes=GetUnique(node,"./Notes")
	if(type(notes)!=type("str")):
		notes=""
	mynode["Notes"]=notes
	nsource=node.find("Source")
	mynode["Source"]=nsource.text
	mynode["SourceType"]=nsource.attrib.get("type")

	uncertainty="100000%"
	mynode["EstimatedUncertainty"]=True
	
	if(len(node.findall("uncertainty"))>0):
		uncertainty=node.find("uncertainty").text+"%"
		mynode["EstimatedUncertainty"]=False
	else:
		if("uncertainty" in node.find("Cross").keys()):
			uncertainty=node.find("Cross").attrib.get("uncertainty")
			mynode["EstimatedUncertainty"]=False
	mynode["Uncertainty"]=uncertainty
	
	try:
		threshold=(node.attrib["threshold"])
	except:
		threshold="0"
	mynode["Threshold"]=threshold
	try:
		emin=(node.find("./Emin")).text
	except:
		emin=threshold
	
	try:
		emax=(node.find("./Emax")).text
	except:
		emax="-1"

	mynode["Recommended"]=False
	if(len(node.findall("./Recommended"))):
			mynode["Recommended"]=True

	if(len(node.findall("./EstimatedUncertainty"))):
			mynode["EstimatedUncertainty"]=True
			
	mynode["Extrapolated"]=False
	if(len(node.findall("./Extrapolated"))):
			mynode["Extrapolated"]=True
			
			
	mynode["Range"]=emin+":"+emax	



	return mynode




def ReadNode(node,arr):
	processus=GetUnique(node,"./Proc")
	#if(len(arr[processus]))==0:
	#	arr[processus]=[]
	if not processus in arr.keys():
		arr[processus]=[]
	arr[processus].append(ReadNode2(node))
	

def ReadFile(file,total,elastic,ionization,dissociation,excitation,emission,recommended):
	""" Reads the file 'file', and fill the sections considering the data inside
	Returns the name of the species and the name of the collider to check the correctness of the database
	"""
	root=ET.parse(file).getroot()

	specie=GetUnique(root,".//Name").strip()
	collider=GetUnique(root,".//Collider").strip()

	if(len(root.findall("./Final"))>0):
		return specie,collider

	recofile=False
	if(len(root.findall("./RecommendedFile"))>0):
		recofile=True
	
	processlist=root.findall(".//Process")
	print "Nous avons trouve ",len(processlist),"processus"
	for proc in processlist:
		if recofile:
			ReadNode(proc,recommended)
		else:
			section=GetUnique(proc,"./Section")
			if section=="ionization":
				ReadNode(proc,ionization)
			elif section=="dissociation":
				ReadNode(proc,dissociation)
			elif section=="excitation":
				ReadNode(proc,excitation)
			elif section=="emission":
				ReadNode(proc,emission)
			else:
				print("Impossible to find you process section: "+section+" in file "+file)
				sys.exit()


	processlist2=root.findall(".//ElasticCrs")
	for proc in processlist2:
		if recofile:
			recommended["Elastic"]=[ReadNode2(proc)]
		else:
			elastic.append(ReadNode2(proc))

	processlist3=root.findall(".//TotalCrs")
	for proc in processlist3:
		if recofile:
			recommended["Total"]=[ReadNode2(proc)]
		else:
			total.append(ReadNode2(proc))
	return specie,collider
	
	
def printnodedico(process,sp,coll):
	print("^ Process    ^   Reference   ^ Threshold        ^ Range of energy  ^   Uncertainty      ^    Notes    ^  Properties  ^")
	for i in process.keys():
		cnt=0
		procname = sp+" + "+coll+" -> "+i
		for j in process[i]:
			prc=""
			if cnt==0:
				prc=procname
				cnt=1
			props=""
			if(j["Recommended"]):
				props+="R"
			if(j["EstimatedUncertainty"]):
				props+="U"
				j["Uncertainty"]=colorisestr(j["Uncertainty"],"red")
			if(j["Extrapolated"]):
				props+="E"
			type=""
			if j["SourceType"]=="measurement":
				type="Meas"
			elif  j["SourceType"]=="bratio":
				type="Bran"
			elif  j["SourceType"]=="theory":
				type="Theo"
			elif  j["SourceType"]=="review":
				type="Revi"
			elif  j["SourceType"]=="adaptation":
				type="Adap"
			else:
				type=colorisestr("????","red")
				j["Source"]=colorisestr(j["Source"],"red")

				


			print("^ "+prc+" | "+type+" "+j["Source"]+" | "+j["Threshold"]+" | "+j["Range"]+" | "+j["Uncertainty"]+ " | "+j["Notes"]+" | "+props+" | ")
	print("")
		

def printnodesimple(process):
	print("^ Reference   ^ Threshold        ^ Range of energy  ^   Uncertainty      ^    Notes    ^  Properties  ^")
	for j in process:
		props=""
		if(j["Recommended"]):
			props+="R"
		if(j["EstimatedUncertainty"]):
			props+="U"
			j["Uncertainty"]=colorisestr(j["Uncertainty"],"red")
		if(j["Extrapolated"]):
			props+="E"
		type=""
		if j["SourceType"]=="measurement":
			type="Meas"
		elif  j["SourceType"]=="bratio":
			type="Bran"
		elif  j["SourceType"]=="theory":
			type="Theo"
		elif  (j["SourceType"]=="review"):
			type="Revi"
		elif  j["SourceType"]=="adaptation":
			type="Adap"
		else:
			type=colorisestr("????","red")
			j["Source"]=colorisestr(j["Source"],"red")
		#print j
		print("^ "+type+" "+j["Source"]+" | "+j["Threshold"]+" | "+j["Range"]+" | "+j["Uncertainty"]+ " | "+j["Notes"]+" | "+props+" | ")
	print("")
		

def colorisestr(str,color):
	m=str.replace(">","&gt;").replace(">","&lt;")
	return "<html><font color=\""+color+"\">"+m+"</font></html>"
	


def Printdatabase(specie,collider,total,elastic,ionization,dissociation,excitation,emission,recommended):
	print("====== Cross section of "+collider+" impact with "+specie+" ======")

	if len(total)!=0:
		print("===== Total Cross Section =====")
		printnodesimple(total)
	if len(elastic)!=0:
		print("===== Elastic Cross Section =====")
		printnodesimple(elastic)
	if len(ionization)!=0:
		print("===== Ionization Cross Section =====")
		printnodedico(ionization,specie,collider)
	if len(dissociation)!=0:
		print("===== Dissociation Cross Section =====")
		printnodedico(dissociation,specie,collider)
	if len(excitation)!=0:
		print("===== Excitation Cross Section =====")
		printnodedico(excitation,specie,collider)
	if len(emission)!=0:
		print("===== Emission Cross Section =====")
		printnodedico(emission,specie,collider)
	if len(recommended)!=0:
		print("===== Recommended File Cross Section =====")
		printnodedico(recommended,specie,collider)
def PrintLegend():
	print("==== Legend for the properties ====")
	print(" * R : Recommended cross section for the processus. It is used in the main file. The selection of the recommended cross section is based on the quality of the data (e.g. errorbars, comparison with other experiments), the possibility of extrapolation, and the origin of the work, coupled with the consistency (sum of recommended cross sections ~ Total cross section)")
	print(" * U : Estimated uncertainty: sometimes, the uncertainty is not given, because of theoretical work... The authors of the database have to estimate the uncertainty, but the quality of that estimation can be questionable. Moreover, when data from different sources have been adapted (e.g. for extrapolation), the uncertainty can be modified...")
	print(" * E : Validated for extrapolation: the extrapolation of these cross sections is plausible. For example, when an analytic function has been applied...")

if __name__=="__main__":
	files=glob("*xml")
	print("Your files:", files)
	# Definition of the different sections
	total=[]
	elastic=[]
	ionization={}
	dissociation={}
	excitation={}
	emission={}
	recommended={}
	specie=""
	collider=""
	for i in files:
		print("Processing file "+i)
		(newsp,newcoll)=ReadFile(i,total,elastic,ionization,dissociation,excitation,emission,recommended)
		#print i,newsp,newcoll
		if newsp=="":
			print "error in the file ",i," the specie name is badly set"
		if specie=="":
			specie=newsp
			collider=newcoll
		if specie!= newsp or collider !=newcoll:
			print "error : directory not homogenous. different species and/or colliders"
			print "Fatal error in file :",i
			sys.exit()

	Printdatabase(specie,collider,total,elastic,ionization,dissociation,excitation,emission,recommended)
	PrintLegend()


