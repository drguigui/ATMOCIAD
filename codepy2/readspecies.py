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
import os
from readxmlimproved import *
def cleanstr(str):
	str=str.replace("%","\%")
	str=str.replace("->","$\\rightarrow$")
	str=str.replace("CO2++","CO$_2^{++}$")
	str=str.replace("CO2+","CO$_2^+$")
	str=str.replace("CO2","CO$_2$")
	str=str.replace("O2++","O$_2^{++}$")
	str=str.replace("O2+","O$_2^+$")
	str=str.replace("O2","O$_2$")
	str=str.replace("H4+","H$_4^+$")
	str=str.replace("H4","H$_4$")
	str=str.replace("H3+","H$_3^+$")
	str=str.replace("H3","H$_3$")
	str=str.replace("H2+","H$_2^+$")
	str=str.replace("H2","H$_2$")
	str=str.replace("H+","H$^+$")
	str=str.replace("N2++","N$_2^{++}$")
	str=str.replace("N2+","N$_2^{+}$")
	str=str.replace("N2","N$_2$")
	str=str.replace("N++","N$^{++}$")
	str=str.replace("N+","N$^{+}$")
	str=str.replace("CO(A1Pi)","CO($A^1\Pi$)")
	str=str.replace("CO(a3Pi)","CO($a^3\Pi$)")
	str=str.replace("O(1S)","O($^1$S)")
	str=str.replace("O(1D)","O($^1$D)")
	str=str.replace("O(3S)","O($^3$S)")
	str=str.replace("O(5S)","O($^5$S)")
	str=str.replace("O++","O$^{++}$")
	str=str.replace("O+(2P)","O$^+$($^2$P)")
	str=str.replace("O+(2D)","O$^+$($^2$D)")
	str=str.replace("O+(4S)","O$^+$($^4$S)")
	str=str.replace("O+(2P*)","O$^+$($^2$P$^*$)")
	str=str.replace("O+(4P*)","O$^+$($^4$P$^*$)")
	str=str.replace("O+","O$^+$")
	str=str.replace("C++","C$^{++}$")
	str=str.replace("C+","C$^+$")
	return str

def GetUnique(root,str):
	processlist=root.findall(str)
	if(len(processlist)==1):
		return  processlist[0].text
	print "Error, impossible to find the node "+str+" in your file"
	sys.exit()

def GetSpecie(spr):
	name=spr.attrib.get("name")
	state=spr.attrib.get("state")
	state=strip(state.replace("-NOTOT",""))
	return (name,state)

def fnametoref(stri):
	stri=strip(stri)
	stri=stri.replace("_","00")
	stri=stri.replace(".pdf","ppdf")
	return stri

def Readf(file,myspec):
	root=ET.parse(file).getroot()
	specie=GetUnique(root,".//Name").strip()
	collider=GetUnique(root,".//Collider").strip()
	
	if len(root.findall(".//RecommendedFile"))!=0:
		sp=root.findall(".//Specie")
		ref=""
		if len(root.findall(".//plotname"))!=0:
			ref="\\ref{"+fnametoref(root.find(".//plotname").text)+"}"
		
		if(len(sp)!=0):
		#	myspec={}
			for s in sp:
				thespecie=GetSpecie(s)

				if thespecie in myspec:
					if not (specie,collider,ref) in myspec[thespecie]:
						myspec[thespecie].append((specie,collider,ref))
				else:
					myspec[thespecie]=[(specie,collider,ref)]
	#		print myspec
	#		print specie," + ",collider, " -> ",len(myspec)," species created for all different process in",file



if __name__=="__main__":
	fichiers=glob("*/*/*xml")
	print "You have",len(fichiers),'Cross section files'

	myspec={}
	for f in fichiers:
		Readf(f,myspec)
	print myspec
	str=""

	specielist=[]
	for i in myspec.keys():
		(name,st)=i
		if not name in specielist:
			specielist.append(name)

#	print specielist
#	print sort(specielist)


	prkey=0
	for lis in sort(specielist):
		for i in myspec.keys():
			# We do something complicated to organise the outputs!
			(name,st)=i
			if name == lis:
				prkey+=1
				specie=name+"("+st+")"
		#	print "Sources for the creation of ",specie
				specie=cleanstr(specie)
				str+="\paragraph{Sources for the creation of "+specie+":} \n"
				str+="\\begin{itemize} \n"
				for j in myspec[i]:
					(sp,coll,ref)=j
		#		print "\t\t",sp," + ",coll," (Fig.",ref,")"
					str+= "\t\t\item "+cleanstr(sp)+" + "+coll.replace("ph","$\lambda$")+" (Fig."+ref+")\n"
				str+="\\end{itemize} \n"
	print "Number of keys in specielist :",len(myspec)," Number of printed keys :",prkey

	f=open("finalspecies.tex","w")
	f.writelines(str)
	f.close()


