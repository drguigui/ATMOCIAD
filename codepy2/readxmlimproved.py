#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
from powerlaw import powerlaw
from MathFun import intloglog
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
from itertools import cycle
colo=['b','g','r','c','m','y','k']
das=[(None,None),(4,1,1,1),(2,1), (5,1,2,1),(4,1),(6,1), (10,1)]
lsty=[  dict(zip(('color','dashes'),(col,sty))) for sty in das for col in colo  ]
linestyleskwargs = cycle(lsty)


def enetoang(ene):
	""" Convert the energy in eV to a wavelength in Angstrom"""
	return 12398.42/ene

class NewShiraiCH4:
	""" Class to read the cross sections defined in the CH4 paper of Shirai et al"""
	def __init__(self,emin,emax,threshold,eqtype,dataeq):
		self.emin=emin
		self.emax=emax
		self.thresh=threshold/1000. # the threshold is put in KeV units for the computation
		self.eq=eqtype
		self.avals=dataeq
		self.sigma=1.E-16 # cm2
		self.Er=1.361E-2 # keV
	def f1(self,E,c1,c2):
		return self.sigma*c1*(E/self.Er)**c2
	def f2(self,E,c1,c2,c3,c4):
		return self.f1(E,c1,c2)/(1+(E/c3)**(c2+c4))
	def f3(self,E,c1,c2,c3,c4,c5,c6):
		return self.f1(E,c1,c2)/(1+(E/c3)**(c2+c4)+(E/c5)**(c2+c6))

	def S1(self,Ei):
		E=Ei-self.thresh
		return self.f1(E,self.avals[0],self.avals[1])
	def S2(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])
	def S3(self,Ei):
		E=Ei-self.thresh
		return self.f1(E,self.avals[0],self.avals[1])+ self.f2(E,self.avals[2],self.avals[3],self.avals[4],self.avals[5])
	def S4(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.avals[4]*self.f2(E/self.avals[5],self.avals[0],self.avals[1],self.avals[2],self.avals[3])

	def S5(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f2(E,self.avals[4],self.avals[5],self.avals[6],self.avals[3])

	def S6(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f2(E,self.avals[4],self.avals[5],self.avals[6],self.avals[7])

	def S7(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f2(E,self.avals[4],self.avals[5],self.avals[6],self.avals[3])+self.f2(E,self.avals[7],self.avals[8],self.avals[9],self.avals[3])

	def S8(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f2(E,self.avals[4],self.avals[5],self.avals[6],self.avals[7])+self.f2(E,self.avals[8],self.avals[9],self.avals[10],self.avals[11])

	def S9(self,Ei):
		E=Ei-self.thresh
		return self.f3(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3],self.avals[4],self.avals[5])


	def S10(self,Ei):
		E=Ei-self.thresh
		return self.f1(E,self.avals[0],self.avals[1])+ self.f3(E,self.avals[2],self.avals[3],self.avals[4],self.avals[5],self.avals[6],self.avals[7])

	def S11(self,Ei):
		E=Ei-self.thresh
 		return self.f3(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3],self.avals[4],self.avals[5])+self.avals[6]*self.f3(E/self.avals[7],self.avals[0],self.avals[1],self.avals[2],self.avals[3],self.avals[4],self.avals[5])

	def S12(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f3(E,self.avals[4],self.avals[5],self.avals[6],self.avals[7],self.avals[8],self.avals[3])

	def S13(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f3(E,self.avals[4],self.avals[5],self.avals[6],self.avals[7],self.avals[8],self.avals[9])+self.f2(E,self.avals[10],self.avals[11],self.avals[12],self.avals[9])

	def S14(self,E):
		return self.sigma*self.avals[0]*( log(E/self.thresh)+self.avals[1])/(self.thresh*E*(1+(self.avals[2]/(E-self.thresh))**self.avals[3]))

	
	def ReturnCrs(self,E):
		result={1: lambda x: self.S1(x),
			2: lambda x: self.S2(x),
			3: lambda x: self.S3(x),
			4: lambda x: self.S4(x),
			5: lambda x: self.S5(x),
			6: lambda x: self.S6(x),
			7: lambda x: self.S7(x),
			8: lambda x: self.S8(x),
			9: lambda x: self.S9(x),
			10: lambda x: self.S10(x),
			11: lambda x: self.S11(x),
			12: lambda x: self.S12(x),
			13: lambda x: self.S13(x),
			14: lambda x: self.S14(x)}[self.eq](E)
		return result
class NewShiraiCO2:
	""" Class to read the cross sections defined in the CO2 paper of Shirai et al"""
	def __init__(self,emin,emax,threshold,eqtype,dataeq):
		self.emin=emin
		self.emax=emax
		self.thresh=threshold/1000.
		self.eq=eqtype
		self.avals=dataeq
		self.sigma=1.E-16 # cm2
		self.Er=1.361E-2 # keV
	def f1(self,E,c1,c2):
		return self.sigma*c1*(E/self.Er)**c2
	def f2(self,E,c1,c2,c3,c4):
		return self.f1(E,c1,c2)/(1+(E/c3)**(c2+c4))
	def f3(self,E,c1,c2,c3,c4,c5,c6):
		return self.f1(E,c1,c2)/(1+(E/c3)**(c2+c4)+(E/c5)**(c2+c6))

	def S1(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])
	def S2(self,Ei):
		E=Ei-self.thresh
		return self.f1(E,self.avals[0],self.avals[1])+ self.f2(E,self.avals[2],self.avals[3],self.avals[4],self.avals[5])
	
	def S3(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f2(E,self.avals[4],self.avals[5],self.avals[6],self.avals[7])
	
	def S4(self,Ei):
		E=Ei-self.thresh
		return self.f1(E,self.avals[0],self.avals[1])+ self.f2(E,self.avals[2],self.avals[3],self.avals[4],self.avals[5])+ self.f2(E,self.avals[6],self.avals[7],self.avals[8],self.avals[9])
	
	def S5(self,Ei):
		E=Ei-self.thresh
		return self.f3(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3],self.avals[4],self.avals[5])

	def S6(self,Ei):
		E=Ei-self.thresh
		return self.f1(E,self.avals[0],self.avals[1])+ self.f3(E,self.avals[2],self.avals[3],self.avals[4],self.avals[5],self.avals[6],self.avals[7])

	def S7(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f3(E,self.avals[4],self.avals[5],self.avals[6],self.avals[7],self.avals[8],self.avals[3])

	def S8(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f3(E,self.avals[4],self.avals[5],self.avals[6],self.avals[7],self.avals[8],self.avals[9])

	def S9(self,Ei):
		E=Ei-self.thresh
		return self.f1(E,self.avals[0],self.avals[1])+ self.f2(E,self.avals[2],self.avals[3],self.avals[4],self.avals[5])+self.f3(E,self.avals[6],self.avals[7],self.avals[8],self.avals[9],self.avals[10],self.avals[11])

	def S10(self,E):
		return self.sigma*self.avals[0]*( log(E/self.thresh)+self.avals[1])/(self.thresh*E*(1+(self.avals[2]/(E-self.thresh))**self.avals[3]))

	
	def ReturnCrs(self,E):
		result={1: lambda x: self.S1(x),
			2: lambda x: self.S2(x),
			3: lambda x: self.S3(x),
			4: lambda x: self.S4(x),
			5: lambda x: self.S5(x),
			6: lambda x: self.S6(x),
			7: lambda x: self.S7(x),
			8: lambda x: self.S8(x),
			9: lambda x: self.S9(x),
			10: lambda x: self.S10(x)}[self.eq](E)
		return result



class NewShiraiN2:
	""" Class to read the cross sections defined in the N2 paper of Tabata, Shirai et al"""
	def __init__(self,emin,emax,threshold,eqtype,dataeq):
		self.emin=emin
		self.emax=emax
		self.thresh=threshold/1000.
		self.eq=eqtype
		self.avals=dataeq
		self.sigma=1.E-16 # cm2
		self.Er=1.361E-2 # keV
	def f1(self,E,c1,c2):
		return self.sigma*c1*(E/self.Er)**c2
	def f2(self,E,c1,c2,c3,c4):
		return self.f1(E,c1,c2)/(1+(E/c3)**(c2+c4))
	def f3(self,E,c1,c2,c3,c4,c5,c6):
		return self.f1(E,c1,c2)/(1+(E/c3)**(c2+c4)+(E/c5)**(c2+c6))

	def S1(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])
	
	def S2(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f2(E,self.avals[4],self.avals[5],self.avals[6],self.avals[3])
	def S3(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f2(E,self.avals[4],self.avals[5],self.avals[6],self.avals[7])
	def S4(self,Ei):
		E=Ei-self.thresh
		return self.f2(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3])+self.f2(E,self.avals[4],self.avals[5],self.avals[6],self.avals[7])+self.f2(E,self.avals[8],self.avals[9],self.avals[10],self.avals[11])

	def S5(self,Ei):
		E=Ei-self.thresh
		return self.f3(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3],self.avals[4],self.avals[5])


	def S6(self,Ei):
		E=Ei-self.thresh
		return  self.f3(E,self.avals[0],self.avals[1],self.avals[2],self.avals[3],self.avals[4],self.avals[5])+self.f2(E,self.avals[6],self.avals[7],self.avals[8],self.avals[9])


	def S7(self,E):
		return self.sigma*self.avals[0]*( log(E/self.thresh)+self.avals[1])/(self.thresh*E*(1+(self.avals[2]/(E-self.thresh))**self.avals[3]))

	
	def ReturnCrs(self,E):
		result={1: lambda x: self.S1(x),
			2: lambda x: self.S2(x),
			3: lambda x: self.S3(x),
			4: lambda x: self.S4(x),
			5: lambda x: self.S5(x),
			6: lambda x: self.S6(x),
			7: lambda x: self.S7(x),
			}[self.eq](E)
		return result






def PlotShiraiNode(vNode,bRfile,vIsLambda,vIsEnergyTbl=False,vEtbl=[]):
	""" Function to plot a node as defined in Shirai et al"""
	leg=""
	if(0==len(vNode.findall("legend"))):
		leg=vNode.attrib["name"]
	else:
		leg=vNode.find("legend").text
	
	if bRfile:
		if(0!=len(vNode.findall("Proc"))):
			leg=vNode.find("Proc").text

	try:
		threshold=float(vNode.attrib["threshold"])
	except:
		threshold=0
	uncertainty=float(vNode.find("uncertainty").text)
	Emin=float(vNode.find("Emin").text)
	Emax=float(vNode.find("Emax").text)

	tid=vNode.find("Equation").attrib["article_id"]
	tip=int(vNode.find("Equation").attrib["type"])
	params=loadtxt(StringIO(vNode.find("params").text.replace("\n"," ")))

	if(tid=="CH4"):
		shirai=NewShiraiCH4(Emin,Emax,threshold,tip,params)
	if(tid=="CO2"):
		shirai=NewShiraiCO2(Emin,Emax,threshold,tip,params)
	if(tid=="N2"):
		shirai=NewShiraiN2(Emin,Emax,threshold,tip,params)

	ene=arange(Emin,Emax,(Emax-Emin)/100.)
	ene=array(powerlaw(Emin,Emax,100))
	# For the extrapolation
	if vIsEnergyTbl:
		ene=array(vEtbl)
	cross=shirai.ReturnCrs(ene*1E-3)

	datauncert=cross*uncertainty/100.
	if(vIsLambda):
		errorbar(enetoang(ene),cross,yerr=datauncert,label=leg)
	else:
		errorbar(ene,cross,yerr=datauncert,label=leg)


def PlotStdNode(vNode,bRfile,vIsLambda,vIsEnergyTbl=False,vEtbl=[]):
	""" Function to plot a standard node"""
	leg=""
	if(0==len(vNode.findall("legend"))):
		try:
			leg=vNode.attrib["name"]
		except:
			leg="Elastic"
	else:
		leg=vNode.find("legend").text
	
	if bRfile:
		if(0!=len(vNode.findall("Proc"))):
			leg=vNode.find("Proc").text
	fact=1
	if("fact" in vNode.find("Egrid").keys()):
		fact=float(vNode.find("Egrid").attrib.get("fact"))
#	print "Votre facteur :",fact
#	print loadtxt(StringIO((vNode.find("Egrid").text).replace("\n"," ")))
	dataenergy=loadtxt(StringIO(vNode.find("Egrid").text.replace("\n"," ")))*fact
	fact=1
	if("fact" in vNode.find("Cross").keys()):
		fact=float(vNode.find("Cross").attrib.get("fact"))
	print "Votre facteur :",fact
	datacrs=loadtxt(StringIO(vNode.find("Cross").text.replace("\n"," ")))*fact
#	print datacrs
	
	uncertainty=0
	datauncert=zeros((len(datacrs)))
	if("uncertainty" in vNode.find("Cross").keys()):
		uncertainty=vNode.find("Cross").attrib.get("uncertainty")
		if (uncertainty.find("%")):
			uncert=float(uncertainty.replace("%",""))/100.
			print "Uncertainty factor",uncert
			datauncert=datacrs*uncert
		else:
			value=float(uncertainty)*fact
			print "Uncertainty value"
			datauncert=ones((len(datacrs)))*value
	if(vIsEnergyTbl):
		vEtbl=array(vEtbl)
	if(vIsLambda):
		if(vIsEnergyTbl):
			style = linestyleskwargs.next()
			print style
			testcrs=intloglog(dataenergy,datacrs,vEtbl)
			if("threshold" in vNode.keys()):
				threshold=float(vNode.attrib.get("threshold"))
				for i in range(len(vEtbl)):
					if(vEtbl[i]<threshold):
						testcrs[i]=0
			plot(enetoang(vEtbl),testcrs,**style)
			errorbar(enetoang(dataenergy),datacrs,yerr=datauncert,label=leg,**style)
		else:
			errorbar(enetoang(dataenergy),datacrs,yerr=datauncert,label=leg)

	else:
		if(vIsEnergyTbl):
			style = linestyleskwargs.next()
			print style
			testcrs=intloglog(dataenergy,datacrs,vEtbl)
			if("threshold" in vNode.keys()):
				threshold=float(vNode.attrib.get("threshold"))
				for i in range(len(vEtbl)):
					if(vEtbl[i]<threshold):
						testcrs[i]=0
			plot((vEtbl),testcrs,**style)
			errorbar(dataenergy,datacrs,yerr=datauncert,label=leg,**style)
		else:
			errorbar(dataenergy,datacrs,yerr=datauncert,label=leg)
	




def CheckRoot(root):
	""" Checks if the files contains the markup for automatic plotting
	Extracts the information of these markups when present
	"""
	if len(root.findall(".//title"))==0 :
		print "Please provide the title /title"
		sys.exit()
	titre=root.find(".//title").text

	if (root.findall(".//Emin")==[]) or (root.findall(".//Emax")==[]) or (root.findall(".//Cmin")==[]) or (root.findall(".//Cmax")==[]):
		print "Please provide the extreme value /Emin /Emax /Cmin /Cmax"
		sys.exit()
	emin=float(root.find(".//Emin").text)
	emax=float(root.find(".//Emax").text)
	cmin=float(root.find(".//Cmin").text)
	cmax=float(root.find(".//Cmax").text)

	if root.findall(".//plotname")==[]:
		print "Please provide the name of the plot /plotname"
		sys.exit()
	figname=root.find(".//plotname").text
	lfig="" # Plot lambda
	efig="" # Plot extrapolated
	elfig="" # Plot lambda extrapolated

	figmult=1

	if len(root.findall(".//lambplotname"))!=0:
		lfig=root.find(".//lambplotname").text

	if (len(root.findall(".//exlambplotname"))!=0):
		elfig=root.find(".//exlambplotname").text

	if len(root.findall(".//explotname"))!=0:
		efig=root.find(".//explotname").text

	if len(root.findall(".//figsize"))!=0:
		figmult=float(root.find(".//figsize").text)
	
	
	return titre,emin,emax,cmin,cmax,figname,lfig,efig,elfig,figmult


def PlotFile(filename):
	""" Performs the plot operations for the given file"""
	root=ET.parse(filename).getroot()
	
	titre,emin,emax,cmin,cmax,figname,lfig,efig,elfig,figmult=CheckRoot(root)
	
	if figmult!=1:
		print "multiplication de la taille de votre figure"
		print figmult
		figure(figsize=[8*figmult,6*figmult])
	clf()
	xscale("log")
	yscale("log")
	processlist=root.findall(".//Process")
	processlist2=root.findall(".//ElasticCrs")
	processlist3=root.findall(".//TotalCrs")

	bRfile=False
	if len(root.findall(".//RecommendedFile"))!=0:
		bRfile=True
		
	print "We found ",len(processlist),"processes"
	
	for proc in processlist:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc,bRfile,False,False)
		else:
			PlotShiraiNode(proc,bRfile,False,False)
	for proc in processlist2:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc,bRfile,False,False)
		else:
			PlotShiraiNode(proc,bRfile,False,False)
	for proc in processlist3:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc,bRfile,False,False)
		else:
			PlotShiraiNode(proc,bRfile,False)
	
	title(titre)
	xlabel("Energy (eV)")
	ylabel(r"Cross section (cm$^2$)")
	axis([emin,emax,cmin,cmax])
	if(abs(2*(emax-emin)/(emax+emin))<1):
		xticks(arange(emin,emax+(emax-emin)/4,(emax-emin)/4),arange(emin,emax+(emax-emin)/4,(emax-emin)/4))
	legend(loc="best")
	savefig(figname)


	if lfig!="":
		clf()
		linestyleskwargs = cycle(lsty)
		xscale("linear")
		yscale("log")
		print "creation of the lambda figure"
		for proc in processlist:
			if(0==len(proc.findall("Shirai"))):
				PlotStdNode(proc,bRfile,True)
			else:
				PlotShiraiNode(proc,bRfile,True)
		for proc in processlist2:
			if(0==len(proc.findall("Shirai"))):
				PlotStdNode(proc,bRfile,True)
			else:
				PlotShiraiNode(proc,bRfile,True)
		for proc in processlist3:
			if(0==len(proc.findall("Shirai"))):
				PlotStdNode(proc,bRfile,True)
			else:
				PlotShiraiNode(proc,bRfile,True)
	
		title(titre)
		xlabel(r"Wavelength ($\AA$)")
		ylabel(r"Cross section (cm$^2$)")
		axis([enetoang(emax),enetoang(emin),cmin,cmax])
		legend(loc="best")
		savefig(lfig)

	Etbl=powerlaw(0.1,100000,500)
	if efig!="":
		global linestyleskwargs
		linestyleskwargs = cycle(lsty)
		print "creation of the extrapolated figure"
		clf()
		xscale("log")
		yscale("log")
		print "creation of the lambda figure"
		for proc in processlist:
			if(0==len(proc.findall("Shirai"))):
				PlotStdNode(proc,bRfile,False,True,Etbl)
			else:
				PlotShiraiNode(proc,bRfile,False,True,Etbl)
		for proc in processlist2:
			if(0==len(proc.findall("Shirai"))):
				PlotStdNode(proc,bRfile,False,True,Etbl)
			else:
				PlotShiraiNode(proc,bRfile,False,True,Etbl)
		for proc in processlist3:
			if(0==len(proc.findall("Shirai"))):
				PlotStdNode(proc,bRfile,False,True,Etbl)
			else:
				PlotShiraiNode(proc,bRfile,False,True,Etbl)
	
		title(titre)
		xlabel("Energy (eV)")
		ylabel(r"Cross section (cm$^2$)")
		axis([emin,emax,cmin,cmax])
		if(abs(2*(emax-emin)/(emax+emin))<1):
			xticks(arange(emin,emax+(emax-emin)/4,(emax-emin)/4),arange(emin,emax+(emax-emin)/4,(emax-emin)/4))
		legend(loc="best")
		savefig(efig)
	if elfig!="":
		print "creation of the extrapolated lambda figure"
		global linestyleskwargs
		linestyleskwargs = cycle(lsty)
		clf()
		xscale("linear")
		yscale("log")
		print "creation of the lambda figure"
		for proc in processlist:
			if(0==len(proc.findall("Shirai"))):
				PlotStdNode(proc,bRfile,True,True,Etbl)
			else:
				PlotShiraiNode(proc,bRfile,True,True,Etbl)
		for proc in processlist2:
			if(0==len(proc.findall("Shirai"))):
				PlotStdNode(proc,bRfile,True,True,Etbl)
			else:
				PlotShiraiNode(proc,bRfile,True,True,Etbl)
		for proc in processlist3:
			if(0==len(proc.findall("Shirai"))):
				PlotStdNode(proc,bRfile,True,True,Etbl)
			else:
				PlotShiraiNode(proc,bRfile,True,True,Etbl)
	
		title(titre)
		xlabel(r"Wavelength ($\AA$)")
		ylabel(r"Cross section (cm$^2$)")
		axis([enetoang(emax),enetoang(emin),cmin,cmax])
		legend(loc="best")
		savefig(elfig)


if __name__=="__main__":
	if(len(sys.argv)<2):
		print "Please give a file name"
		sys.exit()

	filename=sys.argv[1]
	print "Your file name : ",filename
	PlotFile(filename)
	show()
	
	

	


