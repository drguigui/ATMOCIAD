#!/usr/bin/env python
# -*- coding:Utf-8 -*-
from __future__ import division
from powerlaw import powerlaw
from MathFun import *
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


class NewShiraiCH4:
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






def PlotShiraiNode(vNode,teste):
	leg=""
	if(0==len(vNode.findall("legend"))):
		leg=vNode.attrib["name"]
	else:
		leg=vNode.find("legend").text
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

#	def __init__(self,emin,emax,threshold,eqtype,dataeq):
	if(tid=="CH4"):
		shirai=NewShiraiCH4(Emin,Emax,threshold,tip,params)
	if(tid=="CO2"):
		shirai=NewShiraiCO2(Emin,Emax,threshold,tip,params)
	if(tid=="N2"):
		shirai=NewShiraiN2(Emin,Emax,threshold,tip,params)

#	ene=arange(Emin,Emax,(Emax-Emin)/100.)
#	ene=array(powerlaw(Emin,Emax,100))
#	print type(ene)
#	cross=shirai.ReturnCrs(ene*1E-3)
	cross=shirai.ReturnCrs(array(teste)*1E-3)
	print "---------------------------------"
	print "NODE ",leg
	print "---------------------------------"
	for i in range(len(teste)):
		if(teste[i]<threshold):
			cross[i]=0
		if(teste[i]>1000):
			print teste[i],cross[i]
	print "---------------------------------"
	print "---------------------------------"
	print "*********************************"
	#newE=array([0.1,0.2,0.3,0.4,.5,0.6,0.7,0.8,0.9,1,1000,1250,1500,1750,2000,2500,3000,4000,5000,7500,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000])
	newE=array([0.1,0.2,0.3,0.4,.5,0.6,0.7,0.8,0.9,1,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500,3000,4000,5000,7500,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000])
	print newE
	print shirai.ReturnCrs(newE*1E-3)
	print "*********************************"
	print "---------------------------------"
	print "---------------------------------"
	
	datauncert=cross*uncertainty/100.
	errorbar(teste,cross,yerr=datauncert,label=leg)


def PlotStdNode(vNode,teste):
	leg=""
	if(0==len(vNode.findall("legend"))):
		try:
			leg=vNode.attrib["name"]
		except:
			leg="Elastic"
	else:
		leg=vNode.find("legend").text
	
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
			print "Facteur d'incertitude",uncert
			datauncert=datacrs*uncert
		else:
			value=float(uncertainty)*fact
			print "Valeur d'incertitude"
			datauncert=ones((len(datacrs)))*value

#	print datauncert

#	errorbar(datacrs,dataenergy,yerr=datauncert,label=leg)
	errorbar(dataenergy,datacrs,yerr=datauncert,label=leg)

	testcrs=intloglog(dataenergy,datacrs,teste)
	if("threshold" in vNode.keys()):
		threshold=float(vNode.attrib.get("threshold"))
		for i in range(len(teste)):
			if(teste[i]<threshold):
				testcrs[i]=0
	#errorbar(teste,testcrs,yerr=datauncert,label="MODIF"+leg)
	plot(teste,testcrs)
	#plot(teste,testcrs,label="MODIF"+leg)
	







if __name__=="__main__":
	if(len(sys.argv)<2):
		print "veuillez donner un nom de fichier"
		sys.exit()

	filename=sys.argv[1]
	print "Votre nom de fichier : ",filename
	
	
	
	
	xscale("log")
	yscale("log")
	
	root=ET.parse(filename).getroot()
	processlist=root.findall(".//Process")
	print "Nous avons trouve ",len(processlist),"processus"
	
	teste=powerlaw(0.1,100000,500)
	for proc in processlist:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc,teste)
		else:
			PlotShiraiNode(proc,teste)
	#	print proc.attrib["name"]
	#	print len(proc.findall("Ionization"))
	#	print proc.find("Cross").text

	processlist2=root.findall(".//ElasticCrs")
	for proc in processlist2:
		if(0==len(proc.findall("Shirai"))):
			PlotStdNode(proc,teste)
		else:
			PlotShiraiNode(proc,teste)
	
	processlist3=root.findall(".//TotalCrs")
	for proc in processlist3:
		PlotStdNode(proc,teste)
	
	
	
	title("Cross sections comparisons")
	xlabel("Energy [eV]")
	ylabel("Cross section [cm$^2$]")
	legend(loc="best")
	show()

	


