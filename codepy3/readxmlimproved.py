#!/usr/bin/env python
# -*- coding:Utf-8 -*-

from powerlaw import powerlaw
from MathFun import intloglog,extrapolate
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
from itertools import cycle
colo=['b','g','r','c','m','y','k']
das=[(None,None),(4,1,1,1),(2,1), (5,1,2,1),(4,1),(6,1), (10,1)]
lsty=[  dict(list(zip(('color','dashes'),(col,sty)))) for sty in das for col in colo  ]
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




def PlotProtonElasticNode(vNode,bRfile,vIsEnergyTbl=False,vEtbl=[]):
    """ Function to plot a proton node based on the N2 data of Kozelov and Ivanov 1992"""
# 	// from Kozelov & Ivanov 1992, eq. (1), Table 1 for N2
#	double a[5];
#	a[0]=5.3746;
#	a[1]=-0.20842;
#	a[2]=0.40561E-2;
#	a[3]=-0.33036E-2;
#	a[4]=-0.67680E-3;
    a = [5.3746,-0.20842,0.40561E-2,-0.33036E-2,-0.67680E-3]

    ene=array(powerlaw(1E3,1E7,100))
    if vIsEnergyTbl:
        ene=array(vEtbl)

    elastic = zeros((len(ene)))
    for i in range(len(ene)):
        kd = 0
        if ene[i] < 2E5:
            for j in range(5):
                kd += a[j] * log(ene[i]*1E-3 )**j
            elastic[i] = 1E-16 * 0.529 **2 * exp(kd)
        else:
        # mElasticCrscm2[i] = 1E-16*0.529*0.529*exp(8.658 -  log(mGrideV[i]*1E-3));
            elastic[i] = 1E-16*0.529*0.529*exp(8.658 -  log(ene[i]*1E-3))

    plot(ene,elastic,label="Elastic")

def PlotHydrogenElasticNode(vNode,bRfile,vIsEnergyTbl=False,vEtbl=[]):
    a = [4.183,-0.39348,0.65156E-2,0.30656E-2,0.59441E-3]
    ene=array(powerlaw(1E3,1E7,100))
    if vIsEnergyTbl:
        ene=array(vEtbl)

    elastic = zeros((len(ene)))
    for i in range(len(ene)):
        kd = 0
        if ene[i] < 2E5:
            for j in range(5):
                kd += a[j]* log(ene[i] *1E-3 )**j
            elastic[i] = 1E-16 * 0.529 **2 * exp(kd)
        else:
            elastic[i] = 1E-16*0.529*0.529*exp(6.147 - 0.7 *log(ene[i]*1E-3))


    plot(ene,elastic,label="Elastic")


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
 #       print("We have a shirai N2", Emin, Emax, threashold, tip, params)
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


def PlotSinghalNode(vCh, bRfile, vIsEnergyTbl=False,vEtbl=[]):
    """Plot energy vs. cross-section line plot for species with Singhal coefficients
    
    Arguments
    ----------
    vNode: flexible container object from ElementTree [output of root.findall function]
    
    Returns
    -------
    
    None
    
    Notes
    -------
    
    The legend, labels, and error bar information are found in the xml database. Plot created is output in a seperate window..
    
    """
    parray = [float(k) for k in vCh.find("params").text.split()]
    leg = vCh.find("Legend").text
    OmegaP = False
    CtypeP = False
    Excite = False
    Ion = False

    if bRfile:
        if(0!=len(vCh.findall("Proc"))):
            leg=vCh.find("Proc").text

    if vCh.find("Omega") is None:
        OmegaP = False 
    else:
        OmegaP = True
        
    if vCh.find("Ctype") is None:
        CtypeP = False
    else:
        CtypeP = True
        
    if vCh.find("Excitation") is None:
        Excite = False
    else:
        Excite = True
    
    if vCh.find("Ionization") is None:
        Ion = False 
    else:
        Ion = True        
    E=array(powerlaw(0.1,10000,100))
    if(vIsEnergyTbl):
        E = array(vEtbl)
    Cross = zeros((len(E)))
    print(type(E))
    q0=6.513E-14
    if Excite:
        if OmegaP:
                ratio = parray[0] / E
                Term1 = (q0 * parray[5])/(parray[0] * parray[0])
                Term2 = (1 - (ratio ** (parray[1]))) ** parray[2]
                Term3 = ratio**parray[4]
                Cross = Term1 * Term2 * Term3
    
        if CtypeP:
                ratio = parray[0] / E
                Term1 = (q0 * parray[5]) / (E * parray[0])
                Term2 = (1-(ratio ** (parray[1]))) ** parray[2]
                Term3 = log(exp(1) + (4.0 * parray[4] / ratio))
                Cross = Term1 * Term2 * Term3
    
    if Ion:
            Sigma_0 = 1E-16
            A_E = (parray[1] / (E + parray[2])) * log((E / parray[3]) + parray[4]+(parray[5] / E))
            Gamma = (parray[6] * E) / (E + parray[7])
            T_0 = parray[8] - (parray[9] / (E + parray[10]))
            T_m = 0.5*(E - parray[0])
            Term1 = A_E
            Term2 = Gamma
            Term3 = arctan((T_m-T_0) / Gamma) + arctan(T_0 / Gamma)
            Cross = Sigma_0 * Term1 * Term2 * Term3
    


    Cross[E < parray[0] ] = 0
    uncertainty=0
    datauncert=zeros((len(E)))
    try:
        uncertainty=float(vCh.find("uncertainty").text)
        if (uncertainty.find("%")):
            uncert=float(uncertainty.replace("%",""))/100
            print("Uncertainty factor",uncert)
            datauncert=Cross * uncert 
        else:
            datauncert = float(uncertainty) * Cross
    except:
        uncertainty = 0.5
        datauncert= Cross * 0.5 
        
    errorbar(E, Cross,yerr=datauncert,label=leg)




def PlotStdNode(vNode,bRfile,vIsLambda,vIsEnergyTbl=False,vEtbl=[], isphoton=False):
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
    print(leg)
    if("fact" in list(vNode.find("Egrid").keys())):
        fact=float(vNode.find("Egrid").attrib.get("fact"))
#   print "Votre facteur :",fact
#   print loadtxt(StringIO((vNode.find("Egrid").text).replace("\n"," ")))
    dataenergy=loadtxt(StringIO(vNode.find("Egrid").text.replace("\n"," ")))*fact
    fact=1
    if("fact" in list(vNode.find("Cross").keys())):
        fact=float(vNode.find("Cross").attrib.get("fact"))
    print("Votre facteur :",fact)
    #print(vNode.find("Cross").text.replace("\n"," ").strip())
    datacrs=loadtxt(StringIO(vNode.find("Cross").text.replace("\n"," ").strip()))*fact
#   print datacrs
    
    uncertainty=0
    datauncert=zeros((len(datacrs)))
    if("uncertainty" in list(vNode.find("Cross").keys())):
        uncertainty=vNode.find("Cross").attrib.get("uncertainty")
        if (uncertainty.find("%")):
            uncert=float(uncertainty.replace("%",""))/100.
            print("Uncertainty factor",uncert)
            datauncert=datacrs*uncert
        else:
            value=float(uncertainty)*fact
            print("Uncertainty value")
            datauncert=ones((len(datacrs)))*value
    if(vIsEnergyTbl):
        vEtbl=array(vEtbl)
    if(vIsLambda):
        if(vIsEnergyTbl):
            style = next(linestyleskwargs)
            print(style)
            if isphoton:
                testcrs=intloglog(dataenergy,datacrs,vEtbl)
            else:            
                testcrs=extrapolate(dataenergy,datacrs,vEtbl)
            if("threshold" in list(vNode.keys())):
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
            style = next(linestyleskwargs)
            print(style)
            if isphoton:
                testcrs=intloglog(dataenergy,datacrs,vEtbl)
            else:            
                testcrs=extrapolate(dataenergy,datacrs,vEtbl)
            #testcrs=extrapolate(dataenergy,datacrs,vEtbl)
            if("threshold" in list(vNode.keys())):
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
        print("Please provide the title /title")
        sys.exit()
    titre=root.find(".//title").text

    if (root.findall(".//Emin")==[]) or (root.findall(".//Emax")==[]) or (root.findall(".//Cmin")==[]) or (root.findall(".//Cmax")==[]):
        print("Please provide the extreme value /Emin /Emax /Cmin /Cmax")
        sys.exit()
    emin=float(root.find(".//Emin").text)
    emax=float(root.find(".//Emax").text)
    cmin=float(root.find(".//Cmin").text)
    cmax=float(root.find(".//Cmax").text)

    if root.findall(".//plotname")==[]:
        print("Please provide the name of the plot /plotname")
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
    isphoton=False

    
    
    titre,emin,emax,cmin,cmax,figname,lfig,efig,elfig,figmult=CheckRoot(root)
    
    if lfig!="":
        isphoton=True
    if figmult!=1:
        print("multiplication de la taille de votre figure")
        print(figmult)
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
        
    print("We found ",len(processlist),"processes")
    
    for proc in processlist:
        if(0 == len(proc.findall("Singhal"))):
            if(0==len(proc.findall("Shirai"))):
                PlotStdNode(proc,bRfile,False,False,[], isphoton)
            else:
                PlotShiraiNode(proc,bRfile,False,False)
        else:
            PlotSinghalNode(proc, bRfile)

    for proc in processlist2:
        if(0 == len(proc.findall("Singhal"))):
            if(0==len(proc.findall("Shirai"))):
                if(0==len(proc.findall("Zero"))):
                    print("We are searching for proton or hydrogen function")
                    if(0==len(proc.findall("use_proton_function"))):
                       PlotProtonElasticNode(proc,bRfile,False)
                    elif(0==len(proc.findall("use_hydrogen_function"))):
                       PlotHydrogenElasticNode(proc,bRfile, False)
                    else:
                       PlotStdNode(proc,bRfile,False,False,[], isphoton)
                else:
                    print("We have a non-defined process")
            else:
                PlotShiraiNode(proc,bRfile,False,False)
        else:
            PlotSinghalNode(proc, bRfile)
    for proc in processlist3:
        if(0 == len(proc.findall("Singhal"))):
            if(0==len(proc.findall("Shirai"))):
                PlotStdNode(proc,bRfile,False,False,[], isphoton)
            else:
                PlotShiraiNode(proc,bRfile,False)
        else:
            PlotSinghalNode(proc, bRfile)
    
    title(titre)
    xlabel("Energy (eV)")
    ylabel(r"Cross section (cm$^2$)")
    axis([emin,emax,cmin,cmax])
    if(abs(2*(emax-emin)/(emax+emin))<1):
        xticks(arange(emin,emax+(emax-emin)/4,(emax-emin)/4),arange(emin,emax+(emax-emin)/4,(emax-emin)/4))
    legend(loc="best")
    savefig(figname)
    isphoton=False

    if lfig!="":
        isphoton=True
        clf()
        linestyleskwargs = cycle(lsty)
        xscale("linear")
        yscale("log")
        print("creation of the lambda figure")
        for proc in processlist:
           if(0 == len(proc.findall("Singhal"))):
            if(0==len(proc.findall("Shirai"))):
                PlotStdNode(proc,bRfile,True,False,[], isphoton)
            else:
                PlotShiraiNode(proc,bRfile,True)
        for proc in processlist2:
            if(0 == len(proc.findall("Singhal"))):
                if(0==len(proc.findall("Shirai"))):
                    if(0==len(proc.findall("Zero"))):
                        print("We are searching for proton or hydrogen function")
                        if(0==len(proc.findall("use_proton_function"))):
                           PlotProtonElasticNode(proc,bRfile,False)
                        elif(0==len(proc.findall("use_hydrogen_function"))):
                           PlotHydrogenElasticNode(proc,bRfile, False)
                        else:
                           PlotStdNode(proc,bRfile,True,False,[], isphoton)
                    else:
                          print("We have a non-defined process")
                else:
                    PlotShiraiNode(proc,bRfile,True)
        for proc in processlist3:
            if(0 == len(proc.findall("Singhal"))):
                if(0==len(proc.findall("Shirai"))):
                    PlotStdNode(proc,bRfile,True,False,[], isphoton)
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
        #global linestyleskwargs
        linestyleskwargs = cycle(lsty)
        print("creation of the extrapolated figure")
        clf()
        xscale("log")
        yscale("log")
        print("creation of the lambda figure")
        for proc in processlist:
            if(0 == len(proc.findall("Singhal"))):
                if(0==len(proc.findall("Shirai"))):
                    PlotStdNode(proc,bRfile,False,True,Etbl, isphoton)
                else:
                    PlotShiraiNode(proc,bRfile,False,True,Etbl)
            else:
                PlotSinghalNode(proc, bRfile, True, Etbl)
        for proc in processlist2:
            if(0 == len(proc.findall("Singhal"))):
                if(0==len(proc.findall("Shirai"))):
                    if(0==len(proc.findall("Zero"))):
                        print("We are searching for proton or hydrogen function")
                        if(0==len(proc.findall("use_proton_function"))):
                           PlotProtonElasticNode(proc,bRfile, False)#True, Etbl)
                        elif(0==len(proc.findall("use_hydrogen_function"))):
                           PlotHydrogenElasticNode(proc,bRfile, False)#True, Etbl)
                        else:
                           PlotStdNode(proc,bRfile,False,True,Etbl, isphoton)
                    else:
                        print("We have a non-defined process")
                else:
                    PlotShiraiNode(proc,bRfile,False,True,Etbl)
            else:
                PlotSinghalNode(proc, bRfile, True, Etbl)
        for proc in processlist3:
            if(0 == len(proc.findall("Singhal"))):
                if(0==len(proc.findall("Shirai"))):
                    PlotStdNode(proc,bRfile,False,True,Etbl, isphoton)
                else:
                    PlotShiraiNode(proc,bRfile,False,True,Etbl)
            else:
                PlotSinghalNode(proc, bRfile, True, Etbl)
    
        title(titre)
        xlabel("Energy (eV)")
        ylabel(r"Cross section (cm$^2$)")
        axis([emin,emax,cmin,cmax])
        if(abs(2*(emax-emin)/(emax+emin))<1):
            xticks(arange(emin,emax+(emax-emin)/4,(emax-emin)/4),arange(emin,emax+(emax-emin)/4,(emax-emin)/4))
        legend(loc="best")
        savefig(efig)
    if elfig!="":
        print("creation of the extrapolated lambda figure")
        #global linestyleskwargs
        linestyleskwargs = cycle(lsty)
        clf()
        xscale("linear")
        yscale("log")
        print("creation of the lambda figure")
        for proc in processlist:
            if(0 == len(proc.findall("Singhal"))):
                if(0==len(proc.findall("Shirai"))):
                    PlotStdNode(proc,bRfile,True,True,Etbl, isphoton)
                else:
                    PlotShiraiNode(proc,bRfile,True,True,Etbl)
        for proc in processlist2:
            if(0 == len(proc.findall("Singhal"))):
                if(0==len(proc.findall("Shirai"))):
                    PlotStdNode(proc,bRfile,True,True,Etbl, isphoton)
                else:
                    PlotShiraiNode(proc,bRfile,True,True,Etbl)
        for proc in processlist3:
            if(0 == len(proc.findall("Singhal"))):
                if(0==len(proc.findall("Shirai"))):
                    PlotStdNode(proc,bRfile,True,True,Etbl, isphoton)
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
        print("Please give a file name")
        sys.exit()

    filename=sys.argv[1]
    print("Your file name : ",filename)
    PlotFile(filename)
    show()
    
    

    


