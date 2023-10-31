#!/usr/bin/env python
# -*- coding:Utf-8 -*-

""" Creates plot of different cross-sections contained in xml file. Unlike angextrapole.py, plots cross sections as a function of energy (instead of wavelength). Extrapolation done using power law extrapolation of energies.
 
The purpose of this function is to create a plot of all cross sections for a certain species contained in the xml file from the Atmociad database. The function outputs the plot in a separate window after running the code.
If available, the plot shows the uncertainties as error bars for each cross section series. 

The extrapolation of each cross-section uses a power law function to create the extrapolation. See the documentation for powerlaw.py for details on the extrapolation function. The extrapolation is done until the cross section
values fall below a threshold value that is defined in the xml file.

Examples
----------

To use function, type in the following on the command line:
    
python extrapolexml.py filename

Arguments
----------
filename : str
    The file location of the xml file that contains the cross-sections

Returns
-------

None

Notes
--------

If an uncertainty estimate is included in the xml files, then error bars are placed on the line plots.
See comments in code for details about each function contained within. You can see the help for each function by running the following on the command line:
    pydoc ./extrapolxml.py

Notes created by B Hegyi, 10/19/20

"""

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
    
from pylab import *
from io import StringIO
from string import *


class NewShiraiCH4:
    """Processes cross-section data using methods from Shirai reference paper
                
    Arguments
    ----------
    emin : float
        Minimum energy of cross sections from xml file
    emax : float
        Maximum energy of cross sections from xml file
    threshold : float
        Minimum energy required for cross-section greater than zero    
    eqtype : float
        Type of equation  
    dataeq : float
        Shirai parameters; These parameters are used in the interpolation/extrapolation functions
    
    Returns
    -------
    
    Cross section values as a function of energy
    
    
    Reference
    -------
    
    Shirai T, Tabata T, Tawara H, and Itikawa Y. ANALYTIC CROSS SECTIONS FOR ELECTRON COLLISIONS WITH HYDROCARBONS: CH4, C2H6, C2H4, C2H2, C3H8, AND C3H6
    Atomic Data and Nuclear Data Tables, Volume 80, Issue 2, 147-204. https://doi.org/10.1006/adnd.2001.0878
    
    Notes
    -------
    
    The methods associated with this class are arranged as follows: init: The user inputs and defined constants. f1-f3: The basic fitting functions. S1-Sn: The combination of fitting functions used, set as a dictionary entry in ReturnCrs method. 
    The S-equation used is based on the self.eq parameter, which comes from the 'Equation' tag in the xml file.
    
    
    """    
    def __init__(self,emin,emax,threshold,eqtype,dataeq):
        
        """Initial parameters for NewShiraiCH4 method
        
        All of these parameters are input into the method


        Parameters
        ----------
        emin : float
            Minimum energy of cross sections from xml file
            emax : float
            Maximum energy of cross sections from xml file
            threshold : float
            Minimum energy required for cross-section greater than zero    
            eqtype : float
            Type of equation  
            dataeq : float
            Shirai parameters; These parameters are used in the interpolation/extrapolation functions
        """
        
        self.emin=emin #minimum energy to define cross sections
        self.emax=emax #maximum energy to define cross sections
        self.thresh=threshold/1000. # the threshold is put in KeV units for the computation
        self.eq=eqtype
        self.avals=dataeq
        self.sigma=1.E-16 # cm2, this sigma_0 in the paper
        self.Er=1.361E-2 # keV, this is the Rydberg constant
    def f1(self,E,c1,c2):
        return self.sigma*c1*(E/self.Er)**c2
    def f2(self,E,c1,c2,c3,c4):
        return self.f1(E,c1,c2)/(1+(E/c3)**(c2+c4))
    def f3(self,E,c1,c2,c3,c4,c5,c6):
        return self.f1(E,c1,c2)/(1+(E/c3)**(c2+c4)+(E/c5)**(c2+c6))

    # -- The following methods are fitting functions for the cross section, using a combination of {f1, f2, f3}
    # -- Advantages of these fitting functions vs. polynomials: ability to have smooth function near the data bound end points when extrapolating

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
    
    """Processes cross-section data for an interaction with an electron using methods from reference

    Arguments
    ----------
    emin : float
        Minimum energy of cross sections from xml file
    emax : float
        Maximum energy of cross sections from xml file
    threshold : float
        Minimum energy required for cross-section greater than zero    
    eqtype : float
        Type of equation  
    dataeq : float
        Shirai parameters (a1 - a4)

    Returns
    -------

    Cross section values as a function of wavelength

    Reference
    -------

    Shirai T, Tabata T, and Tawara H. Analytic Cross Sections for Electron Collisions with CO, CO2, and H2O Relevant to Edge Plasma Impurities, 
    Volume 79, Issue 1, 2001. 10.1006/adnd.2001.0866

    Notes
    -------
    
    The methods associated with this class are arranged as follows: init: The user inputs and defined constants. f1-f3: The basic fitting functions. S1-Sn: The combination of fitting functions used, set as a dictionary entry in ReturnCrs method. 
    The S-equation used is based on the self.eq parameter, which comes from the 'Equation' tag in the xml file.

    """
    
    def __init__(self,emin,emax,threshold,eqtype,dataeq):
        """Initial parameters for NewShiraiCH4 method
        
        All of these parameters are input into the method


        Parameters
        ----------
        emin : float
            Minimum energy of cross sections from xml file
            emax : float
            Maximum energy of cross sections from xml file
            threshold : float
            Minimum energy required for cross-section greater than zero    
            eqtype : float
            Type of equation  
            dataeq : float
            Shirai parameters; These parameters are used in the interpolation/extrapolation functions
        """
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
    
    # -- The following methods are fitting functions for the cross section, using a combination of {f1, f2, f3}
    # -- Advantages of these fitting functions vs. polynomials: ability to have smooth function near the data bound end points when extrapolating
    
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
    
    """Processes cross-section data for an interaction with an electron using methods from reference

    Arguments
    ----------
    emin : float
        Minimum energy of cross sections from xml file
    emax : float
        Maximum energy of cross sections from xml file
    threshold : float
        Minimum energy required for cross-section greater than zero    
    eqtype : float
        Type of equation  
    dataeq : float
        Shirai parameters 

    Returns
    -------

    Cross section values as a function of energy

    Reference
    -------

    Shirai T, Tabata T, and Tawara H. Analytic cross sections for electron impact collisions with nitrogen molecules,
    Volume 92, Issue 3, 2006. doi.org/10.1016/j.adt.2006.02.002

    Notes
    -------
    
    The methods associated with this class are arranged as follows: init: The user inputs and defined constants. f1-f3: The basic fitting functions. S1-Sn: The combination of fitting functions used, set as a dictionary entry in ReturnCrs method. 
    The S-equation used is based on the self.eq parameter, which comes from the 'Equation' tag in the xml file.
    
    

    """
    
    def __init__(self,emin,emax,threshold,eqtype,dataeq):
        """Initial parameters for NewShiraiCH4 method
        
        All of these parameters are input into the method


        Parameters
        ----------
        emin : float
            Minimum energy of cross sections from xml file
            emax : float
            Maximum energy of cross sections from xml file
            threshold : float
            Minimum energy required for cross-section greater than zero    
            eqtype : float
            Type of equation  
            dataeq : float
            Shirai parameters; These parameters are used in the interpolation/extrapolation functions
        """
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
    
    
    # -- The following methods are fitting functions for the cross section, using a combination of {f1, f2, f3}
    # -- Advantages of these fitting functions vs. polynomials: ability to have smooth function near the data bound end points when extrapolating
    
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






#def PlotShiraiNode(vNode,teste):
def PlotShiraiNode(vNode,bRfile,vIsLambda,vIsEnergyTbl=False,vEtbl=[]):
    """Plotenergy vs. cross-section line plot for species with Shirai coefficients
    
    
    Arguments
    ----------
    vNode: flexible container object from ElementTree [output of root.findall function]
    #ateste: array of energy values that will be used to calculate the cross sections
    
    Returns
    -------
    
    None
    
    Notes
    -------
    
    The legend, labels, and error bar information are found in the xml database. Plot created is output in a seperate window.
    
    """
    
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

   # cross=shirai.ReturnCrs(array(teste)*1E-3)
   # print("---------------------------------")
   # print("NODE ",leg)
   # print("---------------------------------")
    for i in range(len(ene)):
        if(teste[i]<threshold):
            cross[i]=0
   #     if(teste[i]>1000):
   #         print(teste[i],cross[i])
   # print("---------------------------------")
   # print("---------------------------------")
   # print("*********************************")
   # newE=array([0.1,0.2,0.3,0.4,.5,0.6,0.7,0.8,0.9,1,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500,3000,4000,5000,7500,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000])
   # print(newE)
   # print(shirai.ReturnCrs(newE*1E-3))
   # print("*********************************")
   # print("---------------------------------")
   # print("---------------------------------")
    
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





#def PlotStdNode(vNode,teste):
def PlotStdNode(vNode,bRfile,vIsLambda,vIsEnergyTbl=False,vEtbl=[]):
    """Plot wavelength vs. cross-section line plot for species with standard cross section values contained in the xml file
    
    
    Arguments
    ----------
    vNode: flexible container object from ElementTree [output of root.findall function]
    teste: array of energy values that will be used to calculate the cross sections
    
    Returns
    -------
    
    None
    
    Notes
    -------
    
    The legend, labels, and error bar information are found in the xml database. Plot created is output in a seperate window..
    """
    
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
    if("fact" in list(vNode.find("Egrid").keys())):
        fact=float(vNode.find("Egrid").attrib.get("fact"))
    dataenergy=loadtxt(StringIO(vNode.find("Egrid").text.replace("\n"," ")))*fact
    fact=1
    if("fact" in list(vNode.find("Cross").keys())):
        fact=float(vNode.find("Cross").attrib.get("fact"))
        print("Your factor :",fact)
    datacrs=loadtxt(StringIO(vNode.find("Cross").text.replace("\n"," ")))*fact

    
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
            print("Uncertainty Value")
            datauncert=ones((len(datacrs)))*value
    if(vIsEnergyTbl):
        vEtbl=array(vEtbl)
    if(vIsLambda):
        if(vIsEnergyTbl):
            #style = next(linestyleskwargs)
            #print(style)
            testcrs=extrapolate(dataenergy,datacrs,vEtbl)
            if("threshold" in list(vNode.keys())):
                threshold=float(vNode.attrib.get("threshold"))
                for i in range(len(vEtbl)):
                    if(vEtbl[i]<threshold):
                        testcrs[i]=0
            plot(enetoang(vEtbl),testcrs)
            errorbar(enetoang(dataenergy),datacrs,yerr=datauncert,label=leg)
        else:
            errorbar(enetoang(dataenergy),datacrs,yerr=datauncert,label=leg)

    else:
        if(vIsEnergyTbl):
#            style = next(linestyleskwargs)
#            print(style)
            testcrs=extrapolate(dataenergy,datacrs,vEtbl)
            if("threshold" in list(vNode.keys())):
                threshold=float(vNode.attrib.get("threshold"))
                for i in range(len(vEtbl)):
                    if(vEtbl[i]<threshold):
                        testcrs[i]=0
            plot((vEtbl),testcrs)
            errorbar(dataenergy,datacrs,yerr=datauncert,label=leg)
        else:
            errorbar(dataenergy,datacrs,yerr=datauncert,label=leg)
    


#    errorbar(dataenergy,datacrs,yerr=datauncert,label=leg)
#    testcrs=intloglog(dataenergy,datacrs,teste)
#    if("threshold" in list(vNode.keys())):
#        threshold=float(vNode.attrib.get("threshold"))
#        for i in range(len(teste)):
#            if(teste[i]<threshold):
#                testcrs[i]=0

 #   plot(teste,testcrs)






#--This is where the file name is taken from the command line when this .py script is run

if __name__=="__main__":
    if(len(sys.argv)<2):
        print("Please include a file name on the command line")
        exit()

    filename=sys.argv[1]
    print("File name : ",filename)
    
    
    
# -- This is where the plot is created  
    xscale("log")
    yscale("log")
    
    root=ET.parse(filename).getroot()
    processlist=root.findall(".//Process")
    print("We have found ",len(processlist),"processes")   #These are all of the processes that we find in the xml file
 
    bRfile=False
    if len(root.findall(".//RecommendedFile"))!=0:
        bRfile=True
           
    teste=powerlaw(0.1,100000,500)
    for proc in processlist:
        if(0 == len(proc.findall("Singhal"))):
            if(0==len(proc.findall("Shirai"))):
                PlotStdNode(proc,bRfile,False,True, teste)
            else:
                PlotShiraiNode(proc,bRfile,False,True, teste)
        else:
            PlotSinghalNode(proc, bRfile,True, teste)


    processlist2=root.findall(".//ElasticCrs")
    for proc in processlist2:
        if(0 == len(proc.findall("Singhal"))):
            if(0==len(proc.findall("Shirai"))):
                PlotStdNode(proc,bRfile,False,True, teste)
            else:
                PlotShiraiNode(proc,bRfile,False,True, teste)
        else:
            PlotSinghalNode(proc, bRfile,True, teste)
    processlist3=root.findall(".//TotalCrs")
    for proc in processlist3:
        if(0 == len(proc.findall("Singhal"))):
            if(0==len(proc.findall("Shirai"))):
                PlotStdNode(proc,bRfile,False,True, teste)
            else:
                PlotShiraiNode(proc,bRfile,False,True, teste)
        else:
            PlotSinghalNode(proc, bRfile,True, teste)
    
    if len(root.findall(".//title"))==0 :
        titre = "Cross section comparisons"
    else:
        titre= "Cross section comparisons" + root.find(".//title").text


    title(titre)
    xlabel("Energy [eV]")
    ylabel("Cross section [cm$^2$]")
    legend(loc="best")
    show()

    


