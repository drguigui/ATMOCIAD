#!/usr/bin/env python
# -*- coding:Utf-8 -*-

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
from io import StringIO
from string import *
from glob import glob
import sys
import os
from readxmlimproved import *

def GetUnique(root,str):
    processlist=root.findall(str)
    if(len(processlist)==1):
        return  processlist[0].text
    print("Error, impossible to find the node "+str+" in your file")
    sys.exit()


def ReadNode2(node,plotlist):
    mynode={}
    notes=GetUnique(node,"./Notes")
    if(type(notes)!=type("str")):
        notes=""
    mynode["Notes"]=notes
    nsource=node.find("Source")
    mynode["Source"]=nsource.text
    mynode["SourceType"]=nsource.attrib.get("type")

    uncertainty="???%"
    mynode["EstimatedUncertainty"]=True


    if(len(node.findall("texnote"))==1):
        texnote=node.find("texnote")
        mynode["TexNote"]=texnote.text
        mynode["TexTitle"]=texnote.attrib.get("title")
    else:
        mynode["TexNote"]=""
        mynode["TexTitle"]=""


    
    if(len(node.findall("uncertainty"))>0):
        uncertainty=node.find("uncertainty").text+"%"
        mynode["EstimatedUncertainty"]=False
    else:
        if("uncertainty" in list(node.find("Cross").keys())):
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

    stringplot="Fig. "
    for (i,r,m) in plotlist:
        stringplot+="\\ref{"+fnametoref(i)+"} "
    mynode["Plots"]=stringplot

    return mynode




def ReadNode(node,arr,plotlist):
    processus=GetUnique(node,"./Proc")
    #if(len(arr[processus]))==0:
    #   arr[processus]=[]
    if not processus in list(arr.keys()):
        arr[processus]=[]
    arr[processus].append(ReadNode2(node,plotlist))
    

def ReadFile(file,total,elastic,ionization,dissociation,excitation,emission,recommended,plotlist):
    """ Reads the file 'file', and fill the sections considering the data inside
    Returns the name of the species and the name of the collider to check the correctness of the database
    """
    print("READ THE FILE", file)
    root=ET.parse(file).getroot()
        
    specie=GetUnique(root,".//Name").strip()
    collider=GetUnique(root,".//Collider").strip()
    if(len(root.findall("./Final"))>0):
        return specie,collider
    recofile=False
    if(len(root.findall("./RecommendedFile"))>0):
        recofile=True
    if len(root.findall(".//title"))==0 :
        print("Please provide the title /title")
        exit()
    description_img="Cross sections for "+root.find(".//title").text
    figlist=[]
    if len(root.findall(".//plotname"))!=0:
                print("We are plotting", (root.find(".//plotname").text))
                PlotFile(file)
                figlist.append(   (root.find(".//plotname").text,description_img,recofile))
    if len(root.findall(".//lambplotname"))!=0:
        figlist.append(   (root.find(".//lambplotname").text,description_img+" (wavelength version)",recofile))
    if (len(root.findall(".//exlambplotname"))!=0):
        figlist.append(   (root.find(".//exlambplotname").text,description_img+" (with extrapolation version)",recofile))
    if len(root.findall(".//explotname"))!=0:
        figlist.append(   (root.find(".//explotname").text,description_img+" (wavelength with extrapolation version)",recofile))

    plotlist.append(figlist)


    
    processlist=root.findall(".//Process")
    print("Nous avons trouve ",len(processlist),"processus")
    ppp = 0
    for proc in processlist:
        print(ppp)
        ppp += 1
        if recofile:
            ReadNode(proc,recommended,figlist)
        else:
            section=GetUnique(proc,"./Section")
            if section=="ionization":
                ReadNode(proc,ionization,figlist)
            elif section=="dissociation":
                ReadNode(proc,dissociation,figlist)
            elif section=="excitation":
                ReadNode(proc,excitation,figlist)
            elif section=="emission":
                ReadNode(proc,emission,figlist)
            else:
                print(("Impossible to find you process section: "+section+" in file "+file))
                exit()


    processlist2=root.findall(".//ElasticCrs")
    for proc in processlist2:
        if recofile:
            recommended["Elastic"]=[ReadNode2(proc,figlist)]
        else:
            elastic.append(ReadNode2(proc,figlist))
    processlist3=root.findall(".//TotalCrs")
    for proc in processlist3:
        if recofile:
            recommended["Total"]=[ReadNode2(proc,figlist)]
        else:
            total.append(ReadNode2(proc,figlist))
    return specie,collider
    
    
def printnodedico(process,sp,coll,cap):
    
    description=""
    str="\\begin{sidewaystable}\n\centering \small\n\\begin{tabular}{|c||c|c|c|c|c|c|}\n \hline \n"
    #str+=("Process & Reference & Threshold       & Range of energy  &   Uncertainty      &    Notes    &  Properties \\\\ \n \hline \n")
    str+=("Process & Reference & Threshold       & Range of energy  &   Uncertainty         &  Properties & Plots \\\\ \n \hline \n")
#   print("^ Process    ^   Reference   ^ Threshold        ^ Range of energy  ^   Uncertainty      ^    Notes    ^  Properties  ^")
    for i in list(process.keys()):
        cnt=0
        if coll=="ph":
            coll="$\lambda$"
        procname = sp+" + "+coll+" $\\rightarrow$ "+i
        description+="\subsubsection{"+procname+"}\n"
        for j in process[i]:
            prc=""
            if cnt==0:
                prc=procname
                cnt=1
            props=""

            if ""!=j["TexNote"] and ""!= j["TexTitle"]:
                description+="\paragraph{"+j["TexTitle"]+"}\n"+j["TexNote"]+"\n\n"

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
            #str+=(" "+prc+" & "+type+"  "+j["Source"]+" & "+j["Threshold"]+" & "+j["Range"]+" & "+j["Uncertainty"]+ " & "+j["Notes"]+" & "+props+"\\\\ \n ")
            str+=cleanstr(" "+prc+" & "+type+"  "+j["Source"]+" & "+j["Threshold"]+" & "+j["Range"]+" & "+j["Uncertainty"]+" & "+props+" & ")+j["Plots"]+"\\\\ \n "
        str+="\hline\n"
    str+="\end{tabular}\n \caption{"+cap+" for "+coll +" impact on "+cleanstr(sp)+"}\n\end{sidewaystable}\n"
    str+=("\n")
    description+="\n"
    return description+str  

def printnodesimple(process,specie,collider,cap):
    str="\\begin{sidewaystable}\n\centering \small\n\\begin{tabular}{|c|c|c|c|c|c|}\n \hline \n"
    description=""
    #str+=(" Reference & Threshold       & Range of energy  &   Uncertainty      &    Notes    &  Properties \\\\ \n \hline \n")
    if collider=="ph":
        collider="$\lambda$"
    str+=(" Reference & Threshold       & Range of energy  &   Uncertainty     &  Properties & Fig \\\\ \n \hline \n")
    for j in process:
        props=""
        if ""!=j["TexNote"] and ""!= j["TexTitle"]:
            description+="\subsubsection{"+j["TexTitle"]+"}\n"+j["TexNote"]+"\n\n"
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
        #str+=(" "+type+"  "+j["Source"]+" & "+j["Threshold"]+" & "+j["Range"]+" & "+j["Uncertainty"]+ " & "+j["Notes"]+" & "+props+"\\\\ \n ")
        str+=cleanstr(" "+type+"  "+j["Source"]+" & "+j["Threshold"]+" & "+j["Range"]+" & "+j["Uncertainty"]+" & "+props+" & ")+j["Plots"]+"\\\\ \n "
        str+="\hline\n"
    str+="\end{tabular}\n \caption{"+cap+" for "+collider +" impact on "+cleanstr(specie)+"} \n\end{sidewaystable}\n"
    str+=("\n")

    return str
        

def colorisestr(str,color):
    m=str.replace(">","$>$").replace("<","$<$")
    return "\color{"+color+"}{"+m+"}"
    
def fnametoref(stri):
    stri=stri.strip()
    stri=stri.replace("_","00")
    stri=stri.replace(".pdf","ppdf")
    return stri

def Printdatabase(specie,collider,total,elastic,ionization,dissociation,excitation,emission,recommended,plotlist):
    specie=cleanstr(specie)
    collider=cleanstr(collider)
    str=("\n\clearpage\n\section{Cross section of "+collider+" impact with "+specie+"}\n")

    figstr="\clearpage"
    figcount=0
    figrecomm=""
    for f in plotlist:
        for (fig,descr,recomm) in f:
            #figstr+=" \\begin{figure}\n \centering\n \includegraphics[width=10cm]{"+fig+"}\n  \label{"+fnametoref(fig)+"} \caption{Cross section for species "+specie+"} \n \end{figure}\n"
            if not recomm:
                figstr+=" \\begin{figure}\n\centering\n\includegraphics[width=12cm]{"+fig+"} \n\caption{"+descr+" }\n\label{"+fnametoref(fig)+"}\n\end{figure}\n"
                figcount+=1
                if(figcount==2):
                    figstr+="\clearpage\n"
                    figcount=0
            else:
                figrecomm+=" \\begin{figure}\n\centering\n\includegraphics[width=23cm,angle=-90]{"+fig+"} \n\caption{"+descr+" }\n\label{"+fnametoref(fig)+"}\n\end{figure}\n"
                figrecomm+="\clearpage\n"
    figstr+=figrecomm
    if len(total)!=0:
        str+=("\subsection{Total Cross Section}\n")
        str+=printnodesimple(total,specie,collider,"Total cross section")
    if len(elastic)!=0:
        str+=("\subsection{Elastic Cross Section}\n")
        str+=printnodesimple(elastic,specie,collider,"Elastic cross section")
    str+="\subsection{Inelastic Cross Sections}"
    if len(ionization)!=0:
        str+=("\subsubsection{Ionization Cross Sections}\n")
        str+=printnodedico(ionization,specie,collider,"Ionization Cross section")
    if len(dissociation)!=0:
        str+=("\subsubsection{Dissociation Cross Sections}\n")
        str+=printnodedico(dissociation,specie,collider,"Dissociation Cross section")
    if len(excitation)!=0:
        str+=("\subsubsection{Excitation Cross Sections}\n")
        str+=printnodedico(excitation,specie,collider,"Excitation Cross section")
    if len(emission)!=0:
        str+=("\subsection{Emission Cross Sections}\n")
        str+=printnodedico(emission,specie,collider,"Emission Cross section")
    if len(recommended)!=0:
        str+=("\subsection{Recommended data set}\n")
        str+=printnodedico(recommended,specie,collider,"Recommended Cross section")
    return str,figstr
def PrintLegend():
    str=("\subsubsection{Legend for the properties}\n")
    str+=("\paragraph{ R }: Recommended cross section for the processus. It is used in the main file. The selection of the recommended cross section is based on the quality of the data (e.g. errorbars, comparison with other experiments), the possibility of extrapolation, and the origin of the work, coupled with the consistency (sum of recommended cross sections ~ Total cross section)\n")
    str+=("\paragraph{ U }: Estimated uncertainty: sometimes, the uncertainty is not given, because of theoretical work... The authors of the database have to estimate the uncertainty, but the quality of that estimation can be questionable. Moreover, when data from different sources have been adapted (e.g. for extrapolation), the uncertainty can be modified...\n")
    str+=("\paragraph{ E }: Validated for extrapolation: the extrapolation of these cross sections is plausible. For example, when an analytic function has been applied...\n")
    return str


def PrintFolder(folder):
    files=glob(folder+"/*xml")
    print(("Your files:", files))
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
    plotlist=[]
    for i in files:
        print(("Processing file "+i))
    #   PlotFile(i)
        (newsp,newcoll)=ReadFile(i,total,elastic,ionization,dissociation,excitation,emission,recommended,plotlist)
        print(i,newsp,newcoll)
        if newsp=="":
            print("error in the file ",i," the specie name is badly set")
        if specie=="":
            specie=newsp
            collider=newcoll
        if specie!= newsp or collider !=newcoll:
            print("error : directory not homogenous. different species and/or colliders")
            print("Fatal error in file :",i)
            sys.exit()
    print("We are printing the database")
    str,figstr=Printdatabase(specie,collider,total,elastic,ionization,dissociation,excitation,emission,recommended,plotlist)
    str+="\n\n"
    str+=PrintLegend()

    str+="\n"+figstr+"\n\n\n"

    fna=os.path.basename(folder)
    filename=fna+".tex"
    f=open(filename,"w")
    f.writelines(str)
    f.close()

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
if __name__=="__main__":

    dirs=[i for i in glob("../*") if os.path.isdir(i) and i!="../resultat" and i!="../tmp"]
    print(dirs)

    for d in dirs:
        print("We are working on the folder",d)
        PrintFolder(d)
        print("Done")








