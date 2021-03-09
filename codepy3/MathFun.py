#!/usr/bin/env python
# -*-Coding:Utf-8-*-

""" Various interpolation and math functions needed for the other Python files in this package

The code contains various math and interpolation functions needed to successfully run the other Python code in this directory

Examples
----------

To use function, type in the following on the command line:

    python MathFun.py 

Make sure this Python code is in the same directory as the other Python scripts in this package. It is used by some of these scripts. 

Arguments
----------

None

Returns
-------

None

Notes
--------
    
Modules needed: pylab, powerlaw (defined in this package), scipy.optimize (leastsq)
Notes created by B Hegyi, 10/26/20

"""

from pylab import *
from powerlaw import powerlaw
try:
	from scipy.optimize import leastsq
except:
	print("leastsq not working!")
#from scipy import interpolate
#from MathFlux import indirac, ingauss, inmaxw, normflux
#from MathGrid import gridpow,gridexp,widthgrid,gridcst,gridpolo
#from PrintVal import print_array,str_array



#def lininterp1(oldx,oldy,newx):
#	""" Interpolation with the scipy function. When values are out of borders, the interpolation is filled with 0, not extrapolated"""
#	f=interpolate.interp1d(oldx,oldy,bounds_error=False)
#	return f(newx)



def array_reverse(a):
    """To reverse an array
    
    Arguments
    ----------

    a: Numpy array
        Any array of values
    
    Returns
    -------
    
    Returns array, reversed from input 
    
    """
    
    #This is an array reversal function
    
    return a[::-1]



def intlin(goldx,goldy,gnewx):
    """Linear interpolation and extrapolation used in Trans*. 
	When the values are out of border, there is a linear extrapolation
	
    Arguments
    ----------

    goldx: Numpy array
        Any 1D array of x-values [input]
    goldy: Numpy array
        Any 1D array of y-values [input]
    gnewx: Numpy array
        Any 1D array of desired of x-values [output]

    
    Returns
    -------
    
    newy: A 1D array of interpolated/extrapolated y-values 
    
    """
    oldx=goldx
    oldy=goldy
    newx=gnewx
    reorder=False
    newy=zeros(len(newx))
    
    if(len(oldx)==0 or len(newx)==0):
        return
    
    if oldx[0]>oldx[-1]:
        oldx=array_reverse(oldx)
        oldy=array_reverse(oldy)
	
    if newx[0]>newx[-1]:
        newx=array_reverse(newx)
        reorder=True

    nin=len(oldx)
    for i in range(len(newx)):
        j=0
        pasfini=True
        
        while ((j<nin) and pasfini):
            if newx[i]<oldx[0]:
                newy[i]=oldy[1]-(oldx[1]-newx[i])*(oldy[1]-oldy[0])/(oldx[1]-oldx[0])
                pasfini=False
            elif  newx[i]>oldx[-1]:
                newy[i]=oldy[-1]-(oldx[-1]-newx[i])*(oldy[-1]-oldy[-2])/(oldx[-1]-oldx[-2])
                pasfini=False
            elif newx[i]==oldx[j]:
                newy[i]=oldy[j]
                pasfini=False
            elif newx[i]<oldx[j]:
                if j==0:
                    newy[i]=oldy[j]-(oldx[j]-newx[i])*(oldy[1]-oldy[0])/(oldx[1]-oldx[0])
                else:
                    newy[i]=oldy[j]-(oldx[j]-newx[i])*(oldy[j]-oldy[j-1])/(oldx[j]-oldx[j-1])
                pasfini=False
            j+=1
            
    if reorder:
        return array_reverse(newy)
    return newy


def intlog(oldx,oldy,newx):
	"""Logarithmic interpolation and extrapolation
    
    Arguments
    ----------

    oldx: Numpy array
        Any 1D array of x-values [input]
    oldy: Numpy array
        Any 1D array of y-values [input]
    newx: Numpy array
        Any 1D array of desired of x-values [output]

    
    Returns
    -------
    
    A 1D array of interpolated/extrapolated y-values, using a base-10 log interpolation (log-10 taken of x-values) 
    
    
    """
	return 10**(intlin(oldx,log10(oldy),newx))

def intloglog(oldx,oldy,newx):
	"""Logarithmic interpolation and extrapolation
    
    Arguments
    ----------

    oldx: Numpy array
        Any 1D array of x-values [input]
    oldy: Numpy array
        Any 1D array of y-values [input]
    newx: Numpy array
        Any 1D array of desired of x-values [output]

    
    Returns
    -------
    
    A 1D array of interpolated/extrapolated y-values, using a base-10 log interpolation (log-10 taken of both x and y values) 
    
    """
    
	return 10**(intlin(log10(oldx),log10(oldy),log10(newx)))


def splinterp(oldx,oldy,newx,k=3):
	""" Scipy spline interpolation
    
    Arguments
    ----------

    oldx: Numpy array
        Any 1D array of x-values [input]
    oldy: Numpy array
        Any 1D array of y-values [input]
    newx: Numpy array
        Any 1D array of desired of x-values [output]

    
    Returns
    -------
    
    New y-values from a spline interpolation, using new x-values  
    
    """
    
	sp=UnivariateSpline(oldx,oldy,k=k)
	return sp(newx)


def BetheOp(pvals,z):
    
    """ Function using the Bethe Oppenheimer process for interpolation 
    
    Arguments
    ----------

    pvals: Numpy array (3 entries)
    
    z: Numpy array


    
    Returns
    -------
    
    Bethe Oppenheimer values (tm1*tm2)  
    
    """
    
    A=pvals[0]
    B=pvals[1]
    C=pvals[2]
	#return A/z*log(B*z)
    tm1=A/z**C
    tm2=log(B*z)
	#print A,len(tm1),B,len(tm2)a
    if(len(tm2)!=len(tm1)):
        tm2=zeros((len(tm1)))
        for i in range(len(tm1)):
            tm2[i]=log(B*z[i])
        return tm1*tm2

def BetheMin(pvals,z,compar):
    
    """ Comparison between Bethe Oppenheimer values and another array
    
    Arguments
    ----------

    pvals: Numpy array (3 entries)
    
    z: Numpy array
    
    compar: Numpy array
    
        comparison array


    
    Returns
    -------
    
    Difference between comparison array and Bethe Oppenheimer array 
    
    """
    
    return (compar)-(BetheOp(pvals,z))


def intbetheoppenheimer(oldx,oldy,newx):
	"""We perform a least square fit to retrieve the A and B parameters of the 
	Q=A E-1 ln(B E) equation. And then we interpolate the data through this equation
    
    
    Arguments
    ----------

    oldx: Numpy array
        Any 1D array of x-values [input]
    oldy: Numpy array
        Any 1D array of y-values [input]
    newx: Numpy array
        Any 1D array of desired of x-values [output]


    
    Returns
    -------
    
    A new array of y-values using the Bethe Oppenheimer interpolation 
     
    
    """
	pretrieve=array([1E-14,1E-3,1.])
	pleast=leastsq(BetheMin,pretrieve,args=(oldx,oldy))
	print("Bethe Oppenheimer parameters : ",pleast)
	pvals=pleast[0]
	return  BetheOp(pvals,newx) #A/newx*ln(B*newx)
	#return 10**(intlin(log10(oldx),log10(oldy),log10(newx)))



def gaussangles(m):
	"""Compute gaussian angles
    
    Arguments
    ----------

    m : integer
        number of Gaussian angles
    
        
    Returns
    -------
    
    A dictionary with Gaussian angles and weights
    
    """
    
	tol=1E-15
	if m<5:
		tol=1E-30 # truc pour la precision
	en=float(m)
	np1=m+1
	nnp1=m*np1
	cona=float(m-1)/float(8*m**3)
	lim=m/2+1 #attention
	gmu=zeros(m,float)
	gwt=zeros(m,float)
	if m<1:
		gmu[0]=0.5
		gwt[0]=1
		return
	for k in range(1,lim):
		t=(4*k-1)*pi/float(4*m+2)
		x=cos(t+cona/tan(t))
		xi =x
		ggtop=-1
		tmp=0.
		pm2=1.
		while (abs(xi-x)>tol or ggtop==-1):
			x=xi
			ggtop=ggtop+1
			pm2=1.
			pm1=x
			p=0.
			for nn in range(2,m+1):
				p=((2*nn-1)*x*pm1-(nn-1)*pm2)/nn
				pm2=pm1
				pm1=p
			tmp=1./(1.-x**2)
			ppr=en*(pm2-x*p)*tmp
			p2pri=(2.*x*ppr-nnp1*p)*tmp
			xi=x-(p/ppr)*(1.+(p/ppr)*p2pri/(2.*ppr))
		gmu[k-1]=-x
		gwt[k-1]=2./(tmp*(en*pm2)**2)
		gmu[m-k]=-gmu[k-1]
		gwt[m-k]=gwt[k-1]
	if m%2!=0:
		gmu[lim]=0
		prod=1.
		for  k in range(3,m,2):
			prod*=float(k)/float(k-1)
		gwt[lim]=2./prod**2
	truc=gmu.copy()
	for i in range(0,len(gmu)):
		gmu[i]=0.5*gmu[i]+0.5
		gwt[i]=0.5*gwt[i]
		truc[i]=acos(gmu[i])*180/pi
	return dict(list(zip(truc,gwt)))


def gaussa(m):
	""" Return gaussian angle, weight
    
        Arguments
    ----------

    m : integer
        number of Gaussian angles
    
        
    Returns
    -------
    
    A tuple with Gaussian angles and weights
    
    """
    
    
	tol=1E-15
	if m<5:
		tol=1E-30 # truc pour la precision
	en=float(m)
	np1=m+1
	nnp1=m*np1
	cona=float(m-1)/float(8*m**3)
	lim=m/2+1 #attention
	gmu=zeros(m,float)
	gwt=zeros(m,float)
	if m<1:
		gmu[0]=0.5
		gwt[0]=1
		return
	for k in range(1,lim):
		t=(4*k-1)*pi/float(4*m+2)
		x=cos(t+cona/tan(t))
		xi =x
		ggtop=-1
		tmp=0.
		pm2=1.
		while (abs(xi-x)>tol or ggtop==-1):
			x=xi
			ggtop=ggtop+1
			pm2=1.
			pm1=x
			p=0.
			for nn in range(2,m+1):
				p=((2*nn-1)*x*pm1-(nn-1)*pm2)/nn
				pm2=pm1
				pm1=p
			tmp=1./(1.-x**2)
			ppr=en*(pm2-x*p)*tmp
			p2pri=(2.*x*ppr-nnp1*p)*tmp
			xi=x-(p/ppr)*(1.+(p/ppr)*p2pri/(2.*ppr))
		gmu[k-1]=-x
		gwt[k-1]=2./(tmp*(en*pm2)**2)
		gmu[m-k]=-gmu[k-1]
		gwt[m-k]=gwt[k-1]
	if m%2!=0:
		gmu[lim]=0
		prod=1.
		for  k in range(3,m,2):
			prod*=float(k)/float(k-1)
		gwt[lim]=2./prod**2
	truc=gmu.copy()
	for i in range(0,len(gmu)):
		gmu[i]=0.5*gmu[i]+0.5
		gwt[i]=0.5*gwt[i]
	return (gmu,gwt)



class GaussianAngle:
	def __init__(self,nbangles):
		self.nbangles=nbangles
		self.nbmid=int(nbangles/2)
		self.inittbl()
	def inittbl(self):
		(gmu,gwt)=gaussa(self.nbmid)
		self.xmu=append(gmu[::-1],gmu*-1)
		self.weight=append(gwt[::-1],gwt)
		self.angzb=arccos(self.xmu)*180./pi







if __name__=="__main__":
	oldx=array(list(range(10)))
	newx=array(list(range(100)))/5-5
	oldy=oldx**2
	newy=intlog(oldx,oldy,newx)

	plot(oldx,oldy,label="old")
	plot(newx,newy,label="new")
	legend()
	show()
