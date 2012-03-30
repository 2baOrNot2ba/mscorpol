#!/usr/bin/python
import numpy as np
from scipy import *
from scipy.special import *
from datetime import datetime, timedelta
from pyrap.measures import measures
from pyrap.quanta import quantity
import pyrap.tables as pt
from parseAntennaField import parseAntennaField

#Matrix transform from linear to circular basis (Hamaker1996b) for exp(i*w*t)
lin2circ =1/sqrt(2)*np.matrix('1.0 1.0*1j; 1.0    -1.0*1j')
#Inverse of previous
lin2circH=1/sqrt(2)*np.matrix('1.0 1.0;   -1.0*1j +1.0*1j')
#Effective lengths of electric dipoles along nominal LOFAR XY directions
LeffXY=np.matrix('1 0 0; 0 1 0')
#The rotation matrix from station coordinates to antenna XY coordinates:   
bisectRot = np.matrix([[-1/sqrt(2), -1/sqrt(2), 0],
                       [ 1/sqrt(2), -1/sqrt(2), 0],
                       [         0,          0, 1]])   

def getDipJones(obsTimes,stnPos,stnRot,srcDirection,
                doInvJ=False,doPolPrec=True,doCirc=False,
                showJones=False,doNormalize=True):

   obsTimesArr=obsTimes.get_value(); obsTimeUnit=obsTimes.get_unit()
   if doInvJ==True:
      InvJonesDip=zeros((len(obsTimesArr),2,2),dtype=complex)
   else:
      JonesDip=zeros((len(obsTimesArr),2,2),dtype=complex)

   me=measures()
   #Set position of reference frame w.r.t. ITRF
   me.doframe(stnPos)
   if doPolPrec:
       #Get sky precession rotation matrix
       #(Assuming no change over data interval)
       me.doframe(me.epoch('UTC',quantity(obsTimesArr[0],obsTimeUnit)))
       precMat=getSkyPrecessionMat(me,srcDirection)
   for ti in range(len(obsTimesArr)):
       #Set current time in reference frame 
       me.doframe(me.epoch('UTC',quantity(obsTimesArr[ti],obsTimeUnit)))
       #Convert phase ref dir to current polar ITRF
       ITRFdirang=me.measure(srcDirection,'ITRF')
       #Convert polar ITRF to cartesian ITRF
       raITRF=ITRFdirang['m0']['value']
       decITRF=ITRFdirang['m1']['value']
       lmnITRF =sph2crt(raITRF,decITRF)
       l=lmnITRF[0,0]
       m=lmnITRF[0,1]
       n=lmnITRF[0,2]

       #Compute polariz comps in spherical sys to cartesian Station coord sys
       polz2cart =1.0/sqrt(l*l+m*m)* np.matrix([
                          [ m, l*n],
                          [-l, m*n],
                          [ .0, l*l+m*m]])

       JonesDipMat=LeffXY*bisectRot * stnRot  * polz2cart;
       if doPolPrec:
          #With precession:
          JonesDipMat=JonesDipMat * precMat
          #else: 
          #Do not apply precession rotation of polarimetric frame.
          #This is then the apparent polarization frame.
       if doCirc:
          #Let the brightness components be in circular basis:
          JonesDipMat=JonesDipMat*lin2circH
          
       if showJones: 
          print JonesDipMat
       if doInvJ==True :
          IJonesDipMat=np.linalg.inv(JonesDipMat)
          if doNormalize:
             #JonesDipMat=JonesDipMat/sqrt(abs(np.linalg.det(JonesDipMat)))
             g=sqrt(np.matrix.trace(IJonesDipMat*IJonesDipMat.H)/2.0)
             IJonesDipMat=IJonesDipMat/g
          InvJonesDip[ti,:,:]=IJonesDipMat
       else:
          JonesDip[ti,:,:]=JonesDipMat
   #plotSVJonesGnuplot(InvJonesDip)
   if doInvJ==True:
      return InvJonesDip
   else:
      return JonesDip

def getSkyPrecessionMat(me,srcDirection):
       #Compute precession matrix. This is the 2D rotation along srcDirection
       #between J2000 and current epoch.

       #At J2000 epoch:
       #compute 2 cartesian vectors orthogonal to srcDirection: alpha & delta.
       #alpha is the orthogonal direction on equator:
       alphaJ2000ra=srcDirection['m0']['value']+pi/2
       alphaJ2000dec=0.0
       alphaJ2000vec=sph2crt(alphaJ2000ra,alphaJ2000dec)
       #delta is the orthogonal direction along the meridian.
       deltaJ2000ra=srcDirection['m0']['value']
       deltaJ2000dec=srcDirection['m1']['value']+pi/2
       deltaJ2000vec=sph2crt(deltaJ2000ra,deltaJ2000dec)

       alpha=me.direction('J2000',
                   str(alphaJ2000ra)+'rad',str(alphaJ2000dec)+'rad')
       delta=me.direction('J2000',
                   str(deltaJ2000ra)+'rad',str(deltaJ2000dec)+'rad')
       #Convert alpha & delta to directions in the current epoch
       alphaTru=me.measure(alpha,'JMEAN') #'JTRUE' isn't stable
       deltaTru=me.measure(delta,'JMEAN')
       raA=alphaTru['m0']['value']
       decA=alphaTru['m1']['value']
       alphaTruvec=sph2crt(raA,decA)
       raD=deltaTru['m0']['value']
       decD=deltaTru['m1']['value']
       deltaTruvec=sph2crt(raD,decD)

       #Precession of polarization basis is 
       precMat=np.matrix([ [(alphaTruvec*alphaJ2000vec.T)[0,0],
                            (alphaTruvec*deltaJ2000vec.T)[0,0]],
                           [(deltaTruvec*alphaJ2000vec.T)[0,0],
                            (deltaTruvec*deltaJ2000vec.T)[0,0]] ])
#       print precMat
#       print alphaTruvec*deltaTruvec.T
       return precMat


def sph2crt(azi,ele):
    #Spherical polar angles in azimuth and elevation to cartesian conversion
    x=cos(ele)*cos(azi)
    y=cos(ele)*sin(azi)
    z=sin(ele)
    return(np.matrix([x,y,z]))

def getDipJonesByAntFld(obsTimes,stnName,srcDirection,rcumode,lambda0):
   AntFld=parseAntennaField(stnName)
   stnLoc=stnName[0:2]
   if rcumode==3:
      AntBand='LBA'
   elif rcumode==5:
      if stnLoc=='CS' or stnLoc=='RS':
         AntBand='HBA0'
      else:
         AntBand='HBA'
   stnPos=np.matrix(AntFld[AntBand]['POSITION']).T
   stnRot=np.matrix(AntFld[AntBand]['ROTATION_MATRIX'])
   me=measures()
   stnPos_me=me.position('ITRF',str(stnPos[0,0])+'m',str(stnPos[1,0])+'m',str(stnPos[2,0])+'m')
   return getDipJones(obsTimes,stnPos_me,stnRot,srcDirection)

def plotSVJonesGnuplot(Jn):
   print 'plot "-" with lines'
   for ti in range(0,Jn.shape[0]):
       s=np.linalg.svd(Jn[ti,:,:].squeeze(),compute_uv=False)
       c=s[0]/s[1]
       cc=np.linalg.cond(Jn[ti,:,:].squeeze())
       #print ti, -20*log10((c+1)/(c-1) )
       print ti, s[0]+s[1]

def testJonesByAntFld():
   lambda0=2.20 #2.1
   stnName='UK608'
   beginTime=datetime(2011,10,24,18,0,0)
   endTime=datetime(2011,10,25,6,0,0)
   stepTime=timedelta(minutes=5)
   td=endTime-beginTime
   Times=[]
   for ti in range(0,td.seconds/stepTime.seconds):
       Times.append( (quantity((beginTime+ti*stepTime).isoformat())).get_value() )
       #obsTimes.append(beginTime+ti*stepTime)
   obsTimes=quantity(Times,'d') 
   srcDir=measures().direction('J2000','6.11378655886310rad','1.02191936586355rad')
   Jn=getDipJonesByAntFld(obsTimes,stnName,srcDir,5,lambda0)
   plotSVJonesGnuplot(Jn)

if __name__ == "__main__":
   testJonesByAntFld()
