#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.dates
import optparse
import numpy as np
from scipy import *
from scipy.special import *
from datetime import datetime, timedelta
from pyrap.measures import measures
from pyrap.quanta import quantity
import pyrap.tables as pt
from parseAntennaField import parseAntennaField
import parseParset
import dipoleJones
import HamakerJones

__version__="1.0"

bisectRot = np.matrix([[-1/sqrt(2), -1/sqrt(2), 0],
                       [ 1/sqrt(2), -1/sqrt(2), 0],
                       [         0,          0, 1]])

def sph2crtStn(azi,ele):
    #Spherical polar angles in azimuth and elevation to cartesian conversion
    x=cos(ele)*sin(azi)
    y=cos(ele)*cos(azi)
    z=sin(ele)
    return(np.matrix([x,y,z]))

def sph2crtGeo(lng,lat):
    #Spherical polar angles in azimuth and elevation to cartesian conversion
    x=cos(lat)*cos(lng)
    y=cos(lat)*sin(lng)
    z=sin(lat)
    return(np.matrix([x,y,z]))

def plotJonesGnuplot(startTime,stepTime,Jn):
   print "set title 'Jones values'"
   print "set xlabel 'Minutes since %s'" % startTime
   print "set xtics 60"
   print 'plot "-" with lines'
   timtick=0
   stepTime=stepTime/60
   for ti in range(0,Jn.shape[0]):
       s=np.linalg.svd(Jn[ti,:,:].squeeze(),compute_uv=False)
       c=s[0]/s[1]
       cc=np.linalg.cond(Jn[ti,:,:].squeeze())
       #print timtick, -20*log10((c+1)/(c-1) )
       #print timtick, (s[0]**2+s[1]**2)/2.0 #Stokes I gain
       print timtick, (Jn[ti,1,0]).real #Jones component
       timtick=timtick+stepTime

def getSkyPrecessionMat(me,srcDirection):
       #Compute precession matrix. This is the 2D rotation along srcDirection
       #between J2000 and current epoch.

       #At J2000 epoch:
       #compute 2 cartesian vectors orthogonal to srcDirection: alpha & delta.
       #alpha is the orthogonal direction on equator:
       alphaJ2000ra=srcDirection['m0']['value']+pi/2
       alphaJ2000dec=0.0
       alphaJ2000vec=sph2crtGeo(alphaJ2000ra,alphaJ2000dec)
       #delta is the orthogonal direction along the meridian.
       deltaJ2000ra=srcDirection['m0']['value']
       deltaJ2000dec=srcDirection['m1']['value']+pi/2
       deltaJ2000vec=sph2crtGeo(deltaJ2000ra,deltaJ2000dec)

       alpha=me.direction('J2000',
                   str(alphaJ2000ra)+'rad',str(alphaJ2000dec)+'rad')
       delta=me.direction('J2000',
                   str(deltaJ2000ra)+'rad',str(deltaJ2000dec)+'rad')
       #Convert alpha & delta to directions in the current epoch
       alphaTru=me.measure(alpha,'JTRUE') #'JTRUE' isn't stable
       deltaTru=me.measure(delta,'JTRUE')
       raA=alphaTru['m0']['value']
       decA=alphaTru['m1']['value']
       alphaTruvec=sph2crtGeo(raA,decA)
       raD=deltaTru['m0']['value']
       decD=deltaTru['m1']['value']
       deltaTruvec=sph2crtGeo(raD,decD)

       cosPrecAng=((alphaTruvec*alphaJ2000vec.T)[0,0]
                   +(deltaTruvec*deltaJ2000vec.T)[0,0])/2.0
       sinPrecAng=((alphaTruvec*deltaJ2000vec.T)[0,0]
                   -(deltaTruvec*alphaJ2000vec.T)[0,0])/2.0
       #Precession of polarization basis is 
       #precMat=np.matrix([ [(alphaTruvec*alphaJ2000vec.T)[0,0],
       #                     (alphaTruvec*deltaJ2000vec.T)[0,0]],
       #                    [(deltaTruvec*alphaJ2000vec.T)[0,0],
       #                     (deltaTruvec*deltaJ2000vec.T)[0,0]] ])
       precMat=np.matrix([[ cosPrecAng,sinPrecAng],
                          [-sinPrecAng,cosPrecAng]])
       #print precMat
       return precMat

def getHamakerJones(obsTimes,stnPos,stnRot,srcDirection,freq,
                doInvJ=False,doPolPrec=True,doCirc=False,
                showJones=False,doNormalize=False):

   obsTimesArr=obsTimes.get_value(); obsTimeUnit=obsTimes.get_unit()
   if doInvJ==True:
      InvJonesHam=zeros((len(obsTimesArr),2,2),dtype=complex)
   else:
      JonesHam=zeros((len(obsTimesArr),2,2),dtype=complex)

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
       #print quantity(obsTimesArr[ti],obsTimeUnit)
       azel=me.measure(srcDirection,'AZEL')
       az=azel['m0']['value']
       el=azel['m1']['value']
       phihatcrt=np.matrix([-cos(az),sin(az),0])
       thetahatcrt=np.matrix([sin(el)*sin(az),sin(el)*cos(az),-cos(el)])
       #print "AZ",azel['m0']['value']/pi*180
       #print "EL",azel['m1']['value']/pi*180
       #print azel['m0']['value']/pi*180,azel['m1']['value']/pi*180
       stnNorthOffsetAng=arccos(
                  -(stnRot[0,1]*stnRot[0,2]+stnRot[1,1]*stnRot[1,2])/(
                    sqrt(stnRot[0,1]**2+stnRot[1,1]**2)*
                    sqrt(stnRot[0,2]**2+stnRot[1,2]**2)))
       phi=-az+stnNorthOffsetAng-pi/2.0 #Added -pi/2
       theta=el
       (JXphi,JXtheta,JYphi,JYtheta)=HamakerJones.getJonesByAzEl(freq,phi,theta)
       JonesElemMat=np.matrix([[JXphi,JXtheta],[JYphi,JYtheta]])

       #Compute parallactic rotation

       #Convert phase ref dir to current polar ITRF
       srcITRFpolar=me.measure(srcDirection,'ITRF')
       #Convert polar ITRF to cartesian ITRF
       srcITRFlon=srcITRFpolar['m0']['value']
       srcITRFlat=srcITRFpolar['m1']['value']
       lmnITRF =sph2crtGeo(srcITRFlon,srcITRFlat)

       srcN_ITRFpolar=me.measure(
                  me.direction('J2000',
                      str(srcDirection['m0']['value'])+'rad', 
                      str(srcDirection['m1']['value']+math.pi/2.0)+'rad') 
              ,'ITRF')
       srcE_ITRFpolar=me.measure(
                  me.direction('J2000',
                      str(srcDirection['m0']['value']+math.pi/2.0)+'rad', 
                      str(0.0)+'rad') 
              ,'ITRF')
       north=sph2crtGeo(srcN_ITRFpolar['m0']['value'],srcN_ITRFpolar['m1']['value'])
       east=sph2crtGeo(srcE_ITRFpolar['m0']['value'],srcE_ITRFpolar['m1']['value'])
       #print "N*E", north*east.T
       local_pointing= ((stnRot.T* lmnITRF.T).T)
       local_north= ((stnRot.T* north.T).T)
       local_east= np.cross (local_north.squeeze(), local_pointing.squeeze())
       parallacticRot=np.matrix(
           [[(local_east*phihatcrt.T)[0,0],(local_north*phihatcrt.T)[0,0]],
            [(local_east*thetahatcrt.T)[0,0],(local_north*thetahatcrt.T)[0,0]]])
       #print local_pointing*thetahatcrt.T
       #parallacticRot=np.matrix(
       #    [[(local_east*thetahatcrt.T)[0,0],(local_north*thetahatcrt.T)[0,0]],
       #     [(local_east*phihatcrt.T)[0,0],(local_north*phihatcrt.T)[0,0]]])
       #print parallacticRot
       #bisectRot2D=bisectRot[0:2,0:2]
       #JonesHamMat=JonesElemMat*bisectRot2D * parallacticRot
       JonesHamMat=JonesElemMat * parallacticRot

       if doPolPrec:
          #With precession:
          JonesHamMat=JonesHamMat * precMat
          #else: 
          #Do not apply precession rotation of polarimetric frame.
          #This is then the apparent polarization frame.
       if doCirc:
          #Let the brightness components be in circular basis:
          JonesHamMat=JonesHamMat*lin2circH

       if showJones:
          print JonesHamMat
          pass
       if doInvJ==True :
          IJonesHamMat=np.linalg.inv(JonesHamMat)
          if doNormalize:
             #JonesDipMat=JonesDipMat/sqrt(abs(np.linalg.det(JonesDipMat)))
             g=sqrt(np.matrix.trace(IJonesHamMat*IJonesHamMat.H)/2.0)
             IJonesHamMat=IJonesHamMat/g
          InvJonesHam[ti,:,:]=IJonesHamMat
       else:
          JonesHam[ti,:,:]=JonesHamMat
   if doInvJ==True:
      return InvJonesHam
   else:
      return JonesHam


def getJonesByAntFld(model,obsTimes,stnName,srcDirection,freq=75.0E6):
   if freq > 100.0E6:
      rcumode=5
   else:
      rcumode=3
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
   #print stnPos[0,0],stnPos[1,0],stnPos[2,0]
   #print stnRot
   if model == "dipole":
      return dipoleJones.getDipJones(obsTimes,stnPos_me,stnRot,srcDirection)
   elif model == "Hamaker":
      return getHamakerJones(obsTimes,stnPos_me,stnRot,srcDirection,freq)
   else:
         error("Unknown antenna model (choose from 'dipole','Hamaker')")


def printJones(stnName,bTime,duration,stepTime,ra,dec,freqs,model):
      eTime = bTime+duration
      Times=[]
      for ti in range(0,duration.seconds/stepTime.seconds):
          Times.append( (quantity((bTime+ti*stepTime).isoformat())).get_value() )
          #obsTimes.append(beginTime+ti*stepTime)
      obsTimes=quantity(Times,'d')
      srcDir=measures().direction('J2000', ra,dec)
      print "Frequency","Time","J00","J01","J10","J11"
      for freq in freqs:
        Jn=getJonesByAntFld(model,obsTimes,stnName,srcDir,freq)

        #plotJonesGnuplot(quantity(obsTimes.get_value()[0],obsTimes.get_unit()).formatted("YMD"),stepTime.seconds,Jn)
        for ti in range(0,duration.seconds/stepTime.seconds):
          print freq, quantity(obsTimes.get_value()[ti],obsTimes.get_unit()).formatted("YMD"), Jn[ti,0,0], Jn[ti,0,1],Jn[ti,1,0],Jn[ti,1,1]

def plotJones(stnName,bTime,duration,stepTime,ra,dec,freqs,model):
    #frequencys=np.linspace(0,100e6,512)
    #freqs=frequencys[150:350]
    eTime = bTime+duration
    Times=[]
    print duration
    for ti in range(0,duration.seconds/stepTime.seconds):
        Times.append( (quantity((bTime+ti*stepTime).isoformat())).get_value() )
    obsTimes=quantity(Times,'d')
    srcDir=measures().direction('J2000', ra,dec)
    freq=freqs[0]
    Jn=getJonesByAntFld(model,obsTimes,stnName,srcDir,freq)
    p_ch = np.abs(Jn[:,0,0].squeeze())**2+np.abs(Jn[:,0,1].squeeze())**2
    q_ch = np.abs(Jn[:,1,1].squeeze())**2+np.abs(Jn[:,1,0].squeeze())**2
    plt.figure()
    plt.subplot(211)
    plt.plot(np.asarray(Times), 10*np.log10(p_ch))
    plt.title('p channel')
    #plt.clim(-9, -3)
    plt.subplot(212)
    plt.plot(np.asarray(Times), 10*np.log10(q_ch))
    plt.title('q-channel')
    #plt.clim(-9, -3)
    plt.xlabel('Time')
    plt.show()

def args2inpparms(args):
    stnName=args[0]
    bTime = datetime.strptime(args[1], "%Y-%m-%d %H:%M:%S")
    duration =timedelta(0,float(args[2]))
    stepTime =timedelta(0,float(args[3]))
    ra=args[4]+'rad'
    dec=args[5]+'rad'
    return stnName,bTime,duration,stepTime,ra,dec

if __name__ == "__main__":
   usage = "usage: %prog model [-o ObsID | stnName beginUTC duration timeStep pointingRA pointingDEC frequency]"
   #Example: $ antennaJones.py Hamaker SE607 '2012-04-01 01:02:03' 60 1 0 0 60E6
   #     or: $ antennaJones.py dipole -o L29053
   opt = optparse.OptionParser(usage=usage)

   opt.add_option('-o','--obsID', default="null", help='LOFAR Observation ID')
   options, args = opt.parse_args()
   if len(args) >0:
      model=args.pop(0)
   if options.obsID != "null" :
      args.append(' ')
      args.append(parseParset.getStartTime(options.obsID))
      args.append(parseParset.getDuration(options.obsID))
      args.append(parseParset.getStep(options.obsID))
      args.append(str(parseParset.getRA(options.obsID)))
      args.append(str(parseParset.getDEC(options.obsID)))
      stnName,bTime,duration,stepTime,ra,dec=args2inpparms(args)
      freqs=parseParset.getSubbandFreqs(options.obsID)
      #freqs=freqs[0:2]
      stnNames=parseParset.getStnNames(options.obsID)
      #stnNames=[stnNames[0]]
      for stnName in stnNames:
          print "Station: ", stnName
          printJones(stnName,bTime,duration,stepTime,ra,dec,freqs,model)
   elif len(args) == 7:
      stnName,bTime,duration,stepTime,ra,dec=args2inpparms(args)
      freqs=[float(args[6])]
      printJones(stnName,bTime,duration,stepTime,ra,dec,freqs,model)
      #plotJones(stnName,bTime,duration,stepTime,ra,dec,freqs,model)
   else :
      opt.error("incorrect number of arguments")
