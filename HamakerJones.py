#!/usr/bin/python
from subprocess import Popen,PIPE
import numpy as np

def convertCppCmplx(cmplxStr):
    cmplxStr=cmplxStr[1:-1]
    (realStr,imagStr)=cmplxStr.split(',')
    cmplxFlt=float(realStr) +1j*float(imagStr)
    return(cmplxFlt)

def getJonesByAzEl(freq,az,el):
    p = Popen(["./ElementResponseMain", 
      str(freq),str(az),str(el)],
      stdout=PIPE,close_fds=True)
    JonesOut = p.communicate()[0].splitlines()
    (JonesStr11,JonesStr12)=JonesOut[0].split(';')
    (JonesStr21,JonesStr22)=JonesOut[1].split(';')
    JXtheta=  convertCppCmplx(JonesStr11)
    JXphi=convertCppCmplx(JonesStr12)
    JYtheta=  convertCppCmplx(JonesStr21)
    JYphi=convertCppCmplx(JonesStr22)
    return JXphi,JXtheta,JYphi,JYtheta

if __name__ == "__main__":
   
    freq=3e7
    az=np.arange(0.0,2*np.pi,0.1)
    el=np.arange(0.0,np.pi,0.1)
    #azs,els=np.meshgrid(az,el)
    print "set hidden3d"
    print "splot '-' with lines"
    for azi in az:
      for eli in el:
         JXphi,JXtheta,JYphi,JYtheta=getJonesByAzEl(freq,azi,eli)
         print azi,eli,JXphi.real
      print ""
