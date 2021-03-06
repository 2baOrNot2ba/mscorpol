#!/usr/bin/python
from subprocess import Popen,PIPE
import numpy as np

def convertCppCmplx(cmplxStr):
    cmplxStr=cmplxStr[1:-1]
    (realStr,imagStr)=cmplxStr.split(',')
    cmplxFlt=float(realStr) +1j*float(imagStr)
    return(cmplxFlt)

def getJonesByAzEl(freq,az,el):
    try:
      p = Popen(["../lofar_element_response/ElementResponseMain", 
          str(freq),str(az),str(el)],
          stdout=PIPE,close_fds=True)
    except OSError:
      print "Could not run 'ElementResponseMain'"
      exit(1)
    JonesOut = p.communicate()[0].splitlines()
    (JonesStr11,JonesStr12)=JonesOut[0].split(';')
    (JonesStr21,JonesStr22)=JonesOut[1].split(';')
    JXtheta=  convertCppCmplx(JonesStr11)
    JXphi=convertCppCmplx(JonesStr12)
    JYtheta=  convertCppCmplx(JonesStr21)
    JYphi=convertCppCmplx(JonesStr22)
    return JXphi,JXtheta,JYphi,JYtheta

if __name__ == "__main__":
   
    freq=8e7
    az=np.arange(0.0,2*np.pi,0.1)
    el=np.arange(0.0,np.pi/2,0.1)
    #azs,els=np.meshgrid(az,el)
    print "set hidden3d"
    print "splot '-' with lines"
    for azi in az:
      for eli in el:
         JXphi,JXtheta,JYphi,JYtheta=getJonesByAzEl(freq,azi,eli)
         print azi,eli,JXtheta.real
      print ""
