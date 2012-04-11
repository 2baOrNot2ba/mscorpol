#!/usr/bin/python
#TobiaC 2012-01-26
"""Correct polarization in LOFAR Measurement set (MS) data.

mscorpol can correct for the geometric projection of the incident
electric field onto the LOFAR antennas. The correction is made
for the center of the field and models the LOFAR response as
ideal electric dipoles. The output can be either linear or circular
polarized components. 
"""
import sys
import optparse
import numpy as np
from pyrap.measures import measures
from pyrap.quanta import quantity
import pyrap.tables as pt
from LOFARdipoleJones import getDipJones

__version__="1.4"

def correctMSforDipole(msfile):
  defaultVisConjOrder=True #Default order is conjugate(ANTENNA1)*ANTENNA2
  me=measures()
  mstab = pt.table(msfile,readonly=False,ack=True)
  #Get number of antennas
  NrOfAnts=msNrOfAnts(mstab)
  ta = pt.table(mstab.getkeyword('ANTENNA'),ack=False)
  pos = ta.getcol('POSITION')
  tl = pt.table(mstab.getkeyword('LOFAR_ANTENNA_FIELD'),ack=False)
  tf = pt.table(mstab.getkeyword('FIELD'),ack=False)
  srcDir = np.squeeze(tf.getcol('PHASE_DIR'))
  RA = quantity(srcDir[0],'rad'); dec = quantity(srcDir[1],'rad')
  srcDirection = me.direction('J2000',RA,dec)
  #Get unique list of times
  firstAntID=mstab.getcell("ANTENNA2",0)
  antUniq=mstab.query('ANTENNA2 == %d' % firstAntID)
  timevals = antUniq.getcol('TIME')
  #print timevals
  tt = quantity(timevals,'s')
  antUniq.close()
  antNr=0
  for tant in mstab.iter('ANTENNA1'):
      antID=int(tant.getcol('ANTENNA1')[0])
      print "ANTENNA=",antNr,'/',NrOfAnts
      x = quantity(pos[antID,0],'m');
      y = quantity(pos[antID,1],'m')
      z = quantity(pos[antID,2],'m')
      stnPos = me.position('ITRF',x,y,z)
      stnRot=tl.getcol('COORDINATE_AXES')[antID]
      JI=getDipJones(tt,stnPos,stnRot,srcDirection,
                doCirc=not(options.linear),doInvJ=True,doPolPrec=True,showJones=options.jones)

      #JI=np.ones((len(timevals),2,2))
      #Start processing DATA
      if not options.jones:
         #Apply this antennas Jones to DATA that matches value in ANTENNA2
         if defaultVisConjOrder:
            pass
         else:
            JI=np.conj(JI)
         tant2=mstab.query('ANTENNA2 == %d' % antID)
         data=tant2.getcol('DATA')    
         dataCohXY=np.array([[data[:,:,0],data[:,:,1]],[data[:,:,2],data[:,:,3]]])
         dataCor=np.zeros(dataCohXY.shape,dtype=np.complex)
         #Multiply Inverse Jones with visibility data.
         nrSampsPerAntTim=antNr+1
         for JIt in range(JI.shape[0]):
             for visTimInd in range(nrSampsPerAntTim):
                dataCor[:,:,JIt*nrSampsPerAntTim+visTimInd,:]=np.transpose(
                  np.tensordot(
                    dataCohXY[:,:,JIt*nrSampsPerAntTim+visTimInd,:],JI[JIt,:,:],
                  axes=[[1],[1]]),
                  (0,2,1))
         dataCorOut=np.transpose(
          [dataCor[0,0,:,:],dataCor[0,1,:,:],dataCor[1,0,:,:],dataCor[1,1,:,:]],
                              (1,2,0))
         tant2.putcol('DATA',dataCorOut)
         tant2.close()

         #Apply this antennas Jones matrix to DATA for ANTENNA1
         if defaultVisConjOrder:
            JI=np.conj(JI)
         else:
            JI=np.conj(JI) #Note: Jones inv has been conjugated once already.
         data=tant.getcol('DATA')    
         dataCohXY=np.array([[data[:,:,0],data[:,:,1]],[data[:,:,2],data[:,:,3]]])
         dataCor=np.zeros(dataCohXY.shape,dtype=np.complex)
         #Multiply Inverse Jones with visibility data.
         nrSampsPerAntTim=NrOfAnts-antNr
         for JIt in range(JI.shape[0]):
             for visTimInd in range(nrSampsPerAntTim):
                dataCor[:,:,JIt*nrSampsPerAntTim+visTimInd,:]=np.tensordot(
                    JI[JIt,:,:],dataCohXY[:,:,JIt*nrSampsPerAntTim+visTimInd,:],
                         axes=[[1],[0]])

         dataCorOut=np.transpose(
          [dataCor[0,0,:,:],dataCor[0,1,:,:],dataCor[1,0,:,:],dataCor[1,1,:,:]],
                              (1,2,0))
         tant.putcol('DATA',dataCorOut)
      antNr=antNr+1

  if not options.jones: updateMSmetadata(msfile)

def msNrOfAnts(mainTab):
  antslist=set([])
  for ants in mainTab.iter("ANTENNA1"):
      antslist.add(ants.getcell("ANTENNA1",0))
  for ants in mainTab.iter("ANTENNA2"):
      antslist.add(ants.getcell("ANTENNA2",0))
  return len(antslist)


def updateMSmetadata(msfile):
  #Update history to show that this script has modified original data
  tms = pt.table(msfile,readonly=False,ack=False)
  th = pt.table(tms.getkeyword('HISTORY'), readonly=False, ack=False)
  nr=th.nrows()
  th.addrows(1)
  tr=th.row()
  tr.put(nr,{'TIME': quantity('today').get('s').get_value(),
             'OBSERVATION_ID':0,
             'MESSAGE': 'Applied polarization corrections',
             'PRIORITY': 'NORMAL',
             'ORIGIN': '%s: version = %s' % (__file__,__version__),
             'OBJECT_ID':0, 
             'APPLICATION':__file__,
             'CLI_COMMAND':sys.argv,
             'APP_PARAMS': ['']})

  if not options.linear:
     #Change metadata information to be circular feeds
     feed = pt.table(tms.getkeyword('FEED'),readonly=False,ack=False)
     for tpart in feed.iter('ANTENNA_ID'):
      tpart.putcell('POLARIZATION_TYPE',0,['R','L'])

     polariz = pt.table(tms.getkeyword('POLARIZATION'),readonly=False,ack=False)
     polariz.putcell('CORR_TYPE',0,[5,6,7,8])
     tms.close()

    
if __name__ == "__main__":
   opt = optparse.OptionParser()
   opt.add_option('-f','--msfile',help='MS file')
   opt.add_option('-j','--jones',action="store_true",default=False,
                  help='Just show Jones matrices')
   opt.add_option('-l','--linear',action="store_true",default=False,
                  help='Keep in linear basis')
   options, arguments = opt.parse_args()
   correctMSforDipole(options.msfile)
