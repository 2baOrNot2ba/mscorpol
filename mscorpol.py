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

__version__="0.0"

def correctMSforDipole(msfile):
  me=measures()
  mstab = pt.table(msfile,readonly=False,ack=True)
  ta = pt.table(mstab.getkeyword('ANTENNA'),ack=False)
  NrOfAnts=ta.nrows()
  pos = ta.getcol('POSITION')
  tl = pt.table(mstab.getkeyword('LOFAR_ANTENNA_FIELD'),ack=False)
  tf = pt.table(mstab.getkeyword('FIELD'),ack=False)
  srcDir = np.squeeze(tf.getcol('PHASE_DIR'))
  RA = quantity(srcDir[0],'rad'); dec = quantity(srcDir[1],'rad')
  srcDirection = me.direction('J2000',RA,dec)
  for baselinehalf in ['ANTENNA1', 'ANTENNA2']:
    if baselinehalf == 'ANTENNA1':
       dataMultIdx=0
    else:
       dataMultIdx=1
    for tant in mstab.iter([baselinehalf]):
      antNr=int(tant.getcol(baselinehalf)[0])
      print baselinehalf,"=",antNr,'/',NrOfAnts
      x = quantity(pos[antNr,0],'m');
      y = quantity(pos[antNr,1],'m')
      z = quantity(pos[antNr,2],'m')
      stnPos = me.position('ITRF',x,y,z)
      timevals = tant.getcol('TIME')
      tt = quantity(timevals,'s')
      stnRot=tl.getcol('COORDINATE_AXES')[antNr]
      JI=getDipJones(tt,stnPos,stnRot,srcDirection,
                     doCirc=not(options.linear),doInvJ=True,showJones=options.jones)
      #Visibility data is such that 
      #  DATA(ANTENNA1,ANTENNA2)=V(ANTENNA1)^* x V(ANTENNA2),
      #where V() is the raw voltage of the corresponding antenna.
      if baselinehalf == 'ANTENNA1':
         JI=np.conj(JI)

      if not options.jones:
         #Start processing DATA
         data=tant.getcol('DATA')    
         dataCohXY=np.array([[data[:,:,0],data[:,:,1]],[data[:,:,2],data[:,:,3]]])
         dataCor=np.zeros(dataCohXY.shape,dtype=np.complex)
         #Multiply Inverse Jones with visibility data.
         for JIt in range(JI.shape[0]):
             dataCor[:,:,JIt,:]=np.tensordot(JI[JIt,:,:],dataCohXY[:,:,JIt,:],
                         axes=[[1],[dataMultIdx]])

         dataCorOut=np.transpose(
          [dataCor[0,0,:,:],dataCor[0,1,:,:],dataCor[1,0,:,:],dataCor[1,1,:,:]],
                              (1,2,0))
         #N.B. the output matrix is transposed incorrectly. (Fix it!)
         tant.putcol('DATA',dataCorOut)
  if not options.jones: updateMSmetadata(msfile)

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
