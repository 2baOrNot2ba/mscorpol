parsetDirectory='/globalhome/lofarsystem/log/' #Directory on CEP1 where parsets live

#parsetData={'stations': [] }

def getStnNames(obsID):
    return getParam(obsID,"Observation.existingStations","list")

def getStartTime(obsID):
    return getParam(obsID,"Observation.startTime","epoch")

def getDuration(obsID):
    return getParam(obsID,"Observation.Beam[0].duration","float")

def getStep(obsID):
    return getParam(obsID,"OLAP.Correlator.integrationTime","float")

def getRA(obsID):
    return getParam(obsID,"Observation.Beam[0].angle1","float")

def getDEC(obsID):
    return getParam(obsID,"Observation.Beam[0].angle2","float")

def getSubbandList(obsID):
    sblstGaps=getParam(obsID,"Observation.subbandList","list")
    sblst=[]
    for sblstElem in sblstGaps :
        sbStart,sbStop=sblstElem.split('..',2)
        sblst.extend(range(int(sbStart),int(sbStop)+1))
    return sblst

def getSubbandFreqs(obsID):
    sbwidth=getParam(obsID,"Observation.subbandWidth","float")
    sblst=getSubbandList(obsID)
    freqs=[]
    for sbNr in sblst:
      freqs.append(sbNr*sbwidth*1000.0) #in Hz
    return freqs

def getParam(obsID,paramName,dataType):

    basename=obsID
    filename=parsetDirectory+'/'+basename+'/'+basename+'.parset'
    f = open(filename)
    line=f.readline()

    while line:

        if paramName in line:
            param, valStrFull = line.split('=', 1)
            valStrFull=valStrFull[1:-1]
            if dataType == "list":
               valStrip=valStrFull[1:-1]
               #parsetData['stations']=valStrip.split(',')
               paramValue=valStrip.split(',')
            elif dataType == 'epoch':
               valStrip=valStrFull[1:-1]
               paramValue=valStrip
            elif dataType == 'float':
               valStrip=valStrFull
               paramValue=float(valStrip)
            else:
               paramValue=valStrFull
        line=f.readline()
    #print paramValue
    return paramValue

if __name__ == "__main__":
   sbf=getSubbandFreqs('L29053')
   print sbf
