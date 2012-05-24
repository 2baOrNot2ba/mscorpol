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
            else:
               paramValue=valStrFull
        line=f.readline()
    #print paramValue
    return paramValue
