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

def getBandFilter(obsID):
    filterbandStr=getParam(obsID, "Observation.bandFilter", "string")
    bandlabel,filter_low_freq_str,filter_hi_freq_str=filterbandStr.split("_")
    return bandlabel,float(filter_low_freq_str),float(filter_hi_freq_str)

def getSubbandFreqs(obsID):
    sbwidth=getParam(obsID,"Observation.subbandWidth","float")
    sblst=getSubbandList(obsID)
    bandlabel,filter_low_freq,filter_hi_freq=getBandFilter(obsID)
    CLOCKRES=sbwidth/1000.0
#/* VLAD: 17.10.2012 Should work for all possible filters and antennas */
#float extra_half, lower_edge;
    if CLOCKRES == 0.1953125:
       #// 200 MHz clock
       if filter_low_freq >= 200:
          lower_edge = 200.0
       elif filter_low_freq < 200 and filter_low_freq >= 100:
          lower_edge = 100.0
       else:
          lower_edge = 0.0
    else : #// 160 MHz clock
       if filter_low_freq >= 160:
          lower_edge = 160.0
       elif filter_low_freq < 160 and filter_low_freq >= 80:
          lower_edge = 80.0
       else:
          lower_edge = 0.0
#// this line takes care if we use 2nd PPF or not as in the case when we bypass 2nd PPF (1chan/sub)
#// we do not need to subtract half of the channel width for frequency calculation
    #if NCHANNELS > 1:
    #   extra_half = 0.5
    #else:
    #   extra_half = 0.0

    freqs=[]
    for sbNr in sblst:
    #  the_freq=lower_edge+(SUBBLIST[isub]+float(ichan)/float(NCHANNELS)-extra_half)*CLOCKRES
      the_freq=lower_edge+(sbNr)*CLOCKRES
      #freqs.append(sbNr*sbwidth*1000.0) #in Hz
      freqs.append(the_freq*1E6) #in Hz
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
