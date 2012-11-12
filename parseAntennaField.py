AntennaFieldDirectory='/opt/lofar/etc/StaticMetaData/' #'/opt/cep/lofar/share/AntennaFields/'
COMMENT_CHAR = '#'

def parseAntennaField(stationName,AFfileNameType=2):

    if AFfileNameType==2:
       basename=stationName+'-'+'AntennaField'+'.conf'
    else:
       basename='AntennaField'+stationName+'.conf'
    filename=AntennaFieldDirectory+'/'+basename
    AntFldData={'LBA': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA0': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA1': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]}
}
    f = open(filename)
    line=f.readline()

    while line:

        if COMMENT_CHAR in line:
           line, comment = line.split(COMMENT_CHAR, 1)
        if "HBA0" in line:
              AntBand="HBA0"
        elif "HBA1" in line:
              AntBand="HBA1"
        elif "HBA" in line:
              AntBand="HBA"
        elif "LBA" in line:
              AntBand="LBA"
        else:
	      line=f.readline()
	      continue
        where, rest = line.split(AntBand, 1)
  	where=where.strip()
        if where == '':
	      #Read absolute position of station origin
              line=f.readline()
              elementposLine=line.split()[2:5]
	      position=[float(v) for v in elementposLine]
              AntFldData[AntBand]['POSITION']=position
              if AntBand!="HBA0"  and AntBand!="HBA1":
	      #Read relative position of each element
                line=f.readline()
                dimstr,rest=line.split('[',1)
	        shp=dimstr.split('x')
                for elementNr in range(0,int(shp[0])):
                  line=f.readline()
                  vals=line.split()
		  xpos=[float(v) for v in vals[0:3]]
		  ypos=[float(v) for v in vals[3:6]]
		  AntFldData[AntBand]['REL_POS'].append([xpos,ypos])
	        #Read ending ']' line
                line=f.readline()
        elif where=='NORMAL_VECTOR':
              line=f.readline()
              elementposLine=line.split()[2:5]
	      nrmv=[float(v) for v in elementposLine]
              AntFldData[AntBand][where]=nrmv
        elif where=='ROTATION_MATRIX':
              line=f.readline()
              dimstr,rest=line.split('[',1)
	      shp=dimstr.split('x')
	      for xyz in range(3):
	          line=f.readline()
                  rowstr=line.split()
		  row=[float(v) for v in rowstr]
		  AntFldData[AntBand][where].append(row)
	      #Read ending ']' line
              line=f.readline()
        line=f.readline()
    return AntFldData

def getArrayElemRelPos(stnName,rcumode=5):
   AntFld=parseAntennaField(stnName)
   stnLoc=stnName[0:2]
   if rcumode==3:
      AntBand='LBA'
   elif rcumode==5:
      if stnLoc=='CS' or stnLoc=='RS':
         AntBand='HBA0'
      else:
         AntBand='HBA'
   #stnPos=np.matrix(AntFld[AntBand]['POSITION']).T
   #stnRot=np.matrix(AntFld[AntBand]['ROTATION_MATRIX'])
   stnRelPos=np.matrix(AntFld[AntBand]['REL_POS'])
   return stnRelPos

if __name__ == '__main__':
  AFD=parseAntennaField('UK608')
  print AFD['HBA']['POSITION']
