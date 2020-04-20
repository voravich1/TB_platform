import sys
from drug import Drug

#### Grade idicator will have 4 indicator
# 1.high  2.moderate  3.low  4.indeterminate
#########

class svcluster:

    def __init__(self):
        self.__cluster_name = ""


    def add_svevent(self,filename,record):
        return None

class svevent:

    def __init__(self,filename: str, record):
        samplename = filename.split('_')[0].split('.')[0]
        self.__sample__name = samplename
        self.file_path = filename
        self.record = record

    def getRecord(self):
        return self.record

    def getSampleName(self):
        return self.__sample__name

    def getFilePath(self):
        return self.file_path