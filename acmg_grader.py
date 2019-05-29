import sys

class Grader:
    __drugName=""
    __graderMap={{}}

    def __init__(self,drugName):
        self.__drugName = drugName

    def getDrugName(self):
        return self.__drugName

    def addGraderData(self,input_geneName,input_hgvs,input_grade):
        inner_graderMap={}


        if input_geneName in self.__graderMap:
            inner_graderMap = self.__graderMap.get(input_geneName)
            if input_hgvs in inner_graderMap:
            else:
                inner_graderMap


    def grading(self,input_geneName,input_hgvs):
