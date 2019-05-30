import sys

class Drug:

    __drugName = ""

    def __init__(self, name):
        self.__drugName = name

    def getDrugName(self):
        return self.__drugName

    def addGradeRule(self, gene, hgvs, grade):
        return None
