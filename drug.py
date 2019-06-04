import sys

class Drug:

    def __init__(self, name):
        self.__geneGradeMap = {}
        self.__drugName = name

    def getDrugName(self):
        return self.__drugName

    def addGradeRule(self, gene, hgvs, grade):

        if gene in self.__geneGradeMap:
            hgvsGradeMap = self.__geneGradeMap.get(gene)
        else:
            #geneVar = GeneVariant(gene)
            hgvsGradeMap = {}

        if hgvs in hgvsGradeMap:
            print("It should not happen. This hgvs and grade is not unique.\n")
            print("hgvs is = %s \t grade is = %s\n" % (hgvs, grade))
        else:
            hgvsGradeMap[hgvs] = grade

        self.__geneGradeMap[gene] = hgvsGradeMap

    def getGrade(self,gene,hgvs):
        hgvsGradeMap = {}
        if gene in self.__geneGradeMap:
            hgvsGradeMap = self.__geneGradeMap.get(gene)

            if hgvs in hgvsGradeMap:
                grade = hgvsGradeMap.get(hgvs)
                return grade
            else:
                return "No grade"
        else:
            return "No grade"

