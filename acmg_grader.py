import sys
from drug import Drug

class Grader:

    def __init__(self,gradeFile):
        self.__drugName = ""
        self.__drug_map = {}
        self.__gradeFile = gradeFile
        self.createGraderBackBone()

    def getDrugName(self):
        return self.__drugName

    #def addGraderData(self,input_geneName,input_hgvs,input_grade):
        #inner_graderMap={}


        #if input_geneName in self.__graderMap:
            #inner_graderMap = self.__graderMap.get(input_geneName)
            #if input_hgvs in inner_graderMap:
            #else:
                #inner_graderMap

    def createGraderBackBone(self):
        firstFlag = True
        with open(self.__gradeFile) as file:
            for line in file:
                line_noEnd = line.strip()
                if firstFlag==True:
                    firstFlag=False
                    continue

                info = line_noEnd.split("\t")
                drug_name = info[0]
                gene_name = info[1]

                # check existing drug object
                if drug_name in self.__drug_map:
                    drug_object = self.__drug_map[drug_name]
                else:
                    drug_object = Drug(drug_name)

                # begin add info to drug object
                if info[2] != ".":
                    hgvs_list = info[2].split(",")
                    for hgvs in hgvs_list:
                        drug_object.addGradeRule(gene_name, hgvs, "high")

                if info[3] != ".":
                    hgvs_list = info[3].split(",")
                    for hgvs in hgvs_list:
                        drug_object.addGradeRule(gene_name, hgvs, "moderate")

                if info[4] != ".":
                    hgvs_list = info[4].split(",")
                    for hgvs in hgvs_list:
                        drug_object.addGradeRule(gene_name, hgvs, "minimal")

                self.__drug_map[drug_name] = drug_object

    def getGrade(self,drug,gene,hgvs):

        if drug in self.__drug_map:
            drug_object = self.__drug_map.get(drug)
            grade = drug_object.getGrade(gene,hgvs)
            return grade
        else:
            return "No grade"