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
                drug_name = info[0].lower()
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

    def getGrade(self, drug, gene, hgvs):

        if drug in self.__drug_map:
            drug_object = self.__drug_map.get(drug)
            grade = drug_object.getGrade(gene, hgvs)
            return grade
        else:
            return "No grade"

    def getGradeSpecialCase(self, drug, gene, hgvs, totalGeneHgvs: dict):

        if drug in self.__drug_map:
            drug_object = self.__drug_map.get(drug)
            grade = drug_object.getGrade(gene,hgvs)
            grade = self.checkSpecialCase(drug, gene, hgvs, totalGeneHgvs, grade)
            return grade
        else:
            return "No grade"

    def checkSpecialCase(self, drug, gene, hgvs, totalGeneHgvs: dict, grade):
        '''This fuction is a fixed fuction that suitable only for TB grading
            If we have any update on grading criteria, this function should be re write
            In order to match with the updated criteria

            For now, grading criteria have two special case
            '''

        if drug == "isoniazid":
            if gene == "katG" and hgvs == "p.Ser315=":
                if "inhA" in totalGeneHgvs :
                    if "p.Ser315=" in totalGeneHgvs.get("inhA"):
                        return "high"

            elif gene == "inhA" and hgvs == "p.Ser315=":
                if "katG" in totalGeneHgvs :
                    if "p.Ser315=" in totalGeneHgvs.get("katG"):
                        return "high"

        elif drug == "rifampicin":
            if gene == "rpoB" and hgvs == "p.Asp516Gly":
                if "rpoB" in totalGeneHgvs :
                    if "p.Leu533Pro" in totalGeneHgvs.get("rpoB"):
                        return "high"

            elif gene == "rpoB" and hgvs == "p.Leu533Pro":
                if "rpoB" in totalGeneHgvs :
                    if "p.Asp516Gly" in totalGeneHgvs.get("rpoB"):
                        return "high"
                    else:
                        return "moderate"

        return grade
