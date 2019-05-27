#! /Users/worawich/miniconda3/envs/tbprofiler/bin/python

import sys

"""
    This program will grade drug resistance follow the acmg grading criteria
    first input is acmg grading database file
    second input is tbprofiler result txt format
"""
acmg_grading_db = sys.argv[1]




file = open(gene_file,'r')

for line in file:
    line_noend = line.strip()
    info = line_noend.split("\t")




def convertOne2ThreeAmino(oneLetter):

    if(oneLetter == "A"):
        threeLetter="Ala"
        return threeLetter
    elif(oneLetter == "R"):
        threeLetter = "Arg"
        return threeLetter
    elif (oneLetter == "N"):
        threeLetter = "Asn"
        return threeLetter
    elif (oneLetter == "D"):
        threeLetter = "Asp"
        return threeLetter
    elif (oneLetter == "C"):
        threeLetter = "Cys"
        return threeLetter
    elif (oneLetter == "E"):
        threeLetter = "Glu"
        return threeLetter
    elif (oneLetter == "Q"):
        threeLetter = "Gln"
        return threeLetter
    elif (oneLetter == "G"):
        threeLetter = "Gly"
        return threeLetter
    elif (oneLetter == "H"):
        threeLetter = "His"
        return threeLetter
    elif (oneLetter == "I"):
        threeLetter = "Ile"
        return threeLetter
    elif (oneLetter == "L"):
        threeLetter = "Leu"
        return threeLetter
    elif (oneLetter == "K"):
        threeLetter = "Lys"
        return threeLetter
    elif (oneLetter == "M"):
        threeLetter = "Met"
        return threeLetter
    elif (oneLetter == "F"):
        threeLetter = "Phe"
        return threeLetter
    elif (oneLetter == "P"):
        threeLetter = "Pro"
        return threeLetter
    elif (oneLetter == "U"):
        threeLetter = "Gip"
        return threeLetter
    elif (oneLetter == "S"):
        threeLetter = "Ser"
        return threeLetter
    elif (oneLetter == "T"):
        threeLetter = "Thr"
        return threeLetter
    elif (oneLetter == "W"):
        threeLetter = "Trp"
        return threeLetter
    elif (oneLetter == "Y"):
        threeLetter = "Tyr"
        return threeLetter
    elif (oneLetter == "V"):
        threeLetter = "Val"
        return threeLetter




