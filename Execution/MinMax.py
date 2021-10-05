
#········································································································································#
#········································································································································#
#                                                 How to run                                                                             #
#                             python3 MinMax.py --inputFolderPath ../Preparation/descargas --outputPath .                                #
#········································································································································#
#········································································································································#


#This script was written in Python 3 and utilizes the statistics, random, pandas, and matplotlib.pyplot modules
#Imports happen throughout the script to minimize long pause times and to solve errors that occur when imports happen at the beginning

import statistics as s
import pandas as pd #from line 411
import random #from 355 This import is included here to lower the initial wait time for the script
import os
import argparse

#############################################################################################

__author__ = ''

# Parameters:
# 1) --inputFolderPath Path para leer la carpeta con archivos a evaluar
# 2) --outputPath Path para colocar la carpeta con todas las figuras generadas.

# Ejemplos de ejecucion
# python MinMax_modified.py --inputFolderPath .\descargas --outputPath .

# python MinMax_modified.py --inputFolderPath ..\descargas --outputPath .

if __name__ == "__main__":
    # Parameter definition
    parser = argparse.ArgumentParser(description='Calcular y visualizar valores de MinMax')
    parser.add_argument("--inputFolderPath", dest="inputFolderPath",
                      help="Path to read input folder", metavar="PATH")
    parser.add_argument("--outputPath", dest="outputPath",
                      help="Path to place output folder", metavar="PATH")
    args = parser.parse_args()


    # Printing parameter values
    print('-------------------------------- PARAMETROS --------------------------------')
    print("Path para leer el archivo de entrada: " + str(args.inputFolderPath))
    print("Path para el folder de figuras: " + str(args.outputPath))


    ##############################################################################################

    freqDict = dict() #dictionary mapping codons to their frequencies
    mapDict = dict() #dictionary mapping codons to amino acid
    aaFreqDict = dict() #dictionary mapping each amino acid to a list of the frequencies of possible codons for that amino acid
    aaMapDict = dict() #dictionary from amino acid to list of codons with frequencies for it (for RRTs)



    mapDict = {'TCA': 'S', 'AAT': 'N', 'TGG': 'W', 'GAT': 'D', 'GAA': 'E', 'TTC': 'F', 'CCG': 'P',
            'ACT': 'T', 'GGG': 'G', 'ACG': 'T', 'AGA': 'R', 'TTG': 'L', 'GTC': 'V', 'GCA': 'A',
            'TGA': '*', 'CGT': 'R', 'CAC': 'H', 'CTC': 'L', 'CGA': 'R', 'GCT': 'A', 'ATC': 'I',
            'ATA': 'I', 'TTT': 'F', 'TAA': '*', 'GTG': 'V', 'GCC': 'A', 'GAG': 'E', 'CAT': 'H',
            'AAG': 'K', 'AAA': 'K', 'GCG': 'A', 'TCC': 'S', 'GGC': 'G', 'TCT': 'S', 'CCT': 'P',
            'GTA': 'V', 'AGG': 'R', 'CCA': 'P', 'TAT': 'Y', 'ACC': 'T', 'TCG': 'S', 'ATG': 'M',
            'TTA': 'L', 'TGC': 'C', 'GTT': 'V', 'CTT': 'L', 'CAG': 'Q', 'CCC': 'P', 'ATT': 'I',
            'ACA': 'T', 'AAC': 'N', 'GGT': 'G', 'AGC': 'S', 'CGG': 'R', 'TAG': '*', 'CGC': 'R',
            'AGT': 'S', 'CTA': 'L', 'CAA': 'Q', 'CTG': 'L', 'GGA': 'G', 'TGT': 'C', 'TAC': 'Y',
            'GAC': 'D'}

    aaDict = {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'],
            'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
            'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'],
            'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'],
            'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'],
            'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'],
            'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}


    #The below dictionary contains the usage files for a handful of selected species. The files are taken from HIVE-CUT
    speciesDict = {
            'Escherichia_coli' : {'TTT': 22.38, 'TCT': 8.61, 'TAT': 16.36, 'TGT': 5.19, 'TTC': 16.21, #1
                                'TCC': 8.81, 'TAC': 12.15, 'TGC': 6.34, 'TTA': 13.83, 'TCA': 7.57,
                                'TAA': 2.03, 'TGA': 1.04, 'TTG': 13.37, 'TCG': 8.79, 'TAG': 0.25,
                                'TGG': 15.21, 'CTT': 11.44, 'CCT': 7.22, 'CAT': 12.84, 'CGT': 20.7,
                                'CTC': 10.92, 'CCC': 5.56, 'CAC': 9.44, 'CGC': 21.48, 'CTA': 3.93,
                                'CCA': 8.44, 'CAA': 15.1, 'CGA': 3.67, 'CTG': 52.1, 'CCG': 22.65,
                                'CAG': 29.21, 'CGG': 5.72, 'ATT': 30.21, 'ACT': 9.02, 'AAT': 18.26,
                                'AGT': 9.08, 'ATC': 24.6, 'ACC': 22.88, 'AAC': 21.47, 'AGC': 15.89,
                                'ATA': 4.88, 'ACA': 7.63, 'AAA': 33.94, 'AGA': 2.43, 'ATG': 27.59,
                                'ACG': 14.47, 'AAG': 10.7, 'AGG': 1.48, 'GTT': 18.39, 'GCT': 15.54,
                                'GAT': 32.43, 'GGT': 24.45, 'GTC': 15.07, 'GCC': 25.45, 'GAC': 19.14,
                                'GGC': 28.65, 'GTA': 10.97, 'GCA': 20.61, 'GAA': 39.55, 'GGA': 8.44,
                                'GTG': 25.9, 'GCG': 32.79, 'GAG': 18.24, 'GGG': 11.29},


            'Homo_sapien' : {'TTT': 21.42, 'TCT': 16.96, 'TAT': 17.11, 'TGT': 10.99, 'TTC': 23.05,  #2
                            'TCC': 10.61, 'TAC': 13.48, 'TGC': 8.76, 'TTA': 9.31, 'TCA': 21.14,
                            'TAA': 1.07, 'TGA': 0.72, 'TTG': 19.58, 'TCG': 12.81, 'TAG': 0.39,
                            'TGG': 10.68, 'CTT': 20.92, 'CCT': 9.16, 'CAT': 14.3, 'CGT': 11.38,
                            'CTC': 14.54, 'CCC': 4.46, 'CAC': 9.07, 'CGC': 5.03, 'CTA': 7.69,
                            'CCA': 27.25, 'CAA': 27.88, 'CGA': 12.5, 'CTG': 12.05, 'CCG': 10.28,
                            'CAG': 14.88, 'CGG': 4.84, 'ATT': 31.51, 'ACT': 19.3, 'AAT': 30.06,
                            'AGT': 12.28, 'ATC': 18.5, 'ACC': 10.32, 'AAC': 17.93, 'AGC': 8.32,
                            'ATA': 9.03, 'ACA': 20.55, 'AAA': 36.4, 'AGA': 15.25, 'ATG': 26.09,
                            'ACG': 9.22, 'AAG': 25.58, 'AGG': 3.8, 'GTT': 24.14, 'GCT': 22.94,
                            'GAT': 36.73, 'GGT': 11.12, 'GTC': 13.58, 'GCC': 12.74, 'GAC': 17.27,
                            'GGC': 6.78, 'GTA': 9.78, 'GCA': 20.35, 'GAA': 41.61, 'GGA': 32.03,
                            'GTG': 14.53, 'GCG': 8.56, 'GAG': 25.1, 'GGG': 4.32},


            'Mus_musculus' : {'TTT': 15.94, 'TCT': 17.39, 'TAT': 11.15, 'TGT': 10.68, 'TTC': 18.81, #3
                            'TCC': 18.32, 'TAC': 14.42, 'TGC': 10.95, 'TTA': 7.29, 'TCA': 13.31,
                            'TAA': 0.39, 'TGA': 0.76, 'TTG': 13.27, 'TCG': 4.29, 'TAG': 0.35,
                            'TGG': 11.44, 'CTT': 13.45, 'CCT': 20.06, 'CAT': 11.23, 'CGT': 4.64,
                            'CTC': 18.83, 'CCC': 18.34, 'CAC': 15.23, 'CGC': 8.6, 'CTA': 8.04,
                            'CCA': 19.05, 'CAA': 13.11, 'CGA': 6.9, 'CTG': 37.31, 'CCG': 6.11,
                            'CAG': 36.71, 'CGG': 10.46, 'ATT': 14.66, 'ACT': 13.92, 'AAT': 15.9,
                            'AGT': 14.11, 'ATC': 20.33, 'ACC': 18.16, 'AAC': 19.75, 'AGC': 20.77,
                            'ATA': 7.28, 'ACA': 16.67, 'AAA': 23.84, 'AGA': 13.03, 'ATG': 21.7,
                            'ACG': 5.73, 'AAG': 33.95, 'AGG': 12.94, 'GTT': 10.81, 'GCT': 20.19,
                            'GAT': 22.33, 'GGT': 11.09, 'GTC': 14.52, 'GCC': 25.16, 'GAC': 26.3,
                            'GGC': 19.81, 'GTA': 7.48, 'GCA': 16.8, 'GAA': 30.33, 'GGA': 16.77,
                            'GTG': 26.58, 'GCG': 5.86, 'GAG': 41.48, 'GGG': 14.91},

            'Caenorhabditis_elegans' : {'TTT': 17.06, 'TCT': 16.58, 'TAT': 12.04, 'TGT': 10.54, 'TTC': 17.87, #4
                                        'TCC': 17.44, 'TAC': 13.7, 'TGC': 11.15, 'TTA': 8.55, 'TCA': 13.89,
                                        'TAA': 0.46, 'TGA': 0.83, 'TTG': 13.3, 'TCG': 4.18, 'TAG': 0.36,
                                        'TGG': 11.77, 'CTT': 13.95, 'CCT': 18.88, 'CAT': 11.74, 'CGT': 4.54,
                                        'CTC': 18.06, 'CCC': 19.19, 'CAC': 14.76, 'CGC': 9.06, 'CTA': 7.39,
                                        'CCA': 18.45, 'CAA': 13.83, 'CGA': 6.36, 'CTG': 36.75, 'CCG': 6.36,
                                        'CAG': 35.3, 'CGG': 10.88, 'ATT': 16.36, 'ACT': 14.12, 'AAT': 18.16,
                                        'AGT': 13.72, 'ATC': 18.97, 'ACC': 17.95, 'AAC': 18.36, 'AGC': 19.74,
                                        'ATA': 7.98, 'ACA': 16.33, 'AAA': 27.15, 'AGA': 13.09, 'ATG': 21.4,
                                        'ACG': 5.72, 'AAG': 31.89, 'AGG': 12.15, 'GTT': 11.59, 'GCT': 18.77,
                                        'GAT': 23.68, 'GGT': 10.75, 'GTC': 13.58, 'GCC': 26.18, 'GAC': 24.49,
                                        'GGC': 20.23, 'GTA': 7.56, 'GCA': 16.89, 'GAA': 33.04, 'GGA': 17.02,
                                        'GTG': 26.24, 'GCG': 6.26, 'GAG': 39.88, 'GGG': 15.53},


        'Saccharomyces_cerevisiae' : {'TTT': 26.18, 'TCT': 23.35, 'TAT': 19.05, 'TGT': 7.82, 'TTC': 17.88, #5
                                    'TCC': 14.07, 'TAC': 14.6, 'TGC': 4.75, 'TTA': 26.33, 'TCA': 19.05,
                                    'TAA': 0.95, 'TGA': 0.6, 'TTG': 26.5, 'TCG': 8.71, 'TAG': 0.46,
                                    'TGG': 10.35, 'CTT': 12.27, 'CCT': 13.57, 'CAT': 13.89, 'CGT': 6.26,
                                    'CTC': 5.52, 'CCC': 6.91, 'CAC': 7.74, 'CGC': 2.63, 'CTA': 13.52,
                                    'CCA': 17.81, 'CAA': 27.1, 'CGA': 3.1, 'CTG': 10.65, 'CCG': 5.42,
                                    'CAG': 12.42, 'CGG': 1.82, 'ATT': 30.1, 'ACT': 20.24, 'AAT': 36.61,
                                    'AGT': 14.6, 'ATC': 16.99, 'ACC': 12.48, 'AAC': 24.8, 'AGC': 9.96,
                                    'ATA': 18.29, 'ACA': 18.18, 'AAA': 42.83, 'AGA': 21.05, 'ATG': 20.68,
                                    'ACG': 8.15, 'AAG': 30.52, 'AGG': 9.45, 'GTT': 21.47, 'GCT': 20.28,
                                    'GAT': 38.09, 'GGT': 22.59, 'GTC': 11.23, 'GCC': 12.14, 'GAC': 20.39,
                                    'GGC': 9.78, 'GTA': 12.07, 'GCA': 16.26, 'GAA': 45.81, 'GGA': 11.19,
                                    'GTG': 10.72, 'GCG': 6.17, 'GAG': 19.55, 'GGG': 6.06},

        'Homo_sapien_last' : {'TTT': 17.22, 'TCT': 16.96, 'TAT': 12.15, 'TGT': 10.47, 'TTC': 17.51, #6
                         'TCC': 17.29, 'TAC': 13.47, 'TGC': 10.83, 'TTA': 8.82, 'TCA': 14.19,
                         'TAA': 0.44, 'TGA': 0.80, 'TTG': 13.47, 'TCG': 4.05, 'TAG': 0.35,
                         'TGG': 11.66, 'CTT': 14.17, 'CCT': 19.15, 'CAT': 11.90, 'CGT': 4.56,
                         'CTC': 17.81, 'CCC': 18.86, 'CAC': 14.62, 'CGC': 8.71, 'CTA': 7.48,
                         'CCA': 18.79, 'CAA': 14.17, 'CGA': 6.42, 'CTG': 36.01, 'CCG': 6.15,
                         'CAG': 35.28, 'CGG': 10.62, 'ATT': 16.58, 'ACT': 14.30, 'AAT': 18.52,
                         'AGT': 14.06, 'ATC': 18.68, 'ACC': 17.77, 'AAC': 18.27, 'AGC': 19.66,
                         'ATA': 8.16, 'ACA': 16.57, 'AAA': 27.77, 'AGA': 13.40, 'ATG': 21.48,
                         'ACG': 5.59, 'AAG': 31.79, 'AGG': 12.17, 'GTT': 11.80, 'GCT': 18.89,
                         'GAT': 24.09, 'GGT': 10.82, 'GTC': 13.44, 'GCC': 25.64, 'GAC': 24.25,
                         'GGC': 19.73, 'GTA': 7.71, 'GCA': 17.09, 'GAA': 33.85, 'GGA': 17.19,
                         'GTG': 25.75, 'GCG': 5.94, 'GAG': 39.40, 'GGG': 15.26},

        'Danio_rerio':{'TTT': 16.38    , 'TCT': 19.06    , 'TAT': 11.30    , 'TGT': 11.21    , #7
                        'TTC': 18.00    , 'TCC': 15.84    , 'TAC': 14.73    , 'TGC': 10.00    ,
                        'TTA': 7.19    , 'TCA': 16.25    , 'TAA': 0.50    , 'TGA': 0.80    ,
                        'TTG': 12.23    , 'TCG': 5.20    , 'TAG': 0.28    , 'TGG': 10.37    ,
                        'CTT': 13.23    , 'CCT': 18.74    , 'CAT': 11.51    , 'CGT': 6.77    ,
                        'CTC': 16.48    , 'CCC': 13.00    , 'CAC': 14.86    , 'CGC': 9.00    ,
                        'CTA': 6.94    , 'CCA': 18.13    , 'CAA': 14.41    , 'CGA': 6.78    ,
                        'CTG': 35.90    , 'CCG': 7.25    , 'CAG': 34.71    , 'CGG': 6.38    ,
                        'ATT': 15.68    , 'ACT': 15.80    , 'AAT': 16.89    , 'AGT': 15.40    ,
                        'ATC': 21.62    , 'ACC': 15.96    , 'AAC': 22.84    , 'AGC': 19.17    ,
                        'ATA': 7.87    , 'ACA': 19.21    , 'AAA': 30.88    , 'AGA': 15.63    ,
                        'ATG': 23.59    , 'ACG': 6.79    , 'AAG': 29.86    , 'AGG': 10.67    ,
                        'GTT': 13.87    , 'GCT': 20.71    , 'GAT': 25.97    , 'GGT': 13.20    ,
                        'GTC': 14.08    , 'GCC': 18.23    , 'GAC': 26.91    , 'GGC': 15.73    ,
                        'GTA': 7.11    , 'GCA': 17.16    , 'GAA': 27.85    , 'GGA': 20.85    ,
                        'GTG': 26.12    , 'GCG': 6.86    , 'GAG': 44.24    , 'GGG': 9.85},

        'Gallus_gallus': {'TTT':16.78 ,'TCT':17.35 ,'TAT':11.72 ,'TGT':9.71 , #8
                        'TTC':17.05 ,'TCC':16.01 ,'TAC':14.77 ,'TGC':11.89 ,
                        'TTA':8.42 ,'TCA':14.50 ,'TAA':0.45 ,'TGA':0.73 ,
                        'TTG':13.64 ,'TCG':4.65 ,'TAG':0.33 ,'TGG':11.19 ,
                        'CTT':13.89 ,'CCT':17.82 ,'CAT':11.29 ,'CGT':5.47 ,
                        'CTC':15.51 ,'CCC':15.93 ,'CAC':14.10 ,'CGC':8.89 ,
                        'CTA':6.80 ,'CCA':18.56 ,'CAA':14.55 ,'CGA':5.63 ,
                        'CTG':35.12 ,'CCG':7.06 ,'CAG':34.46 ,'CGG':9.03 ,
                        'ATT':16.63 ,'ACT':14.62 ,'AAT':18.49 ,'AGT':14.02 ,
                        'ATC':18.18 ,'ACC':14.68 ,'AAC':19.93 ,'AGC':21.06 ,
                        'ATA':9.42 ,'ACA':17.93 ,'AAA':29.21 ,'AGA':14.18 ,
                        'ATG':22.01 ,'ACG':6.88 ,'AAG':31.70 ,'AGG':12.55 ,
                        'GTT':13.94 ,'GCT':21.67 ,'GAT':26.79 ,'GGT':11.94 ,
                        'GTC':12.43 ,'GCC':19.75 ,'GAC':22.69 ,'GGC':17.41 ,
                        'GTA':8.53 ,'GCA':21.12 ,'GAA':34.36 ,'GGA':18.67 ,
                        'GTG':25.10 ,'GCG':7.21 ,'GAG':38.91 ,'GGG':14.69},

        'Xenopus_tropicalis': {'TTT':21.50 ,'TCT':18.31 ,'TAT':15.33 ,'TGT':10.52 , #9
                        'TTC':16.61 ,'TCC':14.95 ,'TAC':14.85 ,'TGC':10.96 ,
                        'TTA':10.43 ,'TCA':13.47 ,'TAA':0.94 ,'TGA':0.89 ,
                        'TTG':14.81 ,'TCG':4.03 ,'TAG':0.45 ,'TGG':11.44 ,
                        'CTT':17.41 ,'CCT':16.95 ,'CAT':12.78 ,'CGT':5.92 ,
                        'CTC':13.01 ,'CCC':12.42 ,'CAC':12.54 ,'CGC':7.11 ,
                        'CTA':9.63 ,'CCA':19.34 ,'CAA':16.30 ,'CGA':6.10 ,
                        'CTG':29.07 ,'CCG':5.12 ,'CAG':30.09 ,'CGG':7.04 ,
                        'ATT':20.70 ,'ACT':15.34 ,'AAT':21.93 ,'AGT':14.13 ,
                        'ATC':17.01 ,'ACC':14.16 ,'AAC':20.02 ,'AGC':16.67 ,
                        'ATA':12.11 ,'ACA':18.80 ,'AAA':33.06 ,'AGA':14.59 ,
                        'ATG':24.12 ,'ACG':4.80 ,'AAG':31.71 ,'AGG':11.96 ,
                        'GTT':16.25 ,'GCT':20.56 ,'GAT':29.80 ,'GGT':11.86 ,
                        'GTC':11.77 ,'GCC':17.94 ,'GAC':22.23 ,'GGC':14.43 ,
                        'GTA':10.89 ,'GCA':20.51 ,'GAA':36.03 ,'GGA':20.43 ,
                        'GTG':22.08 ,'GCG':5.09 ,'GAG':34.65 ,'GGG':14.06 },

        'Thermus_thermophilus':{'TTT':7.07 ,'TCT':1.30 ,'TAT':1.45 ,'TGT':0.28, #10
                        'TTC':30.36 ,'TCC':15.71 ,'TAC':26.71 ,'TGC':4.02 ,
                        'TTA':1.27,'TCA':0.64 ,'TAA':0.66 ,'TGA':1.64 ,
                        'TTG':10.11 ,'TCG':4.19 ,'TAG':1.27 ,'TGG':14.87 ,
                        'CTT':17.55 ,'CCT':5.46 ,'CAT':1.12 ,'CGT':1.84 ,
                        'CTC':70.39 ,'CCC':47.53 ,'CAC':17.46 ,'CGC':28.34 ,
                        'CTA':3.57 ,'CCA':1.97 ,'CAA':3.45 ,'CGA':1.47 ,
                        'CTG':40.61 ,'CCG':11.92 ,'CAG':21.07 ,'CGG':36.52 ,
                        'ATT':2.47 ,'ACT':0.68 ,'AAT':0.51 ,'AGT':0.56,
                        'ATC':22.02 ,'ACC':27.19 ,'AAC':14.79 ,'AGC':13.73 ,
                        'ATA':1.24 ,'ACA':0.68 ,'AAA':3.25 ,'AGA':1.20,
                        'ATG':14.37 ,'ACG':9.60 ,'AAG':32.23 ,'AGG':16.93 ,
                        'GTT':2.87 ,'GCT':3.54 ,'GAT':1.98 ,'GGT':2.89 ,
                        'GTC':26.18 ,'GCC':82.82 ,'GAC':33.09 ,'GGC':36.54 ,
                        'GTA':1.72 ,'GCA':2.17 ,'GAA':11.21 ,'GGA':6.40 ,
                        'GTG':49.10 ,'GCG':25.60 ,'GAG':72.62 ,'GGG':48.01},


        'Oryctolagus_cuniculus': {'TTT':15.98 ,'TCT':15.55 ,'TAT':10.53 ,'TGT':9.96 , #11
                        'TTC':18.79 ,'TCC':18.51 ,'TAC':14.68 ,'TGC':12.01 ,
                        'TTA':7.91 ,'TCA':12.21 ,'TAA':0.54 ,'TGA':1.25 ,
                        'TTG':12.95 ,'TCG':5.52 ,'TAG':0.46 ,'TGG':11.78 ,
                        'CTT':12.72 ,'CCT':18.21 ,'CAT':10.33 ,'CGT':4.32 ,
                        'CTC':18.63 ,'CCC':20.40 ,'CAC':15.77 ,'CGC':10.09 ,
                        'CTA':6.35 ,'CCA':17.05 ,'CAA':12.87 ,'CGA':6.33 ,
                        'CTG':38.91 ,'CCG':8.54 ,'CAG':36.33 ,'CGG':11.39 ,
                        'ATT':15.04 ,'ACT':12.89 ,'AAT':16.29 ,'AGT':13.10 ,
                        'ATC':19.62 ,'ACC':18.18 ,'AAC':19.15 ,'AGC':20.82 ,
                        'ATA':7.34 ,'ACA':14.70 ,'AAA':25.40 ,'AGA':13.02 ,
                        'ATG':20.93 ,'ACG':7.42 ,'AAG':32.30 ,'AGG':12.61 ,
                        'GTT':10.93 ,'GCT':18.35 ,'GAT':21.48 ,'GGT':9.98 ,
                        'GTC':14.16 ,'GCC':28.02 ,'GAC':26.27 ,'GGC':21.70 ,
                        'GTA':6.73 ,'GCA':16.15 ,'GAA':31.38 ,'GGA':16.88 ,
                        'GTG':27.57 ,'GCG':8.56 ,'GAG':40.36 ,'GGG':15.77},

        'Rattus_norvegicus':{'TTT':15.99 ,'TCT':17.10 ,'TAT':11.21 ,'TGT':11.02 , #12
                        'TTC':19.47 ,'TCC':18.20 ,'TAC':14.90 ,'TGC':11.46 ,
                        'TTA':7.28 ,'TCA':13.30 ,'TAA':0.57 ,'TGA':1.18 ,
                        'TTG':13.39 ,'TCG':4.28 ,'TAG':0.49 ,'TGG':11.87 ,
                        'CTT':13.54 ,'CCT':19.67 ,'CAT':11.16 ,'CGT':4.76 ,
                        'CTC':19.07 ,'CCC':18.17 ,'CAC':15.47 ,'CGC':8.57 ,
                        'CTA':8.15 ,'CCA':18.59 ,'CAA':12.96 ,'CGA':6.78 ,
                        'CTG':37.51 ,'CCG':6.17 ,'CAG':36.13 ,'CGG':10.30 ,
                        'ATT':14.67 ,'ACT':13.97 ,'AAT':15.50 ,'AGT':13.98 ,
                        'ATC':20.39 ,'ACC':18.21 ,'AAC':19.65 ,'AGC':20.25 ,
                        'ATA':7.39 ,'ACA':16.64 ,'AAA':23.75 ,'AGA':13.27 ,
                        'ATG':21.70 ,'ACG':5.77 ,'AAG':34.05 ,'AGG':13.08 ,
                        'GTT':10.86 ,'GCT':19.89 ,'GAT':21.52 ,'GGT':11.23 ,
                        'GTC':14.65 ,'GCC':25.26 ,'GAC':26.43 ,'GGC':19.70 ,
                        'GTA':7.58 ,'GCA':16.70 ,'GAA':29.67 ,'GGA':16.82 ,
                        'GTG':26.79 ,'GCG':6.11 ,'GAG':40.83 ,'GGG':14.98},

        'Vibrio_harveyi':{'TTT':22.57 ,'TCT':16.15 ,'TAT':11.87 ,'TGT':6.51 , #13
                        'TTC':19.54 ,'TCC':4.50 ,'TAC':19.17 ,'TGC':3.75 ,
                        'TTA':17.65 ,'TCA':12.42 ,'TAA':2.23 ,'TGA':0.56 ,
                        'TTG':21.64 ,'TCG':8.18 ,'TAG':0.62 ,'TGG':12.54 ,
                        'CTT':18.80 ,'CCT':11.83 ,'CAT':10.39 ,'CGT':19.17 ,
                        'CTC':8.96 ,'CCC':2.36 ,'CAC':11.60 ,'CGC':12.68 ,
                        'CTA':16.91 ,'CCA':16.87 ,'CAA':31.71 ,'CGA':5.93 ,
                        'CTG':17.19 ,'CCG':7.41 ,'CAG':13.23 ,'CGG':0.92 ,
                        'ATT':28.96 ,'ACT':14.47 ,'AAT':16.53 ,'AGT':11.13 ,
                        'ATC':27.90 ,'ACC':14.60 ,'AAC':26.15 ,'AGC':13.67 ,
                        'ATA':4.54 ,'ACA':12.71 ,'AAA':34.89 ,'AGA':3.81 ,
                        'ATG':26.93 ,'ACG':12.55 ,'AAG':18.75 ,'AGG':1.19 ,
                        'GTT':23.20 ,'GCT':21.58 ,'GAT':32.99 ,'GGT':33.93 ,
                        'GTC':11.88 ,'GCC':11.46 ,'GAC':22.65 ,'GGC':20.89 ,
                        'GTA':16.00 ,'GCA':28.13 ,'GAA':41.59 ,'GGA':6.78 ,
                        'GTG':20.76 ,'GCG':23.93 ,'GAG':23.68 ,'GGG':5.95},

        'Pseudomonas_aeruginosa':{'TTT':2.13 ,'TCT':1.19 ,'TAT':5.40 ,'TGT':1.15 , #14
                        'TTC':33.08 ,'TCC':12.09 ,'TAC':19.78 ,'TGC':9.15 ,
                        'TTA':0.43 ,'TCA':0.93 ,'TAA':0.33 ,'TGA':2.47 ,
                        'TTG':9.09 ,'TCG':13.21 ,'TAG':0.39 ,'TGG':15.06 ,
                        'CTT':3.66 ,'CCT':2.61 ,'CAT':6.48 ,'CGT':8.05 ,
                        'CTC':27.35 ,'CCC':13.20 ,'CAC':15.22 ,'CGC':48.25 ,
                        'CTA':1.69 ,'CCA':2.62 ,'CAA':6.68 ,'CGA':2.84 ,
                        'CTG':80.60 ,'CCG':32.55 ,'CAG':35.84 ,'CGG':14.24 ,
                        'ATT':3.30 ,'ACT':2.02 ,'AAT':3.99 ,'AGT':2.79 ,
                        'ATC':37.28 ,'ACC':32.18 ,'AAC':22.34 ,'AGC':25.48 ,
                        'ATA':1.13 ,'ACA':1.17 ,'AAA':4.11 ,'AGA':0.71 ,
                        'ATG':20.32 ,'ACG':6.98 ,'AAG':24.98 ,'AGG':2.35 ,
                        'GTT':3.20 ,'GCT':5.56 ,'GAT':11.08 ,'GGT':8.49 ,
                        'GTC':28.37 ,'GCC':66.17 ,'GAC':41.90 ,'GGC':60.20 ,
                        'GTA':4.11 ,'GCA':5.68 ,'GAA':23.35 ,'GGA':4.51 ,
                        'GTG':32.74 ,'GCG':38.45 ,'GAG':37.12 ,'GGG':10.16},

        'Salmonella_enterica':{'TTT':23.35 ,'TCT':7.52 ,'TAT':17.21 ,'TGT':5.05 , #15
                        'TTC':15.22 ,'TCC':10.04 ,'TAC':11.47 ,'TGC':6.82 ,
                        'TTA':13.57 ,'TCA':6.56 ,'TAA':2.10 ,'TGA':1.29 ,
                        'TTG':12.71 ,'TCG':9.61 ,'TAG':0.40 ,'TGG':15.31 ,
                        'CTT':12.01 ,'CCT':7.37 ,'CAT':13.38 ,'CGT':18.54 ,
                        'CTC':10.50 ,'CCC':6.99 ,'CAC':9.51 ,'CGC':23.09 ,
                        'CTA':5.13 ,'CCA':6.04 ,'CAA':12.96 ,'CGA':3.85 ,
                        'CTG':52.22 ,'CCG':24.32 ,'CAG':30.38 ,'CGG':7.02 ,
                        'ATT':29.34 ,'ACT':6.99 ,'AAT':18.06 ,'AGT':7.42 ,
                        'ATC':24.29 ,'ACC':23.04 ,'AAC':20.02 ,'AGC':17.38 ,
                        'ATA':5.69 ,'ACA':6.04 ,'AAA':31.70 ,'AGA':2.62 ,
                        'ATG':27.31 ,'ACG':18.78 ,'AAG':11.31 ,'AGG':1.87 ,
                        'GTT':15.69 ,'GCT':12.98 ,'GAT':31.54 ,'GGT':17.37 ,
                        'GTC':18.06 ,'GCC':28.79 ,'GAC':19.99 ,'GGC':35.08 ,
                        'GTA':11.44 ,'GCA':13.13 ,'GAA':34.82 ,'GGA':8.86 ,
                        'GTG':24.71 ,'GCG':41.89 ,'GAG':20.31 ,'GGG':11.93},

        'Haemophilus_influenzae':{'TTT':32.54 ,'TCT':16.46 ,'TAT':24.57 ,'TGT':6.81 , #16
                        'TTC':12.15 ,'TCC':4.50 ,'TAC':7.11 ,'TGC':3.69 ,
                        'TTA':49.05 ,'TCA':12.63 ,'TAA':2.77 ,'TGA':0.59 ,
                        'TTG':18.26 ,'TCG':4.19 ,'TAG':0.63 ,'TGG':11.47 ,
                        'CTT':19.91 ,'CCT':12.31 ,'CAT':13.29 ,'CGT':23.49 ,
                        'CTC':5.28 ,'CCC':2.89 ,'CAC':7.19 ,'CGC':10.05 ,
                        'CTA':6.90 ,'CCA':16.37 ,'CAA':38.78 ,'CGA':5.39 ,
                        'CTG':4.65 ,'CCG':5.19 ,'CAG':7.72 ,'CGG':1.31 ,
                        'ATT':49.98 ,'ACT':15.92 ,'AAT':36.91 ,'AGT':12.10 ,
                        'ATC':14.38 ,'ACC':11.23 ,'AAC':12.66 ,'AGC':8.80 ,
                        'ATA':6.32 ,'ACA':15.84 ,'AAA':56.53 ,'AGA':4.07 ,
                        'ATG':23.64 ,'ACG':8.88 ,'AAG':7.30 ,'AGG':0.61 ,
                        'GTT':20.49 ,'GCT':20.96 ,'GAT':42.19 ,'GGT':27.93 ,
                        'GTC':6.44 ,'GCC':10.79 ,'GAC':7.83 ,'GGC':19.38 ,
                        'GTA':18.63 ,'GCA':32.50 ,'GAA':53.77 ,'GGA':10.85 ,
                        'GTG':20.17 ,'GCG':16.71 ,'GAG':10.81 ,'GGG':7.25},

        'Lactococcus_lactis':{'TTT':36.72 ,'TCT':16.63 ,'TAT':27.73 ,'TGT':3.88 , #17
                        'TTC':11.74 ,'TCC':3.18 ,'TAC':7.80 ,'TGC':1.33 ,
                        'TTA':32.61 ,'TCA':21.83 ,'TAA':2.97 ,'TGA':1.27 ,
                        'TTG':21.06 ,'TCG':3.70 ,'TAG':0.78 ,'TGG':10.17 ,
                        'CTT':25.09 ,'CCT':11.54 ,'CAT':13.10 ,'CGT':14.27 ,
                        'CTC':7.59 ,'CCC':2.61 ,'CAC':4.52 ,'CGC':4.04 ,
                        'CTA':7.90 ,'CCA':14.97 ,'CAA':30.81 ,'CGA':5.55 ,
                        'CTG':5.68 ,'CCG':2.79 ,'CAG':6.24 ,'CGG':2.22 ,
                        'ATT':52.51 ,'ACT':20.74 ,'AAT':41.04 ,'AGT':14.66 ,
                        'ATC':15.23 ,'ACC':6.84 ,'AAC':10.92 ,'AGC':5.86 ,
                        'ATA':9.68 ,'ACA':22.39 ,'AAA':61.63 ,'AGA':8.64 ,
                        'ATG':24.70 ,'ACG':6.58 ,'AAG':12.93 ,'AGG':1.70 ,
                        'GTT':31.89 ,'GCT':29.75 ,'GAT':38.20 ,'GGT':24.11 ,
                        'GTC':11.23 ,'GCC':10.80 ,'GAC':13.78 ,'GGC':8.10 ,
                        'GTA':13.15 ,'GCA':23.27 ,'GAA':56.68 ,'GGA':24.33 ,
                        'GTG':8.72 ,'GCG':7.71 ,'GAG':12.14 ,'GGG':7.84},

        'Actinobacillus_succinogenes':{'TTT':25.35, 'TCT':7.69 ,'TAT':21.03 ,'TGT':5.43, #18
                        'TTC':18.54, 'TCC':12.60 ,'TAC':10.01 ,'TGC':4.96,
                        'TTA':35.70 ,'TCA':6.14 ,'TAA':2.31 ,'TGA':0.46,
                        'TTG':26.20 ,'TCG':9.29 ,'TAG':0.46 ,'TGG':11.68,
                        'CTT':10.22, 'CCT':7.03 ,'CAT':12.74 ,'CGT':19.17,
                        'CTC':5.41 ,'CCC':5.08 ,'CAC':7.69 ,'CGC':14.20,
                        'CTA':3.70 ,'CCA':2.45 ,'CAA':26.95 ,'CGA':4.56,
                        'CTG':21.97 ,'CCG':23.88 ,'CAG':16.69 ,'CGG':7.27,
                        'ATT':39.34 ,'ACT':9.68 ,'AAT':27.32 ,'AGT':9.56,
                        'ATC':22.72 ,'ACC':21.10 ,'AAC':18.76 ,'AGC':10.77,
                        'ATA':4.83 ,'ACA':7.46 ,'AAA':50.69 ,'AGA':2.66,
                        'ATG':24.61 ,'ACG':15.02 ,'AAG':7.00 ,'AGG':0.52,
                        'GTT':15.05 ,'GCT':13.68 ,'GAT':34.19 ,'GGT':26.59,
                        'GTC':11.34 ,'GCC':24.69 ,'GAC':17.37 ,'GGC':24.76,
                        'GTA':17.34 ,'GCA':18.65 ,'GAA':52.12 ,'GGA':10.91,
                        'GTG':25.70 ,'GCG':32.15 ,'GAG':8.97 ,'GGG':7.59},

        'Geobacillus_stearothermophilus':{'TTT':27.17 ,'TCT':3.75 ,'TAT':18.23 ,'TGT':2.51 , #19
                        'TTC':14.37 ,'TCC':8.28 ,'TAC':15.35 ,'TGC':6.32 ,
                        'TTA':9.46 ,'TCA':5.04 ,'TAA':1.73 ,'TGA':1.77 ,
                        'TTG':28.54 ,'TCG':16.16 ,'TAG':0.68 ,'TGG':12.82 ,
                        'CTT':18.31 ,'CCT':3.58 ,'CAT':15.12 ,'CGT':7.49 ,
                        'CTC':18.45 ,'CCC':3.54 ,'CAC':8.09 ,'CGC':27.65 ,
                        'CTA':2.91 ,'CCA':6.08 ,'CAA':21.42 ,'CGA':5.62 ,
                        'CTG':21.25 ,'CCG':29.34 ,'CAG':14.87 ,'CGG':17.12 ,
                        'ATT':27.45 ,'ACT':2.89 ,'AAT':10.85 ,'AGT':2.82 ,
                        'ATC':32.03 ,'ACC':10.06 ,'AAC':19.75 ,'AGC':13.06 ,
                        'ATA':2.64 ,'ACA':8.88 ,'AAA':38.77 ,'AGA':2.88 ,
                        'ATG':25.09 ,'ACG':27.77 ,'AAG':17.69 ,'AGG':2.34 ,
                        'GTT':12.17 ,'GCT':11.29 ,'GAT':24.32 ,'GGT':6.58 ,
                        'GTC':30.27 ,'GCC':30.80 ,'GAC':24.36 ,'GGC':33.10 ,
                        'GTA':5.26 ,'GCA':9.23 ,'GAA':44.06 ,'GGA':13.46 ,
                        'GTG':28.30 ,'GCG':40.87 ,'GAG':29.72 ,'GGG':18.22},

        'Pyrococcus_horikoshii':{'TTT':22.52 ,'TCT':7.94 ,'TAT':19.79 ,'TGT':2.89 , #20
                        'TTC':20.87 ,'TCC':7.32 ,'TAC':19.31 ,'TGC':2.61 ,
                        'TTA':21.25 ,'TCA':9.73 ,'TAA':1.41 ,'TGA':1.79 ,
                        'TTG':11.76 ,'TCG':3.80 ,'TAG':0.86 ,'TGG':11.89 ,
                        'CTT':25.03 ,'CCT':10.94 ,'CAT':7.83 ,'CGT':1.05 ,
                        'CTC':16.79 ,'CCC':9.39 ,'CAC':6.86 ,'CGC':0.85 ,
                        'CTA':18.42 ,'CCA':16.86 ,'CAA':8.04 ,'CGA':0.98 ,
                        'CTG':8.64 ,'CCG':4.89 ,'CAG':8.65 ,'CGG':0.70 ,
                        'ATT':25.19 ,'ACT':12.74 ,'AAT':18.02 ,'AGT':9.90 ,
                        'ATC':15.90 ,'ACC':9.81 ,'AAC':16.37 ,'AGC':11.49 ,
                        'ATA':46.60 ,'ACA':10.32 ,'AAA':32.42 ,'AGA':21.04 ,
                        'ATG':22.62 ,'ACG':9.96 ,'AAG':48.43 ,'AGG':31.28 ,
                        'GTT':36.77 ,'GCT':21.39 ,'GAT':32.81 ,'GGT':15.60 ,
                        'GTC':11.58 ,'GCC':16.37 ,'GAC':12.11 ,'GGC':6.82 ,
                        'GTA':18.05 ,'GCA':18.65 ,'GAA':43.38 ,'GGA':33.03 ,
                        'GTG':11.51 ,'GCG':7.33 ,'GAG':44.73 ,'GGG':16.14},

        'Methanosarcina_barkeri':{'TTT':26.58 ,'TCT':14.41 ,'TAT':21.97 ,'TGT':6.52 , #21
                        'TTC':16.11 ,'TCC':11.37 ,'TAC':13.71 ,'TGC':6.24 ,
                        'TTA':11.86 ,'TCA':14.92 ,'TAA':2.35 ,'TGA':1.71 ,
                        'TTG':8.78 ,'TCG':5.99 ,'TAG':0.41 ,'TGG':10.13 ,
                        'CTT':33.50 ,'CCT':16.59 ,'CAT':10.13 ,'CGT':4.76 ,
                        'CTC':14.06 ,'CCC':6.58 ,'CAC':6.50 ,'CGC':4.19 ,
                        'CTA':6.45 ,'CCA':9.18 ,'CAA':8.53 ,'CGA':3.57 ,
                        'CTG':19.74 ,'CCG':6.60 ,'CAG':17.61 ,'CGG':3.97 ,
                        'ATT':33.08 ,'ACT':17.07 ,'AAT':28.03 ,'AGT':12.06 ,
                        'ATC':20.51 ,'ACC':12.78 ,'AAC':19.33 ,'AGC':10.45 ,
                        'ATA':23.64 ,'ACA':17.75 ,'AAA':48.05 ,'AGA':14.74 ,
                        'ATG':23.42 ,'ACG':6.42 ,'AAG':22.98 ,'AGG':12.00 ,
                        'GTT':23.23 ,'GCT':20.81 ,'GAT':31.47 ,'GGT':14.57 ,
                        'GTC':12.35 ,'GCC':14.13 ,'GAC':20.02 ,'GGC':13.34 ,
                        'GTA':20.68 ,'GCA':26.82 ,'GAA':53.57 ,'GGA':29.21 ,
                        'GTG':11.51 ,'GCG':5.81 ,'GAG':22.49 ,'GGG':12.65},

        'Thermotoga_maritima':{'TTT':18.87 ,'TCT':12.32 ,'TAT':9.75 ,'TGT':4.09 , #22
                        'TTC':32.87 ,'TCC':12.00 ,'TAC':25.93 ,'TGC':2.91 ,
                        'TTA':3.44 ,'TCA':7.20 ,'TAA':0.63 ,'TGA':2.27 ,
                        'TTG':12.50 ,'TCG':7.87 ,'TAG':0.29 ,'TGG':11.05 ,
                        'CTT':23.94 ,'CCT':9.86 ,'CAT':5.25 ,'CGT':2.93 ,
                        'CTC':33.09 ,'CCC':10.40 ,'CAC':10.55 ,'CGC':2.21 ,
                        'CTA':2.68 ,'CCA':9.62 ,'CAA':4.82 ,'CGA':3.01 ,
                        'CTG':24.79 ,'CCG':9.85 ,'CAG':15.28 ,'CGG':1.79 ,
                        'ATT':13.90 ,'ACT':7.47 ,'AAT':10.83 ,'AGT':7.94 ,
                        'ATC':29.83 ,'ACC':12.32 ,'AAC':25.19 ,'AGC':8.97 ,
                        'ATA':27.98 ,'ACA':12.85 ,'AAA':42.33 ,'AGA':29.35 ,
                        'ATG':22.94 ,'ACG':12.55 ,'AAG':33.56 ,'AGG':15.75 ,
                        'GTT':27.00 ,'GCT':14.19 ,'GAT':28.00 ,'GGT':19.47 ,
                        'GTC':16.20 ,'GCC':14.29 ,'GAC':21.49 ,'GGC':8.10 ,
                        'GTA':9.71 ,'GCA':14.46 ,'GAA':56.16 ,'GGA':32.35 ,
                        'GTG':33.68 ,'GCG':15.42 ,'GAG':32.70 ,'GGG':8.94},

        'Corynebacterium_diphtheriae':{'TTT':15.42 ,'TCT':11.73 ,'TAT':9.41 ,'TGT':3.62 , #23
                        'TTC':18.04 ,'TCC':14.88 ,'TAC':13.11 ,'TGC':5.56 ,
                        'TTA':6.44 ,'TCA':7.63 ,'TAA':1.77 ,'TGA':0.99 ,
                        'TTG':20.34 ,'TCG':12.57 ,'TAG':1.16 ,'TGG':14.32 ,
                        'CTT':17.46 ,'CCT':11.42 ,'CAT':9.84 ,'CGT':17.91 ,
                        'CTC':21.35 ,'CCC':10.70 ,'CAC':14.37 ,'CGC':23.34 ,
                        'CTA':9.22 ,'CCA':15.36 ,'CAA':16.36 ,'CGA':8.23 ,
                        'CTG':19.15 ,'CCG':11.46 ,'CAG':18.40 ,'CGG':6.26 ,
                        'ATT':22.15 ,'ACT':12.79 ,'AAT':11.68 ,'AGT':5.51 ,
                        'ATC':30.12 ,'ACC':25.82 ,'AAC':18.98 ,'AGC':10.61 ,
                        'ATA':3.48 ,'ACA':10.03 ,'AAA':16.16 ,'AGA':2.34 ,
                        'ATG':22.23 ,'ACG':12.41 ,'AAG':20.03 ,'AGG':2.57 ,
                        'GTT':21.48 ,'GCT':27.51 ,'GAT':29.60 ,'GGT':23.94 ,
                        'GTC':19.22 ,'GCC':28.71 ,'GAC':27.97 ,'GGC':30.48 ,
                        'GTA':13.12 ,'GCA':29.83 ,'GAA':31.13 ,'GGA':14.24 ,
                        'GTG':28.37 ,'GCG':23.56 ,'GAG':26.74 ,'GGG':9.36},

        'Propionibacterium_freudenreichii':{'TTT':1.36 ,'TCT':1.22 ,'TAT':5.70 ,'TGT':1.44 , #24
                        'TTC':27.83 ,'TCC':14.47 ,'TAC':14.10 ,'TGC':6.87 ,
                        'TTA':0.16 ,'TCA':3.69 ,'TAA':0.13  ,'TGA':2.45 ,
                        'TTG':8.23 ,'TCG':20.92 ,'TAG':0.53 ,'TGG':15.02 ,
                        'CTT':3.87 ,'CCT':2.37 ,'CAT':6.44 ,'CGT':10.14 ,
                        'CTC':25.90 ,'CCC':25.82 ,'CAC':17.08 ,'CGC':37.33 ,
                        'CTA':0.89 ,'CCA':3.05 ,'CAA':3.64 ,'CGA':4.12 ,
                        'CTG':56.92 ,'CCG':24.43 ,'CAG':30.27 ,'CGG':15.57 ,
                        'ATT':3.46 ,'ACT':2.25 ,'AAT':5.91 ,'AGT':2.84 ,
                        'ATC':40.36 ,'ACC':39.12 ,'AAC':17.46 ,'AGC':15.31 ,
                        'ATA':0.39 ,'ACA':2.50 ,'AAA':1.04 ,'AGA':0.70 ,
                        'ATG':22.63 ,'ACG':16.12 ,'AAG':23.87 ,'AGG':3.52 ,
                        'GTT':3.89 ,'GCT':5.46 ,'GAT':15.05 ,'GGT':12.20 ,
                        'GTC':28.69 ,'GCC':71.69 ,'GAC':47.87 ,'GGC':50.92 ,
                        'GTA':1.44 ,'GCA':11.46 ,'GAA':10.25 ,'GGA':9.83 ,
                        'GTG':50.98 ,'GCG':34.51 ,'GAG':40.56 ,'GGG':15.73},

        'Shigella_flexneri':{'TTT':21.79 ,'TCT':8.94 ,'TAT':16.14 ,'TGT':5.72 , #25
                        'TTC':16.47 ,'TCC':9.10 ,'TAC':12.40 ,'TGC':7.09 ,
                        'TTA':13.81 ,'TCA':7.91 ,'TAA':2.51 ,'TGA':1.92 ,
                        'TTG':13.62 ,'TCG':9.11 ,'TAG':0.52 ,'TGG':15.69 ,
                        'CTT':11.25 ,'CCT':7.29 ,'CAT':12.98 ,'CGT':20.93 ,
                        'CTC':10.88 ,'CCC':5.61 ,'CAC':9.90 ,'CGC':21.75 ,
                        'CTA':4.23 ,'CCA':8.58 ,'CAA':15.15 ,'CGA':4.03 ,
                        'CTG':51.56 ,'CCG':22.54 ,'CAG':28.69 ,'CGG':6.36 ,
                        'ATT':29.25 ,'ACT':9.27 ,'AAT':17.87 ,'AGT':9.11 ,
                        'ATC':24.58 ,'ACC':22.17 ,'AAC':21.46 ,'AGC':15.83 ,
                        'ATA':5.26 ,'ACA':7.74 ,'AAA':33.94 ,'AGA':2.94 ,
                        'ATG':27.45 ,'ACG':14.30 ,'AAG':11.13 ,'AGG':2.09 ,
                        'GTT':18.23 ,'GCT':15.91 ,'GAT':31.04 ,'GGT':23.93 ,
                        'GTC':14.84 ,'GCC':24.63 ,'GAC':18.83 ,'GGC':28.10 ,
                        'GTA':11.01 ,'GCA':20.41 ,'GAA':38.89 ,'GGA':8.64 ,
                        'GTG':25.63 ,'GCG':32.01 ,'GAG':17.91 ,'GGG':11.11},

        'Bacillus_subtilis':{'TTT':30.56 ,'TCT':12.68 ,'TAT':22.53 ,'TGT':3.70 , #26
                        'TTC':14.33 ,'TCC':8.08 ,'TAC':12.04 ,'TGC':4.47 ,
                        'TTA':19.01 ,'TCA':14.78 ,'TAA':2.36 ,'TGA':1.16 ,
                        'TTG':15.42 ,'TCG':6.45 ,'TAG':0.59 ,'TGG':10.31 ,
                        'CTT':23.02 ,'CCT':10.33 ,'CAT':15.20 ,'CGT':7.51 ,
                        'CTC':10.93 ,'CCC':3.29 ,'CAC':7.48 ,'CGC':8.61 ,
                        'CTA':4.82 ,'CCA':6.87 ,'CAA':19.53 ,'CGA':4.13 ,
                        'CTG':23.42 ,'CCG':16.09 ,'CAG':18.76 ,'CGG':6.65 ,
                        'ATT':36.68 ,'ACT':8.50 ,'AAT':21.98 ,'AGT':6.56 ,
                        'ATC':27.07 ,'ACC':8.72 ,'AAC':17.16 ,'AGC':14.17 ,
                        'ATA':9.41 ,'ACA':22.07 ,'AAA':49.32 ,'AGA':10.81 ,
                        'ATG':27.01 ,'ACG':14.78 ,'AAG':21.04 ,'AGG':3.98 ,
                        'GTT':18.88 ,'GCT':18.77 ,'GAT':32.59 ,'GGT':12.57 ,
                        'GTC':17.44 ,'GCC':16.02 ,'GAC':18.63 ,'GGC':23.48 ,
                        'GTA':13.19 ,'GCA':21.31 ,'GAA':48.79 ,'GGA':21.63 ,
                        'GTG':17.72 ,'GCG':20.50 ,'GAG':22.93 ,'GGG':11.16},

        'Thermoproteus_tenax':{'TTT':9.33 ,'TCT':8.31 ,'TAT':16.86 ,'TGT':3.35 , #27
                        'TTC':24.96 ,'TCC':11.97 ,'TAC':24.67 ,'TGC':4.57 ,
                        'TTA':7.99 ,'TCA':5.72 ,'TAA':1.02 ,'TGA':1.35 ,
                        'TTG':23.56 ,'TCG':13.78 ,'TAG':1.27 ,'TGG':13.22 ,
                        'CTT':10.94 ,'CCT':9.14 ,'CAT':4.02 ,'CGT':2.68 ,
                        'CTC':29.71 ,'CCC':19.55 ,'CAC':10.14 ,'CGC':8.69 ,
                        'CTA':12.14 ,'CCA':5.92 ,'CAA':7.81 ,'CGA':1.80 ,
                        'CTG':24.37 ,'CCG':15.89 ,'CAG':13.19 ,'CGG':4.58 ,
                        'ATT':8.39 ,'ACT':10.39 ,'AAT':7.38 ,'AGT':2.67 ,
                        'ATC':18.31 ,'ACC':10.81 ,'AAC':18.34 ,'AGC':13.92 ,
                        'ATA':35.79 ,'ACA':7.83 ,'AAA':14.25 ,'AGA':22.08 ,
                        'ATG':19.90 ,'ACG':13.40 ,'AAG':33.03 ,'AGG':30.19 ,
                        'GTT':11.48 ,'GCT':14.47 ,'GAT':16.53 ,'GGT':5.22 ,
                        'GTC':29.02 ,'GCC':44.48 ,'GAC':28.84 ,'GGC':41.04 ,
                        'GTA':14.57 ,'GCA':12.42 ,'GAA':9.78 ,'GGA':14.02 ,
                        'GTG':34.79 ,'GCG':29.13 ,'GAG':57.04 ,'GGG':18.01},

        'Saccharolobus_solfataricus':{'TTT':25.95 ,'TCT':14.99 ,'TAT':30.90 ,'TGT':4.41 , #28
                        'TTC':17.47 ,'TCC':7.53 ,'TAC':16.45 ,'TGC':2.16 ,
                        'TTA':41.20 ,'TCA':15.46 ,'TAA':2.57 ,'TGA':1.67 ,
                        'TTG':15.50 ,'TCG':4.80 ,'TAG':1.26 ,'TGG':10.56 ,
                        'CTT':15.13 ,'CCT':12.18 ,'CAT':8.30 ,'CGT':1.83 ,
                        'CTC':6.94 ,'CCC':5.39 ,'CAC':4.58 ,'CGC':0.64 ,
                        'CTA':19.24 ,'CCA':16.12 ,'CAA':15.77 ,'CGA':1.47 ,
                        'CTG':5.58 ,'CCG':4.17 ,'CAG':5.51 ,'CGG':0.64 ,
                        'ATT':33.31 ,'ACT':20.05 ,'AAT':32.51 ,'AGT':16.75 ,
                        'ATC':11.13 ,'ACC':6.83 ,'AAC':16.29 ,'AGC':7.27 ,
                        'ATA':49.70 ,'ACA':13.92 ,'AAA':40.07 ,'AGA':25.67 ,
                        'ATG':20.78 ,'ACG':6.41 ,'AAG':36.63 ,'AGG':17.26 ,
                        'GTT':27.14 ,'GCT':22.06 ,'GAT':34.16 ,'GGT':21.50 ,
                        'GTC':7.29 ,'GCC':7.21 ,'GAC':12.28 ,'GGC':6.80 ,
                        'GTA':28.16 ,'GCA':19.31 ,'GAA':38.80 ,'GGA':26.04 ,
                        'GTG':11.92 ,'GCG':7.20 ,'GAG':29.23 ,'GGG':9.94},

        'Bacillus_anthracis':{'TTT':32.61 ,'TCT':15.20 ,'TAT':28.06 ,'TGT':6.49 , #29
                        'TTC':14.20 ,'TCC':3.37 ,'TAC':9.47 ,'TGC':2.42 ,
                        'TTA':48.84 ,'TCA':14.56 ,'TAA':3.33 ,'TGA':1.16 ,
                        'TTG':10.08 ,'TCG':4.68 ,'TAG':1.10 ,'TGG':10.37 ,
                        'CTT':17.86 ,'CCT':8.96 ,'CAT':16.32 ,'CGT':13.65 ,
                        'CTC':4.25 ,'CCC':1.14 ,'CAC':4.90 ,'CGC':4.68 ,
                        'CTA':10.69 ,'CCA':15.93 ,'CAA':29.93 ,'CGA':5.60 ,
                        'CTG':3.55 ,'CCG':7.30 ,'CAG':7.09 ,'CGG':1.50 ,
                        'ATT':49.71 ,'ACT':12.35 ,'AAT':33.10 ,'AGT':14.47 ,
                        'ATC':13.27 ,'ACC':2.76 ,'AAC':13.31 ,'AGC':5.93 ,
                        'ATA':17.55 ,'ACA':27.30 ,'AAA':56.05 ,'AGA':9.93 ,
                        'ATG':27.85 ,'ACG':13.22 ,'AAG':18.69 ,'AGG':2.74 ,
                        'GTT':25.17 ,'GCT':20.32 ,'GAT':36.93 ,'GGT':23.99 ,
                        'GTC':5.69 ,'GCC':3.95 ,'GAC':8.68 ,'GGC':8.23 ,
                        'GTA':30.27 ,'GCA':28.65 ,'GAA':55.89 ,'GGA':23.88 ,
                        'GTG':10.73 ,'GCG':12.57 ,'GAG':18.18 ,'GGG':9.34},

        'Chlorobaculum_tepidum':{'TTT':13.52 ,'TCT':3.82 ,'TAT':9.37 ,'TGT':2.18 , #30
                        'TTC':29.94 ,'TCC':9.91 ,'TAC':18.49 ,'TGC':9.28 ,
                        'TTA':1.62 ,'TCA':5.95 ,'TAA':1.20 ,'TGA':2.04 ,
                        'TTG':8.86 ,'TCG':21.46 ,'TAG':0.32 ,'TGG':10.86 ,
                        'CTT':15.98 ,'CCT':4.97 ,'CAT':7.22 ,'CGT':7.58 ,
                        'CTC':36.11 ,'CCC':10.97 ,'CAC':13.98 ,'CGC':26.75 ,
                        'CTA':1.34 ,'CCA':5.15 ,'CAA':5.34 ,'CGA':4.13 ,
                        'CTG':35.57 ,'CCG':23.93 ,'CAG':25.88 ,'CGG':12.41 ,
                        'ATT':16.58 ,'ACT':4.13 ,'AAT':9.56 ,'AGT':3.91 ,
                        'ATC':42.64 ,'ACC':26.09 ,'AAC':24.50 ,'AGC':17.68 ,
                        'ATA':3.39 ,'ACA':4.80 ,'AAA':22.58 ,'AGA':4.06 ,
                        'ATG':26.75 ,'ACG':15.06 ,'AAG':29.95 ,'AGG':5.35 ,
                        'GTT':12.08 ,'GCT':12.09 ,'GAT':20.88 ,'GGT':12.77 ,
                        'GTC':26.04 ,'GCC':37.57 ,'GAC':30.64 ,'GGC':43.20 ,
                        'GTA':4.74 ,'GCA':12.12 ,'GAA':29.18 ,'GGA':10.33 ,
                        'GTG':27.96 ,'GCG':30.38 ,'GAG':37.63 ,'GGG':9.22},

        'Archaeoglobus_fulgidus':{'TTT':18.02 ,'TCT':6.43 ,'TAT':8.49 ,'TGT':2.59 , #31
                        'TTC':27.15 ,'TCC':9.13 ,'TAC':27.61 ,'TGC':9.22 ,
                        'TTA':5.82 ,'TCA':9.16 ,'TAA':1.43 ,'TGA':2.25 ,
                        'TTG':11.12 ,'TCG':7.01 ,'TAG':0.69 ,'TGG':10.43 ,
                        'CTT':23.39 ,'CCT':7.56 ,'CAT':4.27 ,'CGT':0.96 ,
                        'CTC':24.02 ,'CCC':11.29 ,'CAC':10.70 ,'CGC':2.43 ,
                        'CTA':4.95 ,'CCA':8.40 ,'CAA':3.04 ,'CGA':1.52 ,
                        'CTG':25.39 ,'CCG':10.99 ,'CAG':14.95 ,'CGG':1.27 ,
                        'ATT':28.18 ,'ACT':8.06 ,'AAT':10.05 ,'AGT':4.78 ,
                        'ATC':20.78 ,'ACC':11.63 ,'AAC':21.83 ,'AGC':18.60 ,
                        'ATA':22.91 ,'ACA':10.04 ,'AAA':24.57 ,'AGA':22.61 ,
                        'ATG':24.90 ,'ACG':11.82 ,'AAG':44.40 ,'AGG':29.77 ,
                        'GTT':38.40 ,'GCT':20.67 ,'GAT':22.36 ,'GGT':13.37 ,
                        'GTC':14.91 ,'GCC':19.54 ,'GAC':26.58 ,'GGC':17.57 ,
                        'GTA':10.58 ,'GCA':20.80 ,'GAA':28.79 ,'GGA':25.09 ,
                        'GTG':21.69 ,'GCG':16.57 ,'GAG':60.42 ,'GGG':16.07},

        'Thermoplasma_acidophilum':{'TTT':15.03 ,'TCT':10.23 ,'TAT':19.89 ,'TGT':1.25 , #32
                        'TTC':31.58 ,'TCC':15.82 ,'TAC':26.15 ,'TGC':4.52 ,
                        'TTA':3.93 ,'TCA':15.99 ,'TAA':0.83 ,'TGA':2.26 ,
                        'TTG':6.48 ,'TCG':12.34 ,'TAG':0.52 ,'TGG':8.53 ,
                        'CTT':20.15 ,'CCT':8.11 ,'CAT':8.63 ,'CGT':3.16 ,
                        'CTC':18.83 ,'CCC':6.32 ,'CAC':7.64 ,'CGC':5.17 ,
                        'CTA':6.23 ,'CCA':12.35 ,'CAA':2.87 ,'CGA':1.97 ,
                        'CTG':28.07 ,'CCG':12.40 ,'CAG':18.62 ,'CGG':2.63 ,
                        'ATT':12.40 ,'ACT':6.81 ,'AAT':18.85 ,'AGT':4.75 ,
                        'ATC':22.14 ,'ACC':13.37 ,'AAC':23.39 ,'AGC':15.86 ,
                        'ATA':55.60 ,'ACA':13.46 ,'AAA':17.51 ,'AGA':15.68 ,
                        'ATG':31.06 ,'ACG':13.98 ,'AAG':39.28 ,'AGG':26.34 ,
                        'GTT':21.05 ,'GCT':12.66 ,'GAT':37.43 ,'GGT':14.97 ,
                        'GTC':16.07 ,'GCC':20.76 ,'GAC':19.98 ,'GGC':24.94 ,
                        'GTA':14.11 ,'GCA':22.39 ,'GAA':28.36 ,'GGA':22.37 ,
                        'GTG':20.54 ,'GCG':13.63 ,'GAG':31.78 ,'GGG':9.98}

    }


    #Below are the lists of the species contained in the included usage files
    speciesList = ['Escherichia coli', 'Caenorhabditis elegans', 'Mus musculus', 'Homo sapien',
                    'Saccharomyces cerevisiae', 'Homo sapien (last)', 'Danio rerio',
                   'Gallus gallus', 'Xenopus tropicalis', 'Thermus thermophilus',
                   'Oryctolagus cuniculus', 'Rattus norvegicus', 'Vibrio harveyi',
                   'Pseudomonas aeruginosa', 'Salmonella enterica', 'Haemophilus influenzae',
                   'Lactococcus lactis', 'Actinobacillus succinogenes', 'Geobacillus stearothermophilus',
                   'Pyrococcus horikoshii', 'Methanosarcina barkeri', 'Thermotoga maritima',
                   'Corynebacterium diphtheriae', 'Propionibacterium freudenreichii', 'Shigella flexneri',
                   'Bacillus subtilis', 'Thermoproteus tenax', 'Saccharolobus solfataricus',
                   'Bacillus anthracis', 'Chlorobaculum tepidum', 'Archaeoglobus fulgidus',
                   'Thermoplasma acidophilum']
    speciesList2 = ['escherichia coli', 'caenorhabditis elegans', 'mus musculus',
                    'homo sapien', 'saccharomyces cerevisiae', 'homo sapien (last)',
                    'danio rerio', 'gallus gallus', 'xenopus tropicalis',
                    'thermus thermophilus', 'oryctolagus cuniculus', 'rattus norvegicus',
                    'vibrio harveyi', 'pseudomonas aeruginosa', 'salmonella enterica',
                    'haemophilus influenzae', 'lactococcus lactis', 'actinobacillus succinogenes',
                    'geobacillus stearothermophilus', 'pyrococcus horikoshii', 'methanosarcina barkeri',
                    'thermotoga maritima', 'corynebacterium diphtheriae', 'propionibacterium freudenreichii',
                    'shigella flexneri', 'bacillus subtilis', 'thermoproteus tenax',
                    'saccharolobus solfataricus', 'bacillus anthracis', 'chlorobaculum tepidum',
                    'archaeoglobus fulgidus', 'thermoplasma acidophilum',
                    '1', '2', '3', '4','5', '6', '7','8','9', '10','11','12','13', '14', '15', '16', '17',
                    '18', '19', '20', '21','22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32']


    #inputDict allows numbers or names to be entered for ease of use
    inputDict = {'escherichia coli':'Escherichia_coli', '1': 'Escherichia_coli',
                'caenorhabditis elegans': 'Caenorhabditis_elegans', '2':'Caenorhabditis_elegans',
                'mus musculus':'Mus_musculus', '3': 'Mus_musculus',
                'homo sapien':'Homo_sapien', '4':'Homo_sapien',
                'saccharomyces cerevisiae':'Saccharomyces_cerevisiae', '5':'Saccharomyces_cerevisiae',
                'homo sapien (last)':'Homo_sapien_last', '6':'Homo_sapien_last',
                'danio rerio':'Danio_rerio', '7':'Danio_rerio',
                'gallus gallus':'Gallus_gallus', '8':'Gallus_gallus',
                'xenopus tropicalis':'Xenopus_tropicalis', '9':'Xenopus_tropicalis',
                'thermus thermophilus':'Thermus_thermophilus', '10':'Thermus_thermophilus',
                'oryctolagus cuniculus':'Oryctolagus_cuniculus', '11':'Oryctolagus_cuniculus',
                'rattus norvegicus':'Rattus_norvegicus', '12':'Rattus_norvegicus',
                'vibrio harveyi':'Vibrio_harveyi', '13':'Vibrio_harveyi',
                'pseudomonas aeruginosa':'Pseudomonas_aeruginosa', '14':'Pseudomonas_aeruginosa',
                'salmonella enterica':'Salmonella_enterica', '15':'Salmonella_enterica',
                 'haemophilus influenzae':'Haemophilus_influenzae', '16':'Haemophilus_influenzae',
                 'lactococcus lactis':'Lactococcus_lactis', '17':'Lactococcus_lactis',
                 'actinobacillus succinogenes':'Actinobacillus_succinogenes', '18':'Actinobacillus_succinogenes',
                 'geobacillus stearothermophilus':'Geobacillus_stearothermophilus', '19':'Geobacillus_stearothermophilus',
                 'pyrococcus horikoshii':'Pyrococcus_horikoshii', '20':'Pyrococcus_horikoshii',
                 'methanosarcina barkeri':'Methanosarcina_barkeri', '21':'Methanosarcina_barkeri',
                 'thermotoga maritima':'Thermotoga_maritima', '22':'Thermotoga_maritima',
                 'corynebacterium diphtheriae':'Corynebacterium_diphtheriae', '23':'Corynebacterium_diphtheriae',
                 'propionibacterium freudenreichii':'Propionibacterium_freudenreichii', '24':'Propionibacterium_freudenreichii',
                 'shigella flexneri':'Shigella_flexneri', '25':'Shigella_flexneri',
                 'bacillus subtilis':'Bacillus_subtilis', '26':'Bacillus_subtilis',
                 'thermoproteus tenax':'Thermoproteus_tenax', '27': 'Thermoproteus_tenax',
                'saccharolobus solfataricus':'Saccharolobus_solfataricus', '28':'Saccharolobus_solfataricus',
                 'bacillus anthracis':'Bacillus_anthracis', '29':'Bacillus_anthracis',
                 'chlorobaculum tepidum':'Chlorobaculum_tepidum', '30':'Chlorobaculum_tepidum',
                 'archaeoglobus fulgidus':'Archaeoglobus_fulgidus', '31':'Archaeoglobus_fulgidus',
                 'thermoplasma acidophilum':'Thermoplasma_acidophilum', '32':'Thermoplasma_acidophilum'}




    def MinMaxfun(archivo, ventana, path_figuras):


        with open(str(args.inputFolderPath) + '/' + str(archivo), "r") as f:
            lineas = f.readlines()
            titulo = lineas[0]
            sequence = lineas[1]


        ###
        user_filename = str(archivo)[:-6]
        ###
        #Sequence Input Handling
        ###sequence = input('Nucleotide Sequence (start to stop codon): ')
        sequence = sequence.upper()

        if sequence[0:3] != "ATG":
            ###response = input("Sequence does not begin with ATG, would you re-enter the sequence (yes or no): ")
            ###response = response.lower()
            response = "no"
            while response == "yes":
                sequence = input("Please reenter the codon sequence: ")
                if sequence[0:3] != "ATG":
                    response = input("Sequence does not begin with ATG, would you like to re-enter the sequence (yes or no): ")
                else:
                    response = "no"
        while len(sequence)%3 != 0:
            ###sequence = input("The entered nucleotide sequence is not divisible by three, please reenter the sequence:")
            ###sequence = sequence.upper()
            sequence = sequence[:-1]



        ###windowSize = int(input('Sliding Window Length (odd # only): ')) #size of the sliding window for %MinMax
        windowSize = ventana
        while windowSize%2 != 1:
            windowSize = int(input('ERROR: Please enter an odd numbered window size: '))
        numRRT = 200 #number of random reverse translations to simulate

        #Usage file input handling
        print("")
        ###print("You may either use one of the following preloaded codon frequency tables or input your own.")
        ###for i in range(int(len(speciesList))):
            ###print(i+1, speciesList[i])
        ###usageChoice = input("Would you like to use one of those tables (yes or no): ")
        usageChoice = 'yes'
        usageChoice = usageChoice.lower()
        while usageChoice != "yes" and usageChoice != "no":
            usageChoice = input("Unrecognized input, please respond 'yes' or 'no': ")
            usageChoice = usageChoice.lower()

        #If the user wishes to use one of the included input files
        if usageChoice == "yes":
            print("")
            ###usageFile = input("Input species for usage file (name or number above): ")

            ### Tomar la especie y genero a partir del titulo de cada archivo fasta que se abre.
            x = titulo.split(' ')
            genero = x[1]
            especie = x[2]
            especie2 = genero + " " + especie
            ###
            # Determinar a que especie corresponde dicha secuencia de nucleotidos.
            default = 0
            if especie2.lower() == "escherichia coli":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '1'
            elif especie2.lower() == "caenorhabditis elegans":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '2'
            elif especie2.lower() == "mus musculus":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '3'
            elif especie2.lower() == "homo sapiens":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '6'
            elif especie2.lower() == "saccharomyces cerevisiae":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '5'
            elif especie2.lower() == "danio rerio":
                print("Especie ENOCNTRADA: " + especie2)
                usageFile = '7'
            elif especie2.lower() == "gallus gallus":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '8'
            elif especie2.lower() == "xenopus tropicalis":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '9'
            elif especie2.lower() == "thermus thermophilus":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '10'
            elif especie2.lower() == "oryctolagus cuniculus":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '11'
            elif especie2.lower() == "rattus norvegicus" or especie2.lower() == "rattus sp.":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '12'
            elif especie2.lower() == "vibrio harveyi":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '13'
            elif especie2.lower() == "pseudomonas aeruginosa":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '14'
            elif especie2.lower() == "salmonella enterica":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '15'
            elif especie2.lower() == "haemophilus influenzae":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '16'
            elif especie2.lower() == "lactococcus lactis":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '17'
            elif especie2.lower() == "actinobacillus succinogenes":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '18'
            elif especie2.lower() == "geobacillus stearothermophilus":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '19'
            elif especie2.lower() == "pyrococcus horikoshii":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '20'
            elif especie2.lower() == "methanosarcina barkeri":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '21'
            elif especie2.lower() == "thermotoga maritima":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '22'
            elif especie2.lower() == "corynebacterium diphtheriae" or especie2.lower() == "corynebacterium sp.":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '23'
            elif especie2.lower() == "propionibacterium freudenreichii":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '24'
            elif especie2.lower() == "shigella flexneri":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '25'
            elif especie2.lower() == "bacillus subtilis":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '26'
            elif especie2.lower() == "thermoproteus tenax":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '27'
            elif especie2.lower() == "saccharolobus solfataricus":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '28'
            elif especie2.lower() == "bacillus anthracis":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '29'
            elif especie2.lower() == "chlorobaculum tepidum":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '30'
            elif especie2.lower() == "archaeoglobus fulgidus":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '31'
            elif especie2.lower() == "thermoplasma acidophilum":
                print("Especie ENCONTRADA: " + especie2)
                usageFile = '32'
            else:
                print("Especie NO INCLUIDA: " + especie2 + ". DEFAULT: Escherichia coli")
                usageFile = '1'
                default = 1


            #usageFile = '3' ### Comentar esto!
            #default = 0 ### Comentar esto!

            usageFile = usageFile.lower()
            while usageFile not in speciesList2:
                print("")
                print("Unrecognized input, please enter either species name or index number: ")
                for i in range(int(len(speciesList2)/2)):
                    print(i+1, speciesList[i])
                usageFile = input("Input species for usage file (name or corresponding number): ")
                usageFile = usageFile.lower()
            freqDict = speciesDict[inputDict[usageFile]]

        #Or if they want to input one of their own
        else:
            print("")
            print("You have the option to manually input your own codon frequencies or to import a HIVE-CUT file.")
            manFile = input("'manually' or 'file': ")
            manFile = manFile.lower()
            while manFile != 'manually' and manFile != 'file':
                print("")
                manFile = input("Unrecognized input, please respond 'manually' or 'file': ")
                manFile = manFile.lower()

            #If the user wants to type by hand their own file
            if manFile == 'manually':
                print("")
                print("Please enter codon frequencies in the following format: <Codon><space><frequency><space>")
                print("e.g. ATG .18 CTG .46 GTG .3 TGT .01 ...")
                print("Do not hit enter until every codon has been input")
                print("Frequencies can be in any units (e.g. percent used, frequency per 1000, observed occurences, etc.)")
                inputFreq = input("Frequencies: ")
                inputFreq = inputFreq.upper()
                inputFreq = inputFreq.split()

                #Checks that input "codons" are actually codons and sequence has correct format
                codonFlag = 0
                correctCodons = 0
                while codonFlag == 0:
                    correctCodons = 0

                    while len(inputFreq)%2 !=0:
                        print("")
                        inputFreq = input("Error: Please enter frequencies in the following format: <Codon><space><frequency>")
                        inputFreq = inputFreq.upper()
                        inputFreq = inputFreq.split()

                    for i in range(0,len(inputFreq),2):
                        try:
                            mapDict[inputFreq[i]]
                            correctCodons += 1
                        except KeyError:
                            print("The codon " + inputFreq[i] + " was entered incorrectly.")


                    if correctCodons == int(len(inputFreq)/2):
                        codonFlag = 1
                    else:
                        print("")
                        print("Only enter codons as triplets comprised of ATCG")
                        inputFreq = input("Please reenter all codon frequencies correcting the above errors: ")
                        inputFreq = inputFreq.upper()
                        inputFreq = inputFreq.split()

                #Put input codons in kazusa like format
                newString = ""
                i = 0
                while i < len(inputFreq):
                    newString += str(inputFreq[i]) + " " + str(inputFreq[i+1]) + " " + "junk" + " "
                    i += 2
                frequenciesFile = []
                frequenciesFile.append(newString)

            #If the user wishes to include a file by typing the file path
            else:
                print("")
                print("A standard line for a HIVE-CUT uses NCBI's standard genetic code definition: ")
                print(" TTT 26.18 (76390)  TCT 23.35 (68138)  TAT 19.05 (55593)   TGT  7.82 (22826)")
                print('Ensure this is the format of your file')
                print('Copying and pasting from HIVE-CUT into a text file will preserve this format')
                print("The file should be saved as a '.txt' file")
                flag = 0
                while flag == 0:
                    filePath = input('Input file path: ')
                    try:
                        frequenciesFile = list(open(filePath))
                        flag = 1
                    except FileNotFoundError:
                        print("No file was found at this file path, please try again")


        print(user_filename + '.fasta')
        print("Calculating %MinMax...")

        #Data Cleaning, creates dictionaries (defined in lines 5-8) necessary for %MinMax and harmonization
        if usageChoice == "no":
            for line in frequenciesFile:
                line = line.split()
                i=0
                if len(line)>11:
                    while i < len(line):
                        freqDict[line[i]] = float(line[i+1])
                        aaMapDict[mapDict[line[i]]] = []
                        aaFreqDict[mapDict[line[i]]] = []
                        i+=3

            for line in frequenciesFile:
                line = line.split()
                i = 0
                if len(line)>11:
                    while i<len(line):
                        aaFreqDict[mapDict[line[i]]].append(float(line[i+1]))
                        aaMapDict[mapDict[line[i]]].append(str(line[i] + " " + line[i+1]))
                        i+=3
        else:
            for i in aaDict:#for each amino acid, initialize aaMapDict and aaFreqDict
                aaFreqDict[i] = []
                aaMapDict[i] = []
            for i in aaDict:#for each amino acid
                for j in aaDict[i]:#for each codon in that amino acid
                    aaFreqDict[i].append(freqDict[j])
                    aaMapDict[i].append(j + " " + str(freqDict[j]))

        #For a given input fasta sequence, break into codons and corresponding amino acids
        codonSeq = []

        extras = ""
        for line in sequence:
            line = line.rstrip()
            string = str(extras) + str(line)
            i=0
            j=3
            while j<=len(string):
                codonSeq.append(string[i:j])
                i+=3
                j+=3
            extras = str(string[i:])

        aaSeq = []
        for codon in codonSeq:
            aaSeq.append(mapDict[codon])

        #The function for calculating %MinMax
        def calculateMinMax(sequence, aaFreqDict, freqDict, mapDict, windowSize):
            freqDict = freqDict
            aaFreqDict = aaFreqDict
            windowSize = windowSize
            mapDict = mapDict
            codonSeq = sequence
            minMaxValues = [] #list to be returned of the %MinMax values

            for i in range(int(windowSize/2)): #%MinMax is undefined for the first and last (windowSize/2) condons
                minMaxValues.append(0)

            #Using the specified sliding window size (windowSize/2 - 1 on either side of the central codon), min/max is calculated
            for i in range(len(codonSeq)-windowSize+1):
                window = codonSeq[i:i+windowSize] #list of the codons in the current window

                Actual = 0.0     #average of the actual codon frequencies
                Max = 0.0        #average of the min codon frequencies
                Min = 0.0        #average of the max codon frequencies
                Avg = 0.0        #average of the averages of all the frequencies associated with each amino acid

                #Sum the frequencies
                for codon in window:
                    frequencies = aaFreqDict[mapDict[codon]] #list of all frequencies associated with the amino acid this codon encodes

                    Actual += freqDict[codon]
                    Max += max(frequencies)
                    Min += min(frequencies)
                    Avg += sum(frequencies)/len(frequencies)

                #Divide by the window size to get the averages
                Actual = Actual/windowSize
                Max = Max/windowSize
                Min = Min/windowSize
                Avg = Avg/windowSize

                percentMax = ((Actual-Avg)/(Max-Avg))*100
                percentMin = ((Avg-Actual)/(Avg-Min))*100

                if(percentMax >= 0):
                    minMaxValues.append(round(percentMax,2))
                else:
                    minMaxValues.append(round(-percentMin,2))

            #fills in values for codons where window size makes min/max unable to be calculated
            for i in range(int(windowSize/2)):
                minMaxValues.append(0)

            return minMaxValues



        minMaxValues = calculateMinMax(codonSeq, aaFreqDict, freqDict, mapDict, windowSize )


        print("Calculating Random Reverse Translations...")



        #Calculates RRTs which simulate weighted-random codon sequences for the given AA sequence

        #calculate RRTs
        rrtCodonList = []
        for i in range(numRRT):
            rrtCodonList.append([])
            for j in range(len(aaSeq)):
                aa = aaSeq[j]
                data = aaMapDict[aa]
                codons = []
                freqs = []
                for duo in data:
                    cAndF = duo.split()
                    codons.append(cAndF[0])
                    freqs.append(float(cAndF[1]))
                freqSum = sum(freqs)
                randomNum = random.uniform(0,freqSum)
                total = 0
                flag = 0
                for k in range(len(freqs)):
                    if flag == 0:
                        total += freqs[k]
                        if total > randomNum:
                            rrtCodonList[i].append(codons[k])
                            flag = 1
        rrtList = []
        for i in range(numRRT):
            codonSequence = rrtCodonList[i]
            mmValues = calculateMinMax(codonSequence, aaFreqDict, freqDict, mapDict, windowSize)
            rrtList.append(mmValues)

        #Average the RRTs
        avgRRT = []
        stdRRT = []
        for i in range(int(windowSize/2)):
            avgRRT.append(None)
            stdRRT.append(None)
        for i in range(int(windowSize/2),len(rrtList[0])-int(windowSize/2)):
            avg = 0
            stdlist = []
            for j in range(len(rrtList)):
                avg += rrtList[j][i]
                stdlist.append(rrtList[j][i])
            avg = avg/numRRT
            stdRRT.append(s.stdev(stdlist))
            avgRRT.append(round(avg,2))
        for i in range(int(windowSize/2)):
            avgRRT.append(None)
            stdRRT.append(None)

        num = [] #index for graphs later
        for i in range(1,len(minMaxValues)+1):
            num.append(i)








        #Simple Harmonization Algorithm, choses the codon of equal rank in guest organism, uses similar framework as above
        print("")
        ###harmonize = input("Would you like to harmonize this sequence to another organism? ('Yes' or 'No'): ")
        ###
        harmonize = 'No'
        ###
        harmonize = harmonize.lower()
        while harmonize != "yes" and harmonize != "no":
            harmonize = input("Unrecognized input, please enter 'Yes' or 'No': ")
            harmonize = harmonize.lower()

        #If a harmonized sequence is desired
        if harmonize == "yes":
            print("")
            print("You may either use one of the following preloaded codon frequency tables or input your own.")
            for i in range(int(len(speciesList))):
                print(i+1, speciesList[i])
            usageChoice2 = input("Would you like to use one of those tables (yes or no): ")
            usageChoice2 = usageChoice2.lower()
            while usageChoice2 != "yes" and usageChoice2 != "no":
                usageChoice2 = input("Unrecognized input, please respond 'yes' or 'no': ")
                usageChoice2 = usageChoice2.lower()

            #If the user would like to harmonize to one of the included frequency tables
            if usageChoice2 == "yes":
                print("")
                usageFile2 = input("Input species for usage file (name or number above): ")
                usageFile2 = usageFile2.lower()
                while usageFile2 not in speciesList2:
                    print("")
                    print("Unrecognized input, please enter either species name or index number: ")
                    for i in range(int(len(speciesList2)/2)):
                        print(i+1, speciesList[i])
                    usageFile2 = input("Input species for usage file (name or corresponding number): ")
                    usageFile2 = usageFile2.lower()
                freqDict2 = speciesDict[inputDict[usageFile2]]

            #If the user would like to harmonize to their own frequency table
            else:
                print("")
                print("Like above, you may either enter the path for a codon table, or enter the frequencies by hand.")
                harmonize2 = input("'Manually' or 'File': ")
                harmonize2 = harmonize2.lower()
                while harmonize2 != 'manually' and harmonize2 != 'file':
                    print("")
                    harmonize2 = input("Unrecognized input, please respond 'manually' or 'file': ")
                    harmonize2 = harmonize2.lower()
                if harmonize2 == 'manually':
                    print("")
                    print("Please enter codon frequencies in the following format: <Codon><space><frequency>")
                    print("e.g. ATG .18 CTG .46 GTG .3 TGT .01 ...")
                    print("Do not hit enter until every codon has been input")
                    print("Frequencies can be in any units (e.g. percent used, frequency per 1000, observed occurences, etc.)")
                    harmonize3 = input("Frequencies: ")
                    hamronize3 = harmonize3.upper()
                    harmonize3 = harmonize3.split()


                    #Check that input "codons" are actually codons and sequence has correct format
                    codonFlag = 0
                    correctCodons = 0
                    while codonFlag == 0:
                        correctCodons = 0

                        while len(harmonize3)%2 !=0:
                            print("")
                            harmonize3 = input("Error: Please enter frequencies in the following format: <Codon><space><frequency>")
                            harmonize3 = harmonize3.upper()
                            harmonize3 = harmonize3.split()

                        for i in range(0,len(harmonize3),2):
                            try:
                                mapDict[harmonize3[i]]
                                correctCodons += 1
                            except KeyError:
                                print("The codon " + harmonize3[i] + " was entered incorrectly.")


                        if correctCodons == int(len(harmonize3)/2):
                            codonFlag = 1
                        else:
                            print("")
                            print("Only enter codons as triplets comprised of ATCG")
                            harmonize3 = input("Please reenter all codon frequencies correcting the above errors: ")
                            harmonize3 = harmonize3.upper()
                            harmonize3 = harmonize3.split()

                    #Put input codons in kazusa like format
                    newString = ""
                    i = 0
                    while i < len(harmonize3):
                        newString += str(harmonize3[i]) + " " + str(harmonize3[i+1]) + " " + "junk" + " "
                        i += 2
                    frequenciesFile2 = []
                    frequenciesFile2.append(newString)
                else:
                    print("")
                    print("A standard line for a HIVE-CUT file looks like: ")
                    print("TTT 26.18 (76390)  TCT 23.35 (68138)  TAT 19.05 (55593)   TGT  7.82 (22826)")
                    print('Ensure this is the format of your file (and that the amino acid for each codon is not listed")')
                    print("The file should be saved as a '.txt' file")
                    flag = 0
                    while flag == 0:
                        filePath = input('Input file path: ')
                        try:
                            frequenciesFile2 = list(open(filePath))
                            flag = 1
                        except FileNotFoundError:
                            print("No file was found at this file path, please try again")



            #Data Cleaning
            if usageChoice2 == "no":
                freqDict2 = dict()
                aaMapDict2 = dict()#dictionary from amino acid to list of codons with frequencies for it (for RRTs)
                aaFreqDict2 = dict()
                for line in frequenciesFile2:
                    line = line.split()
                    i=0
                    if len(line)>11:
                        while i < len(line):
                            freqDict2[line[i]] = float(line[i+1])
                            aaMapDict2[mapDict[line[i]]] = []
                            aaFreqDict2[mapDict[line[i]]] = []
                            i+=3

                for line in frequenciesFile2:
                    line = line.split()
                    i = 0
                    if len(line)>11:
                        while i<len(line):
                            aaFreqDict2[mapDict[line[i]]].append(float(line[i+1]))
                            aaMapDict2[mapDict[line[i]]].append(str(line[i] + " " + line[i+1]))
                            i+=3

            else:
                aaMapDict2 = dict()#dictionary from amino acid to list of codons with frequencies for it (for RRTs)
                aaFreqDict2 = dict()
                for i in aaDict:#for each amino acid, initialize aaMapDict and aaFreqDict
                    aaFreqDict2[i] = []
                    aaMapDict2[i] = []
                for i in aaDict:#for each amino acid
                    for j in aaDict[i]:#for each codon in that amino acid
                        aaFreqDict2[i].append(freqDict2[j])
                        aaMapDict2[i].append(j + " " + str(freqDict2[j]))


            hCodonSeq = []
            hMinMax = []
            for i in codonSeq:
                midList1 = []
                midList2 = []
                origList = aaMapDict[mapDict[i]]
                for j in origList:
                    j = j.split()
                    midList1.append((float(j[1]),j[0]))
                midList1 = sorted(midList1)
                newList = aaMapDict2[mapDict[i]]
                for j in newList:
                    j = j.split()
                    midList2.append((float(j[1]),j[0]))
                midList2 = sorted(midList2)
                value = 0
                counter = 0
                for j in midList1:
                    if i == j[1]:
                        value = counter
                    counter = counter + 1
                hCodonSeq.append(midList2[value][1])
            hMinMax = calculateMinMax(hCodonSeq, aaFreqDict2, freqDict2, mapDict, windowSize)


        print("Creating graph, saving data to working directory...")



        #make graph and save it
        stdevbel = []
        stdevabv = []
        stdevabv0 = []
        stdevbel0 = []
        for i in range(len(avgRRT)):
            try:
                stdevbel.append(avgRRT[i]-2*stdRRT[i])
                stdevabv.append(avgRRT[i]+2*stdRRT[i])
                stdevbel0.append(avgRRT[i]-2*stdRRT[i])
                stdevabv0.append(avgRRT[i]+2*stdRRT[i])
            except TypeError:
                stdevbel0.append(0)
                stdevabv0.append(0)
                stdevbel.append(None)
                stdevabv.append(None)


        import matplotlib.pyplot as plt
        if harmonize == 'yes':
            #Change the zeros to "None" values for graphing
            for i in range(int(windowSize/2)):
                minMaxValues[i]=None
                hMinMax[i] = None
            for i in range(1,int(windowSize/2)+1):
                minMaxValues[-i] = None
                hMinMax[-i] = None

            #make csv and save it
            df = pd.DataFrame({"Codon #": num, "Codon": codonSeq, "AA": aaSeq, "%MinMax": minMaxValues, "Avg RRT": avgRRT, "STDev": stdRRT, "2 STDev Above": stdevabv, "2 STDev Below": stdevbel, "Harmonized Sequence": hCodonSeq, "Harmonized %MinMax": hMinMax})
            df = df.set_index("Codon #")
            df = df[["Codon", "AA","%MinMax","Avg RRT", "STDev", "2 STDev Above","2 STDev Below", "Harmonized Sequence", "Harmonized %MinMax"]]
            if default == 1:
                base_filename = 'MinMax_' + user_filename + '_' + especie2 + '_DEFAULT_win' + str(windowSize) + '.csv'
            else:
                base_filename = 'MinMax_' + user_filename + '_' + especie2 + '_win' + str(windowSize) + '.csv'
            df.to_csv(str(path_figuras) + '/' + base_filename)

            plt.close() #included if running in a jupyter notebook or other environment which saves variables
            plt.plot(num,minMaxValues, label = "%MinMax", color = '#0000CC')
            plt.plot(num, avgRRT, label = "Average RRT", color = '#808080')
            plt.plot(num, hMinMax, label = "Harmonized %MinMax", color = '#4C9900')
            plt.plot([0,len(minMaxValues)],[0,0],'r--')
            plt.legend(loc = 4)
            plt.xlim(0,len(minMaxValues))
            plt.ylim(-100,100)
            plt.fill_between(num, stdevbel0,stdevabv0, color = '#808080', alpha = .1)
            plt.xlabel("Center Of Codon Window")
            plt.ylabel("%MinMax")
            if default == 1:
                figure_filename = 'MinMax_' + user_filename + '_' + especie2 + '_DEFAULT_win' + str(windowSize) + '_Figure' + '.png'
            else:
                figure_filename = 'MinMax_' + user_filename + '_' + especie2 + '_win' + str(windowSize) + '_Figure' + '.png'

            plt.savefig(str(path_figuras) + '/' + figure_filename, dpi =300)
            print("Done")

        else:
            #Change the zeros to "None" values for graphing
            for i in range(int(windowSize/2)):
                minMaxValues[i]=None
            for i in range(1,int(windowSize/2)+1):
                minMaxValues[-i] = None

            df = pd.DataFrame({"Codon #": num, "Codon": codonSeq, "AA": aaSeq, "%MinMax": minMaxValues, "Avg RRT": avgRRT, "STDev": stdRRT, "2 STDev Above": stdevabv, "2 STDev Below": stdevbel})
            df = df.set_index("Codon #")
            df = df[["Codon", "AA","%MinMax","Avg RRT", "STDev", "2 STDev Above","2 STDev Below"]]
            if default == 1:
                base_filename = 'MinMax_' + user_filename + '_' + especie2 + '_DEFAULT_win' + str(windowSize) + '.csv'
            else:
                base_filename = 'MinMax_' + user_filename + '_' + especie2 + '_win' + str(windowSize) + '.csv'
            df.to_csv(str(path_figuras) + '/' + base_filename)

            plt.close()
            plt.plot(num,minMaxValues, label = "%MinMax", color = '#0000CC')
            plt.plot(num, avgRRT, label = "Average RRT", color = '#808080')
            plt.plot([0,len(minMaxValues)],[0,0],'r--')
            plt.legend(loc = 4)
            plt.xlim(0,len(minMaxValues))
            plt.ylim(-100,100)
            plt.fill_between(num, stdevbel0,stdevabv0, color = '#808080', alpha = .1)
            plt.xlabel("Center Of Codon Window")
            plt.ylabel("%MinMax")
            if default == 1:
                figure_filename = 'MinMax_' + user_filename + '_' + especie2 + '_DEFAULT_win' + str(windowSize) + '_Figure' + '.png'
            else:
                figure_filename = 'MinMax_' + user_filename + '_' + especie2 + '_win' + str(windowSize) + '_Figure' + '.png'
            plt.savefig(str(path_figuras) + '/' + figure_filename, dpi =300)
            print("Done")

    # Obtencion de la lista de todos los archivos fasta dentro de la carpeta
    filenames = os.listdir(str(args.inputFolderPath))
    #print(filenames)
    print("Cantidad de archivos: " + str(len(filenames)))
    ###user_filename = input('Name your file:')#Allows user to name their output file - will be saved to the working directory
    ventanas = []
    ventanas = list(map(int, input("Ingrese las longitudes de ventana deseadas. Separadas por espacios  (odd # only): ").split()))

    # Creacion de carpeta de figuras (output) si es que aun no existe.
    path = str(args.outputPath) + r'/figuras'
    if not os.path.exists(path):
        os.makedirs(path)





    for user_filename in filenames: # For para cada archivo dentro de la carpeta
        for i in range(len(ventanas)): # For para todas las ventanas que se quieran probar
            MinMaxfun(archivo = user_filename, ventana = int(ventanas[i]), path_figuras = path)
