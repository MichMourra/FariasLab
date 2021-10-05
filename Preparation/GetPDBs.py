
#········································································································································#
#········································································································································#
#                                                 How to run                                                                             #
#      python3 ./GetPDBs.py -i ./descargas -o ../cg_cotrans/input_pdbs_calc_consensus/ -f ../cg_cotrans/input_fastas_calc_consensus/     #
#········································································································································#
#········································································································································#


#··················································································#
#··················································································#
#                                      Modules                                     #
#··················································································#
#··················································································#

import wget
import os
import argparse

#··················································································#
#··················································································#
#                                     Arguments                                    #
#··················································································#
#··················································································#


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--InputFolder", help="Path to the folder that contains the elongation profile files")
parser.add_argument("-o", "--OutputPath", help="Path to the folder that will contains the resulting files")
parser.add_argument("-f", "--FastasPath", help="Path to the folder that will contains the fasta files")
args = parser.parse_args()


InputFolder = args.InputFolder
OutputPath = args.OutputPath
FastasPath = args.FastasPath


#··················································································#
#··················································································#
#                                     Main Code                                    #
#··················································································#
#··················································································#


Files = os.listdir(InputFolder)
Proteins = []

for line in Files:
    FastaFile = "{}/{}".format(InputFolder,line)
    CopyCommand = "cp {} {}".format(FastaFile,FastasPath)
    os.popen(CopyCommand)
    # Obtener una lista de proteinas que queremos descargar
    line = line.replace("d","").replace(".fasta","")
    Prot = line[0:4]
    Proteins.append(Prot)

for protein in Proteins:
    try:
        # Generamos el url para poder descargar los PDBs
        url = "http://www.pdb.org/pdb/files/{}.pdb.gz".format(protein)
        wget.download(url, OutputPath)
    except:
        continue
