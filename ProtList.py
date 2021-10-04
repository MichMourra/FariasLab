
#··················································································#
#··················································································#
#                                    How to run                                    #
#        python3 ProtList.py -i ./ScopDatabaseFile.txt -o ScopeNewList.txt         #
#··················································································#
#··················································································#


#··················································································#
#··················································································#
#                                      Modules                                     #
#··················································································#
#··················································································#


import re
import os
import argparse


#··················································································#
#··················································································#
#                                     Arguments                                    #
#··················································································#
#··················································································#

# Arguments that contains the paths to the
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputfile", help="Path to the file that contains the SCOP database proteins")
parser.add_argument("-o", "--outputfile", help="Path to the file that will contain the protein list")
args = parser.parse_args()

# Get values from the arguments
InputFile = args.inputfile
OutputFile = args.outputfile


#··················································································#
#··················································································#
#                                    Main Code                                     #
#··················································································#
#··················································································#



# Open the file that contains the proteins from SCOP database
File = open(InputFile,"r");
Result = open(OutputFile,"w")

# Find all the proteins in the file with a regular expression
Pattern = re.compile(r'>.+ \(\w:\)')
Count = 0

for line in File:
    if line.startswith(">"):
        try:
            # Once the program found a protein, this save it into a list file
            Match = Pattern.findall(line)
            Protein = Match[0].replace(">","")
            Protein = Protein.split(" ")
            Protein = Protein[0] + " " + Protein[1] + "\n"
            Result.write(Protein)
            Count += 1
        except:
            continue

# Show how many proteins the script fount
print("The script succesfully found {} proteins".format(Count))

# Close files
File.close()
Result.close()
