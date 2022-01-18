'''

|||    |||   |||||||||   |||            |||
|||    |||   |||   |||    |||          |||
||||||||||   |||   |||     ||          ||
||||||||||   |||   |||     |||   ||   |||
|||    |||   |||   |||      ||||||||||||
|||    |||   |||||||||        |||  |||


        |||||||||||  |||||||||
            |||      |||   |||
            |||      |||   |||
            |||      |||   |||
            |||      |||   |||
            |||      |||||||||

|||||||||||    |||     |||  ||||     |||
|||     |||    |||     |||  ||| ||   |||
|||||||||||    |||     |||  |||  ||  |||
|||    |||     |||     |||  |||   || |||
|||     |||    |||     |||  |||    |||||
|||      |||   |||||||||||  |||     ||||


Example:

./DownoloadFastas.py -I ./ENA.txt -O "./SCOPeFastas/" -E cmourra@lcg.unam.mx


||||||||||   ||||||||||  |||        |||    |||  |||       |||  ||||||||||
|||          |||    |||  |||        |||    |||  |||||   |||||  |||
|||          |||    |||  |||        |||    |||  ||| || || |||  ||||||||||
|||          |||    |||  |||        |||    |||  |||  |||  |||  ||||||||||
|||          |||    |||  |||        |||    |||  |||       |||         |||
||||||||||   ||||||||||  |||||||||  ||||||||||  |||       |||  ||||||||||

ENA File Columns

 ·· 01 ·· Test (UniProt Entry)
 ·· 02 ·· SCOPe Name
 ·· 03 ·· SCOPe Id
 ·· 04 ·· Blasth
 ·· 05 ·· V3
 ·· 06 ·· V4
 ·· 07 ·· V5
 ·· 08 ·· V6
 ·· 09 ·· V7
 ·· 10 ·· V8
 ·· 11 ·· V9
 ·· 12 ·· V10
 ·· 13 ·· V11
 ·· 14 ·· V12
 ·· 15 ·· Entry.name
 ·· 16 ·· Status
 ·· 17 ·· Protein.names
 ·· 18 ·· Gene.names
 ·· 19 ·· Organism
 ·· 20 ·· Length
 ·· 21 ·· Gene.names...ordered.locus..
 ·· 22 ·· Pathway
 ·· 23 ·· DNA.binding
 ·· 24 ·· Gene.ontology.GO.
 ·· 25 ·· Sequence.similarities
 ·· 26 ·· Cross.reference..Pfam.
 ·· 27 ·· Cross.reference..PANTHER.
 ·· 28 ·· Cross.reference..Reactome.
 ·· 29 ·· Cross.reference..UniPathway
 ·· 30 ·· Cross.reference..BioCyc
 ·· 31 ·· Ensembl.transcript
 ·· 32 ·· Cross.reference..KEGG.
 ·· 33 ·· Cross.reference..PATRIC.
 ·· 34 ·· EnsemblBacteria.transcript
 ·· 35 ·· Cross.reference..STRING.
 ·· 36 ·· Gene.ontology.IDs
 ·· 37 ·· Cross.reference..eggNOG.
 ·· 38 ·· Cross.reference..OrthoDB.
 ·· 39 ·· Gene.names...ORF..
 ·· 40 ·· Gene.names...primary..
 ·· 41 ·· Sequence
 ·· 42 ·· Function..CC.
 ·· 43 ·· Activity.regulation
 ·· 44 ·· Gene.names...synonym..
 ·· 45 ·· Organism.ID
 ·· 46 ·· Annotation
 ·· 47 ·· EnsemblPlants.transcript
 ·· 48 ·· Cross.reference..RefSeq
 ·· 49 ·· Cross.reference..PlantReactome.
 ·· 50 ·· Cross.reference.GeneID
 ·· 51 ·· Cross.reference..EMBL.
'''

#··················································································#
#··················································································#
#                                      Modules                                     #
#··················································································#
#··················································································#

import Bio
from Bio import Entrez

#··················································································#
#··················································································#
#                                     Arguments                                    #
#··················································································#
#··················································································#


# Arguments for request the minmax file paths, elongation profile paths and output path
parser = argparse.ArgumentParser()
parser.add_argument("-I", "--InputFilePath", help="Path to the file that contains the ENA proteins data")
parser.add_argument("-O", "--OutputFolder", help="Path to the folder that will contains the resulting files")
parser.add_argument("-E", "--Email", help="Path to the folder that contains the elongation profile files")
args = parser.parse_args()

# Get the content of the arguments
filepath = args.InputFilePath
PathMinMax = args.OutputFolder
mail = args.Email

#··················································································#
#··················································································#
#                                     Main Code                                    #
#··················································································#
#··················································································#

# Log in entrez
Entrez.email = mail


first = 1
# open file
with open(filepath,"r") as file:

    for line in file:
        # ignore first line
        if first == 1:
            first = 0
            continue
        # get the columns elements in the list
        line = line.replace("\n","").split("\t")
        # extract elements
        geneName = line[17]
        Organism = line[18].split(" (")[0]
        scope = line[1]
        protein = line[16]
        entry = line[0]
        idEMBL = line[50].split(";")
        idEMBL.remove('')
        # header of the file that contains the protein data
        header = ">{}\t{}\t{}\t{}\n".format(entry,geneName,Organism,protein)

        # loop through the ids on EMBL
        for count,element in enumerate(idEMBL):
            # if we the script could access to the id we write the sequence in a new file
            try:
                handle = Entrez.efetch(db="nucleotide", id=element, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                handle.close()
                sequence = str(record.seq)
                newpath ="{}{}({}).fasta".format(fastasDir,scope,count)
                # open the new file and write the data
                with open(newpath,"w") as newfile:
                    newfile.write(header)
                    for num,char in enumerate(sequence):
                        newfile.write(char)
                        if ((num%69 == 0) and (num!=0)):
                        newfile.write("\n")
            # if we can not access to the id we ignore this data and continue with the next id
            except:
                continue

    print("Archivos fasta descargados exitosamente")
