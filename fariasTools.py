'''
NAME
    FariasLab Tools
VERSION
    1.0

AUTORS
    Carlos Michel Mourra Diaz
    Daniel Alejandro Lopez
    Enrique Rojas

DESCRIPTION
    This module aims to provide a set of functions
    for the bulk protein analysis.

CATEGORY
    PROTEIN ANALYSES

SOFTWARE REQUIREMENTS
    Python 3.10

MODULES
    -- os
    -- re
    -- urllib
    -- wget
    -- subprocess
    -- shutil
    -- xmltramp2
    -- Bio
        |_ Seq.IO
        |_ Entrez

INFORMATION
    This Python module has been designed with the purpose of
    provide all the necesary tools to perform cotranslational
    folding analyses in proteins. In this module we can find
    the following functions:
        1. -- getScopeId -> A function which extract the scope ids from the astral scope file.
        2. -- GetIdsfromfile -> A function which aims to get the PDB ids from the list file.
        3. -- GetAminoFiles -> A function which aims to download the amino acid sequence of each PDB id into fasta format file.
        4. -- GetUniprotAccessions -> A function which aims to generate a dictionary with the following structure {pdb_id : uniprot entry}.
        5. -- GetEnaAccessions -> A function which aims to generate a dictionary with the following structure {pdb_id : ena ids}.
        6. -- Dictfromfile -> A function which aims to return a dictionary with the data extracted from certain file.
        7. -- AllIdsFile -> A function which aims to create a file with all the cross references ids.
        8. -- GetAminoLength -> A function which aims to generate a file that contains the PDB id and the length of the Sequence.
        9. -- GetNucleotideFiles -> A function which aims to download the nucleotide fasta files from the given pdb ids.
        10. -- blastDB -> A function aims to make a blast database for alignments.
        11. -- performBlast -> A function which aims to perform blast process for a selected set of files.
        12. -- GetRemoteAlignments -> A function which aims to download the nucleotide fastas corresponding for each pdb best alignment in the remote server.
        13. -- GetBestLocalAlignments -> A function which aims to download the nucleotide fasta files for the best alignment from a local database.
        14. -- GetPDBs -> A function which aims to download the pdb files from a set of proteins.
        15. -- TranslatedFiles -> A function which aims to translate the nucleotide fasta files from the given pdb ids.
        16. -- ProteinAlignment -> A function which aims to align 2 amino sequence files by using the program emboss_needle.
        17. -- SearchPDBlist -> A function which aims to download all the neccessary files to perform cotranslational analyses from a set of proteins.

USAGE EXAMPLES
    1. getScopeId("astral.fa.txt","scope-pdb.txt")
    2. AllIdsFile('scope-pdb.txt','uniprotTest.txt','enaTest.txt','allIdsTest.txt')
    3. GetIdsfromfile('PDBlist.txt')
    4. GetUniprotAccessions(ids,'uniprotTest.txt',limit=3)
    5. GetEnaAccessions(uniDict,'enaTest.txt',limit=3)
    6. GetAminoFiles(ids,'aminoTest')
    7. Dictfromfile('uniprotTest.txt',key=0,value=1,sep='\t',ignore='>')
    8. GetNucleotideFiles(enadict,'nucfilesTest')
    9. blastDB('nucfilesTest','blastDBtest')
    10. performBlast('aminoTest','CrossReferencesTaxId.txt','alignmentsTest',mode='local',blastdatabase='blastDBtest')
    11. GetBestLocalAlignments('alignmentsTest','nucfilesTest','descargas')
    12. TranslatedFiles('descargas','Translate','Positions.txt')
    13. ProteinAlignment('aminoTest/1eso.fasta','Translate/d1eso_.fasta','cmourra@lcg.unam.mx','1eso.needle')
    14. getPDBs('descargas','PDBS','fastasPDBs')

'''


import os
import re
import urllib
import wget
import Bio.Seq
from Bio import Entrez, SeqIO
import subprocess
import shutil


def getScopeId(path,newfilepath):
    """
    +++ Objective +++
        -- This function aims to extract the scope ids from the astral scope file
    +++ Parameters +++
        -- ::param:: path (string) -> Path to the astral scope file
        -- ::param:: newfilepath (string) -> Path to the file where the SCOPe - PDB association will be stored
    +++ Return +++
        -- scope_pdb (dictionary) -> A dictionary with the association scope - pdb
    """
    # Initialize a list where the scope ids will be stored and a dictionary with scope -> pdb
    scope = []
    scope_pdb = {}
    # Read file and get the scope ids
    with open(path,"r") as scopefile:
        print("+++++ Looking for SCOPE IDS +++++")
        for line in scopefile:
            line = line.replace("\n","")
            # Extract the first element in the list
            if line.startswith(">"):
                line = line.split(" ")
                scopeid = line[0].replace(">","")
                scope.append(scopeid)
        # Write the new file with the scope id - PDB id
        with open(newfilepath,"w") as newfile:
            newfile.write(">SCOPe ID\tPDB ID\n")
            for idscope in scope:
                scope_pdb[idscope] = idscope[1:5]
                newfile.write("{}\t{}\n".format(idscope,scope_pdb[idscope]))
    print("+++++ SCOPE IDs succesfully find it! +++++")
    return scope_pdb


def GetIdsfromfile(path):
    '''
    *************
    +++ Objective +++
        -- This function aims to get the PDB ids from the list file
    +++ Parameters +++
        -- ::param:: path (string) -> Path to the file that contains a list of the PDB ids
    +++ Return +++
        -- ids (list) -> A list with all the PDB ids
    *************
    '''
    # Initialize the list where the pdb ids will be stored
    ids = []
    # Open file to get the pdb ids
    with open(path,"r") as file:
        for line in file:
            # Ignore the header
            if line.startswith("#"):
                continue
            # Delete line break and append the pdb id to the list
            pdb = line.rstrip()
            ids.append(pdb)
    return ids


def GetAminoFiles(idsList,aminoFolder):
    '''
    *************
    +++ Objective +++
        -- This function aims to download the amino acid sequence of each PDB id into fasta format file
    +++ Parameters +++
        -- ::param:: idsList (list) -> List with all the PDB ids
        -- ::param:: aminoFolder (string) -> Path to the folder where all the fasta files will be stored
    +++ Return +++
        -- No return, is a procedure
    *************
    '''
    if os.path.isdir(aminoFolder) == False:
        os.mkdir(aminoFolder)
    print("++++ Downloading amino sequence fasta files ++++")
    # Loop through the ids list to get the pdbs
    for ids in idsList:
        # Generate the path for each protein
        filename = "{}.fasta".format(ids)
        # Generate the URL for each protein
        url = 'https://www.rcsb.org/fasta/entry/{}'.format(ids)
        # If the file is already in the folder, try with another protein
        if os.path.isfile("{}/{}".format(aminoFolder,filename)):
            print("The file {} already exists in this folder!".format(filename))
            continue
        # Download the amino sequence in the selected path with wget
        name = wget.download(url,out="{}/{}".format(aminoFolder,filename))
        print(name)


def GetUniprotAccessions(idsList,newfilepath,limit=100000):
    '''
    *************
    +++ Objective +++
        -- This function aims to generate a dictionary with the following structure {pdb_id : uniprot entry}
        -- and also create a file that will save the same data of this dictionary
    +++ Parameters (Required) +++
        -- ::param:: idsList (list) -> List with all the PDB ids
        -- ::param:: newfilepath (string) -> Path where the file with the pdb ids and uniprot entries will be stored
    +++ Parameters (Optional) +++
        -- ::param:: limit -> The number of proteins the program will search the ena id
    +++ Return +++
        -- uniprot (dictionary) -> A dictionary with the key -> PDB id, value -> uniprot entry
    *************
    '''
    # Initialize the dictionary that will save the association pdb - uniprot entry
    uniprot = {}
    # Open the new file to save the uniprot entries
    with open(newfilepath,"w") as newfile:
        # Write the header in the file
        newfile.write(">PDB ID\tUniprot ID\n")
        print("+++++ Looking for UNIPROT IDS +++++")
        # Loop through the pdb ids
        for count,ids in enumerate(idsList):
            # If the limit of proteins selected by the user is reached stop the program
            # and saved that number of proteins in the selected file
            if count == limit:
                break
            # Generate the pdb URL to get the uniprot Id
            link = 'https://www.rcsb.org/structure/{}'.format(ids)
            try:
                # Open the HTML with urllib
                with urllib.request.urlopen(link) as response:
                    html = response.read()
                    # Use a regular expression to find the uniprot entry
                    unipId = re.findall("www.uniprot.org/uniprot/\w+\" target=",str(html))[0]
                    # Delete extra characters from the string to get only the Uniprot id
                    unipId = unipId.replace("www.uniprot.org/uniprot/","").replace("\" target=","")
                    uniprot[ids] = unipId
                    # Write the pdb id and uniprot entry in a text file
                    textline = "{}\t{}\n".format(ids,uniprot[ids])
                    print(textline,end="")
                    newfile.write(textline)
            except:
                continue
    print("+++++ UNIPROT entries succesfully find it! +++++")
    return uniprot


def GetEnaAccessions(uniprotDict,newfilepath,limit=100000):
    '''
    *************
    +++ Objective +++
        -- This function aims to generate a dictionary with the following structure {pdb_id : ENA accession}
        -- and also create a file that will save the same data of this dictionary
    +++ Parameters (Required) +++
        -- ::param:: uniprotDict (dictionary) -> Dictionary with the PDB ids as keys and uniprot entries as values
        -- ::param:: newfilepath (string) -> Path where the file with pdb and ena ids will be stored
    +++ Parameters (Optionals) +++
        -- ::param:: limit -> The number of proteins the program will search the ena id
    +++ Return +++
        -- pdb_ena (dictionary) -> A dictionary with the key -> PDB id, value -> ENA accessions
    *************
    '''
    # Initialize with empty lists the dictionary to fill it with the association pdb - ena ids
    pdb_ena = {pdb:[] for pdb in uniprotDict}
    # Open the file where the ids will be saved
    with open(newfilepath,"w") as newfile:
        # Write header in the file
        newfile.write(">PDB IDs\tENA IDs\n")
        print("+++++ Looking for ENA IDS +++++")
        # Loop through the dictionary which has the association pdb - uniprot entry
        for count,pdb in enumerate(uniprotDict):
            try:
                # If the limit is reached, stop the program and return the dictionary
                if count == limit:
                    break
                # Generate the link according to each Uniprot Id
                link = "https://www.uniprot.org/uniprotkb/{}.xml".format(uniprotDict[pdb])
                with urllib.request.urlopen(link) as response:
                    html = response.read()
                    # Find the ENA ids with a regular expression
                    enaId = re.findall("<dbReference type=\"EMBL\" id=\"\w*",str(html))
                    # If no ENA ids were found search with another pdb id
                    if len(enaId) == 0:
                        continue
                    # Delete the part of the string that we are not going to use
                    # and append the ENA id according to their pdb
                    for index,ena in enumerate(enaId):
                        enaId[index] = ena.replace("<dbReference type=\"EMBL\" id=\"","")
                        pdb_ena[pdb].append(enaId[index])
                    # Get only the list without repeated ids
                    pdb_ena[pdb] = set(pdb_ena[pdb])
                    # Save the ids in a string separated by spaces and saved it in a tabular file
                    ENAs = " ".join(pdb_ena[pdb])
                    dataline = "{}\t{}\n".format(pdb,ENAs)
                    print(dataline,end="")
                    newfile.write(dataline)
            # If the script woldnt be able to do the connection try with another pdb
            except:
                continue
        # Delete the empty keys in the dictionary
        for pdb in uniprotDict:
            if len(uniprotDict[pdb]) == 0:
                del uniprotDict[pdb]
    print("+++++ ENA IDs succesfully find it! +++++")
    return pdb_ena


def Dictfromfile(filepath,key=0,value=1,sep='\t',ignore='>'):
    '''
    *************
    +++ Objective +++
        -- This function aims to return a dictionary with the data extracted from certain file
    +++ Parameters (Required) +++
        -- ::param:: outfile_path (string) -> Path to the outfile that
    +++ Parameters (Optional) +++
        -- ::param:: key (int) -> The number of the column that will be the keys (By default this argument take the value of 0)
        -- ::param:: value (int) -> The number of the column that will be the values (By default this argument take the value of 1)
        -- ::param:: sep (string) -> The caracter that will be the separator of the columns in the file
        -- ::param:: ignore (string) -> The character that will serves to recognize the header line
    +++ Return +++
        -- ::param:: new_dict (Dictionary) -> Dictionary with the selected key - value association
    *************
    '''
    # Initialize empty dictionary
    new_dict = {}
    # Open the file to extract the information
    with open(filepath,"r") as file:
        for line in file:
            # If the line is from the header ignore it
            if line.startswith(ignore):
                continue
            # Delete line breaks
            line = line.rstrip('\n')
            line = line.split(sep)
            # Get keys and values
            dict_key = line[key]
            dict_value = line[value]
            # Save the association key -> value
            new_dict[dict_key] = dict_value
    return new_dict


def AllIdsFile(PdbScopePath,PdbUniprotPath,PdbEnaPath,NewFilePath):
    '''
    *************
    +++ Objective +++
        -- This function aims to create a file with all the cross references ids
    +++ Parameters +++
        -- ::param:: outfile_path (string) -> Path to the outfile that
        -- ::param:: pdb_ena_dict (dictionary) -> Dictionary that contains the PDB ids as keys and the ENA accession
    +++ Return +++
        -- pdb_ena (dictionary) -> A dictionary with the key -> PDB id, value -> ENA accessions
    *************
    '''
    # Generate all the dictionaries with all the data we need for the file
    PdbScope = Dictfromfile(PdbScopePath,key=1,value=0)
    PdbUniprot = Dictfromfile(PdbUniprotPath)
    PdbEna = Dictfromfile(PdbEnaPath)
    alldict = {}
    print("+++++ Looking for ALL IDS +++++")
    # Open the file where the data will be saved
    with open(NewFilePath,"w") as newfile:
        # Write in the file the header line
        newfile.write(">PDB ID\tSCOPE ID\tUNIPROT ENTRY\tENA IDS\n")
        # Get each database id
        for pdb,scope in PdbScope.items():
            try:
                uniprot = PdbUniprot[pdb]
                ena = PdbEna[pdb]
                alldict[pdb] = [scope,uniprot,ena]
                # Save the information in one line
                dataline = "{}\t{}\t{}\t{}\n".format(pdb,scope,uniprot,ena)
                print(dataline,end="")
                newfile.write(dataline)
            # If the script isn't able to obtain the data try with another pdb
            except:
                continue
    print("+++++ ALL IDs succesfully saved! +++++")
    return alldict


def GetAminoLength(AminoFolderPath,AminoLengthPath):
    '''
    *************
    +++ Objective +++
        -- This function aims to generate a file that contains the PDB id and the length of the Sequence
    +++ Parameters (Required) +++
        -- ::param:: AminoFolderPath (dictionary): A dictionary with pdb ids as keys and Ena accessions as values
        -- ::param:: AminoLengthPath (string): The path to the folder where the nucleotide fastas will be downloaded
    +++ Return +++
        -- No return, it's a procedure
    *************
    '''
    # Get a list with the amino sequence filenames
    Aminofiles = os.listdir(AminoFolderPath)
    with open(AminoLengthPath,'w') as newfile:
        # Write an identifier for the columns
        header = '>PDB Id\tAmino Lenght\tExpected Nucleotide Length\n'
        newfile.write(header)
        for aminoFile in AminoFiles:
            # Generate the path for each amino sequence file
            pdb = aminoFile.replace('.fasta','')
            filepath = AminoFolderPath + '/' + aminoFile
            # Read file
            with open(filepath,'r') as aafile:
                # If the line is a header ignore it
                for line in aafile:
                    if line.startswith('>'):
                        continue
                    # Write the data in the new file
                    line = line.rstrip('\n')
                    seqLen = len(line)
                    nuccount = seqLen*3
                    text = '{}\t{}\t{}\n'.format(pdb,seqLen,nuccount)
                    newfile.write(text)


def GetNucleotideFiles(pdbEnaDict,nucFolder):
    '''
    *************
    +++ Objective +++
        -- This function aims to download the nucleotide fasta files from the given pdb ids
    +++ Parameters (Required) +++
        -- ::param:: pdbEnaDict (dictionary): A dictionary with pdb ids as keys and Ena accessions as values
        -- ::param:: nucFolder (string): The path to the folder where the nucleotide fastas will be downloaded
    +++ Return +++
        -- No return, it's a procedure
    +++ Usage +++
        -- GetNucleotideFiles(pdbEnaDict,'./Allfastas')
    *************
    '''
    if os.path.isdir(nucFolder) == False:
        os.mkdir(nucFolder)
    # Loop through each pdb
    for pdb in pdbEnaDict:
        # Generate the folder path
        folderpath = "{}/{}".format(nucFolder,pdb)
        # Check if the path already exists
        if os.path.isdir(folderpath) == False:
            # Create the folder that will contains the fastas
            os.mkdir(folderpath)
        for access in pdbEnaDict[pdb]:
            # Generate the URL for each ENA accession number
            url = "https://www.ebi.ac.uk/ena/browser/api/fasta/{}".format(access)
            filename = "{}/{}.fasta".format(folderpath,access)
            # If the file already exists in the folder try with another accession number
            if os.path.isfile(filename):
                print("The file {} already exists in this folder!".format(filename.replace(nucFolder + "/","")))
                continue
            # Open the web page to get the fasta file with RE and urllib
            with urllib.request.urlopen(url) as response:
                html = response.read()
                text = str(html).lstrip("b").replace("\'","").split("\\n")
                header = text[0]
                print(header)
                # If the file contains a whole genome ignore it and try with another fasta
                if "genome" in header:
                    continue
                # Write a new file to save the sequence
                with open(filename,"w") as newfile:
                    sequence = "".join(text[1:])
                    newfile.write(header + "\n" + sequence)
                    print("New file {} succesfully created!".format(filename.replace(nucFolder + "/","")))
    if os.path.isfile("{}/.DS_Store".format(nucFolder)):
        os.remove("{}/.DS_Store".format(nucFolder))

def blastDB(FastasFolder,dbFolderPath):
    '''
    *************
    +++ Objective +++
        -- This function aims to make a blast database for alignments
    +++ Parameters (Required) +++
        -- FastasFolder (String): Path to the nucleotide sequence folder
        -- dbFolderPath (String): Path to the folder where the databases will be saved
    +++ Return +++
        -- No return, it's a procedure
    +++ Usage +++
        -- blastDB("./AllFastas","./BlastDB")
    *************
    '''
    if os.path.isdir(dbFolderPath) == False:
        os.mkdir(dbFolderPath)
    # Get a list with the amino sequence filenames
    nucfolders = os.listdir(FastasFolder)
    # Loop through each nucleotide file and save the data in one fasta file
    for nucfolder in nucfolders:
        pdb = nucfolder
        allfastasPath = '{}/{}.fasta'.format(FastasFolder,pdb)
        databaseFolder = dbFolderPath + '/' + pdb
        databasePath = databaseFolder + '/' + pdb + 'db'
        databasefile = databaseFolder + '.fasta'
        fastafiles = os.listdir(FastasFolder + '/' + nucfolder)
        if os.path.isdir(databaseFolder) == False:
            os.mkdir(databaseFolder)
        for nucfile in fastafiles:
            filepath = FastasFolder + '/' + nucfolder + '/' + nucfile
            with open(filepath,'r') as file:
                with open(databasefile,'a') as newfile:
                    for line in file:
                        # Get the information of each fasta header
                        if line.startswith('>') or line.startswith("\""):
                            line = line.split('|')
                            # Write the sequence header
                            text = '>{}\n'.format(line[1])
                            newfile.write(text)
                            continue
                        newfile.write(line + '\n')
        # Once we have the fasta file, create the folder to store the database filenames
        if os.path.isdir(databaseFolder) == False:
            os.mkdir(databaseFolder)
        # Generate the bash command for make the blast database
        bashCommand = 'makeblastdb -in {} -out {} -parse_seqids -dbtype nucl'.format(databasefile,databasePath + '/' + pdb + 'db')
        print(bashCommand)
        # Run the bash command
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

def performBlast(AminoSequencePath,ReferencesFilePath,AlignmentFolderPath,mode='remote',blastdatabase=None):
    '''
    *************
    +++ Objective +++
        -- This function aims to perform blast process for a selected set of files
    +++ Parameters (Required) +++
        -- AminoSequencePath (String): Path to the amino sequence folder
        -- ReferencesFilePath (String): Path to the file with the CrossReferences
        -- AlignmentFolderPath (String): Path to the folder where the alignment files will be stored
    +++ Parameters (Optional) +++
        -- remote (int):
        -- blastdatabase (String):
    +++ Return +++
        -- No return, it's a procedure
    +++ Usage +++
        -- performBlast("./AminoFolder","./CrossReferencesTaxId.txt","./BlastAlignments")
        -- performBlast("./AminoFolder","./CrossReferencesTaxId.txt","./BlastAlignments","local","./BLASTDB")
    *************
    '''
    # Initialize a dictionary to save the association pdb - taxId
    pdb_taxid = {}
    # Generate de directory where Alignments where be save
    if os.path.isdir(AlignmentFolderPath) == False:
        os.mkdir(AlignmentFolderPath)
    # Open and read the CrossReferences file
    with open(ReferencesFilePath,'r') as References:
        for line in References:
            # Ignore the first header line
            if line.startswith('>'):
                continue
            # Delete line breaks and separate the line by tabs
            line = line.rstrip('\n')
            line = line.split('\t')
            try:
                # Get the pdbs and taxIds to fill the dictionary
                pdb = line[0]
                taxid = line[4]
                pdb_taxid[pdb] = taxid
            except:
                # If there isn't any taxId, try with another pdb
                continue
    # Get a list with the amino sequence filenames
    seqfiles = os.listdir(AminoSequencePath)
    # If the user gives a local database get the database for each protein
    if blastdatabase:
        nucdatabase = os.listdir(blastdatabase)
    # Loop through each file
    for aminofile in seqfiles:
        aminofilepath = AminoSequencePath + '/' + aminofile
        pdbfile = aminofile.replace('.fasta','')
        # If the pdb hasn't a taxId, try with another pdb
        if pdbfile not in pdb_taxid.keys():
            continue
        # Generate the alignment file path and perform the BLAST
        pdbfilepath = AlignmentFolderPath + '/' + pdbfile
        # Execute blast in commanda line depending on the option selected by user
        if mode == 'remote':
            bashCommand = "tblastn -query {} -db nt -remote -out {}.out -num_alignments 1 -outfmt 7".format(aminofilepath,pdbfilepath)
            # Run the bash comand line
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
        if mode == 'local':
            # get each database path for local alignment
            for nucdata in nucdatabase:
                if '.ipynb_checkpoints' in nucdata:
                    continue
                if '.fasta' in nucdata:
                    continue
                aminoFile = AminoSequencePath + '/' + '{}.fasta'.format(nucdata)
                BlastDataPath = blastdatabase + '/' + nucdata + '/' + nucdata + 'db' + '/' + nucdata + 'db'
                AlignPath = AlignmentFolderPath + '/' + nucdata
                bashCommand = "tblastn -query {} -db {} -out {}.out -outfmt 7".format(aminoFile,BlastDataPath,AlignPath)
                print(bashCommand)
                # Run the bash comand line
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()

def GetRemoteAlignments(remoteBlastFiles,BestNucFastas,mail):
    '''
    *************
    +++ Objective +++
        -- This function aims to download the nucleotide fastas corresponding for each pdb best remote alignment
    +++ Parameters (Required) +++
        -- ::param:: remoteBlastFiles (String): A path to the file which shows the higher percent of identity
        -- ::param:: BestNucFastas (String): The path to the folder where the nucleotide Fastas will be stored
    +++ Return +++
        -- No return, it's a procedure
    +++ Usage +++
        -- GetRemoteAlignments("./RemoteBlastAlignments","./NucleotideFastas",xxxxx@ccg.unam.mx)
    *************
    '''
    if os.path.isdir(BestNucFastas) == False:
        os.mkdir(BestNucFastas)
    # Get the remote blast alignment files
    filenames = os.listdir(remoteBlastFiles)
    Entrez.email = mail
    for filename in filenames:
        filepath = remoteBlastFiles + '/' + filename
        # Extract the NCBI accession from alignment file
        with open(filepath,'r') as BlastFile:
            for line in BlastFile:
                if line.startswith('#'):
                    continue
                line = line.rstrip('\n')
                elements = line.split('\t')
                print(elements)
                NCBIid = elements[1]
                # Get the nucleotide sequence with entrez
                handle = Entrez.efetch(db="nucleotide", id=NCBIid,rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                # Get the main elements for each nucleotide sequence
                sequence = record.seq
                ID = re.findall('ID:.+\n',str(record))[0].rstrip('\n').replace('ID: ','')
                NAME = re.findall('Name:.+\n',str(record))[0].rstrip('\n').replace('Name: ','')
                Description = re.findall('Description:.+\n',str(record))[0].rstrip('\n').replace('Description: ','')
                header = '>NCBI|{}|{} {}\n'.format(NAME,ID,Description)
                filepath = BestNucFastas + '/' + NCBIid + '.fasta'
                # Write the nucleotide fasta file
                with open(filepath,'w') as newfile:
                    newfile.write(header + str(sequence))

def GetBestLocalAlignments(AlignmentsPath,databasefiles,FinalNucFolder):
    '''
    *************
    +++ Objective +++
        -- This function aims to download the nucleotide fasta files for the best alignment from a local database
    +++ Parameters (Requiresd) +++
        -- ::param:: AlignmentsPath (String): A path to the folder which have the files with each protein alignment
        -- ::param:: databasefiles (String): The path to the folder where the nucleotide fastas from ENA database are stored
        -- ::param:: FinalNucFolder (String): Path to the file where the nucleotide fasta files will be saved
    +++ Return +++
        -- No return, it's a procedure
    +++ Usage +++
        -- TranslatedFiles('./NucleotideFastas','./TranslatedFiles','./PositionsData.txt')
    *************
    '''
    # Verify if the output folder already exist
    if os.path.isdir(FinalNucFolder) == False:
        os.mkdir(FinalNucFolder)
    # Get the filenames for each directory
    filepaths = os.listdir(AlignmentsPath)
    nucfiles = os.listdir(databasefiles)
    # Loop through each file name
    for filepath in filepaths:
        pdb = filepath.replace('.out','')
        # Open the BLAST alignment path for each protein
        with open(AlignmentsPath + '/' + filepath,'r') as alignfile:
            for line in alignfile:
                if line.startswith('#'):
                    continue
                else:
                    # Get the fields from each file
                    elements = line.split('\t')
                    fileid = elements[1]
                    identity = elements[2]
                    # Get the file with higher identity value
                    if float(identity) > 90:
                        filepath = databasefiles + '/' + pdb + '/' + fileid + '.fasta'
                        newfilepath = FinalNucFolder + '/' + 'd' + pdb + '_.fasta'
                        shutil.copyfile(filepath, newfilepath)
                        break
                    else:
                        continue


def TranslatedFiles(NucFolder,TranslateFolder,PositionsFile):
    '''
    *************
    +++ Objective +++
        -- This function aims to translate the nucleotide fasta files from the given pdb ids
    +++ Parameters (Requiresd) +++
        -- ::param:: NucFolder (String): A path to the folder which have the files with the higher percent of identity
        -- ::param:: TranslateFolder (String): The path to the folder where the translated Fastas will be stored
        -- ::param:: PositionsFile (String): Path to the file where the translated fasta information will be saved
    +++ Return +++
        -- No return, it's a procedure
    +++ Usage +++
        -- TranslatedFiles('./NucleotideFastas','./TranslatedFiles','./PositionsData.txt')
    *************
    '''
    # Generate a list with the folder names
    filenames = os.listdir(NucFolder)
    if os.path.isdir(TranslateFolder) == False:
        os.mkdir(TranslateFolder)
    # Open the file that will contains the positions of the file
    with open(PositionsFile,"w") as positionfile:
        # Loop through each pdb folder that contains the nucleotide files
        for filename in filenames:
            path = NucFolder + "/" + filename
            # Open the file and read the information within
            with open(path,"r") as file:
                # Initialize the sequence as an empty string
                sequence = ""
                # Read each line from the file
                for line in file:
                    # Ignore the first line and save the header
                    if line.startswith(">"):
                        header = line
                        continue
                    # Delete line breaks from the sequence
                    sequence += line.replace("\n","")
                aminoSeqs = []
                # Generate a dictionary that will contain the nucleotide sequence for each orf
                orfs = {"orf1":sequence,"orf2":sequence[1:],"orf3":sequence[2:],"orf4":str(Bio.Seq.Seq(sequence).reverse_complement()), "orf5":str(Bio.Seq.Seq(sequence).reverse_complement())[1:], "orf6":str(Bio.Seq.Seq(sequence).reverse_complement())[2:]}
                orfrange = {}
                # Translate each orf
                for orf in orfs:
                    try:
                        # Translate each orf by using a Biopython tool
                        protein = str(Bio.Seq.translate(orfs[orf])).split("*")
                    except:
                        continue
                    # Get the largest protein for each orf
                    value = max(protein,key=len)
                    # Get the positions of the largest peptide in the initial protein for each orf
                    listindex = protein.index(value)
                    # If the translated peptide is at the begining of the protein
                    if listindex == 0:
                        rangevalues = "{}-{}".format(1,len(value))
                    else:
                        # If the peptide is in another position sum the len of the other peptides
                        begin = 0
                        for prot in protein[:listindex]:
                            begin += len(prot)
                        # Write the range wich the longest peptide is
                        rangevalues = "{}-{}".format(begin,begin+len(value))
                    # Save the range of the positions and the translated sequence in a tupple
                    orfrange[orf] = (rangevalues,value)
                    # Save all the amino sequence corresponding for each orf in a list
                    aminoSeqs.append(value)
                print(orfrange)
                # Obtain which orf have the longest protein
                finalProt = max(aminoSeqs,key=len)
                #print(orfrange)
                for orf in orfrange:
                    # Get the information of the longest protein and save it in the file
                    if orfrange[orf][1] == finalProt:
                        text = "{},{},{},{}\n".format(orf,orfrange[orf][0],filename,NucFolder)
                        positionfile.write(text)
                # Save the amino sequence in the new folder for transalted files
                newfilepath = "{}/{}".format(TranslateFolder,filename)
                with open(newfilepath,"w") as newtranslate:
                    # Generate the text that will be write in the file
                    translatetext = header + str(finalProt)
                    newtranslate.write(translatetext)
                    # Indicate to the user that the file was succesfully translated
                    print("File {} succesfully translated!".format(str(newfilepath).replace("./TranslatedFolder/","")))


def ProteinAlignment(aSeq,bSeq,mail,out):
    '''
    *************
    +++ Objective +++
        -- This function aims to align 2 amino sequence files by using the program emboss_needle
    +++ Parameters (Requiresd) +++
        -- ::param:: aSeq (String): Path to the amino sequence file from PDB database
        -- ::param:: bSeq (String): Path to the translated nucleotide file
        -- ::param:: mail (String): The Email that will be use to do the alignment
        -- ::param:: out (String): Path where the output file will be saved
    +++ Return +++
        -- No return, it's a procedure
    +++ Usage +++
        -- ProteinAlignment("./1spd.fasta","./AK312116.fasta",xxxxxx@ccg.unam.mx,"./TEST/1spd.needle")
    *************
    '''
    # Run the bash comand line
    bashCommand = "python3 emboss_needle.py --asequence {} --bsequence {} --gapopen 10.0 --gapext 0.5 --outfile {}  --email {}".format(aSeq,bSeq,out,mail)
    print(bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

def getPDBs(InputFolder,OutputPath,FastasPath):
    '''
    *************
    +++ Objective +++
        -- This function aims to download the pdb files from a set of proteins
    +++ Parameters (Requiresd) +++
        -- ::param:: InputFolder (String): Path to the nucleotide Fastas Folder
        -- ::param:: OutputPath (String): Path to the Folder where the files will be saved
        -- ::param:: FastasPath (String): The Folder where the Fasta files will be copied
    +++ Return +++
        -- No return, it's a procedure
    +++ Usage +++
        -- getPDBs("./descargas","./input_pdbs_calc_consensus","./input_fastas_calc_consensus")
    *************
    '''
    Files = os.listdir(InputFolder)
    # Verify if the paths already exists, if not make it
    if os.path.isdir(OutputPath) == False:
        os.mkdir(OutputPath)
    if os.path.isdir(FastasPath) == False:
        os.mkdir(FastasPath)
    Proteins = []
    # Loop through each file and copy it to the other folder
    for line in Files:
        FastaFile = "{}/{}".format(InputFolder,line)
        CopyCommand = "cp {} {}".format(FastaFile,FastasPath)
        os.popen(CopyCommand)
        # Obtener una lista de proteinas que queremos descargar
        line = line.replace("d","").replace(".fasta","")
        # Get the pdb id from the file name and
        Prot = line[0:4]
        Proteins.append(Prot)
    # Get each protein from the List and download it from the database with wget
    for protein in Proteins:
        try:
            # Generamos el url para poder descargar los PDBs
            url = "http://www.pdb.org/pdb/files/{}.pdb.gz".format(protein)
            wget.download(url, OutputPath)
        except:
            continue


def SearchPDBlist(pdblistpath,email,nucFolder,PdbFiles,mode='remote',crossdatabase=None):
    '''
    *************
    +++ Objective +++
        -- This function aims to download the pdb files from a set of proteins
    +++ Parameters (Required) +++
        -- ::param:: pdblistpath (String): Path to the pdb list file
        -- ::param:: email (String): Email as an identifier to perform the BLAST
        -- ::param:: nucFolder (String): Path to the folder where all the fasta files will be download
        -- ::param:: PdbFiles (String): Path to the folder where the pdb files will be saved
    +++ Parameters (Optional) +++
        -- ::param:: mode (String): Select if the BLAST will perform the aligments in the -remote server or make a local database
        -- ::param:: crossdatabase (String): Path to the nucleotide Fastas Folder
    +++ Return +++
        -- No return, it's a procedure
    +++ Usage +++
        -- SearchPDBlist("./pdblist.txt","cmourra@lcg.unam.mx","./allFastas","./input_pdbs_calc_consensus")
        -- SearchPDBlist("./pdblist.txt","cmourra@lcg.unam.mx","./allFastas","./input_pdbs_calc_consensus","local",'./CrossReferencesTaxId.txt')
    *************
    '''
    # Get all the PDB ids in a list object from the list file
    with open(pdblistpath,"r") as pdblistfile:
        allfilepdbs = []
        for line in pdblistfile:
            if line.startswith('>'):
                continue
            pdb = line.rstrip('\n')
            allfilepdbs.append(pdb)
    # In the default case permorm a remote alignment
    if mode == 'remote':
        if os.path.isdir('./AminoPDBlist') == False:
            os.mkdir('./AminoPDBlist')
        GetAminoFiles(allfilepdbs,'./AminoPDBlist')
        performBlast('./AminoPDBlist','./CrossReferencesTaxId.txt','./PDBlistAlignments')
        GetRemoteAlignments('./PDBlistAlignments','./descargas',email)
        getPDBs('./descargas',PdbFiles,'./input_fastas_calc_consensus')
    # If the user specifies run a local alignment by making databases
    if mode == 'local':
        crosspdb = []
        # Get all the pdbs present in the cross reference databases
        with open(crossdatabase,'r') as crossfile:
            for line in crossfile:
                if line.startswith('>'):
                    continue
                elements = line.split('\t')
                pdb = elements[0]
                crosspdb.append(pdb)
        notincross = []
        incross = []
        # Classify the pdbs if they are or not in the database file
        for filepdb in allfilepdbs:
            if filepdb not in crosspdb:
                notincross.append(filepdb)
            else:
                incross.append(filepdb)
        with open(crossdatabase,'r') as crossfile:
            databaseinfo = {}
            incrossdict = {}
            for line in crossfile:
                if line.startswith('>'):
                    continue
                line = line.split('\t')
                pdb = line[0]
                enaids = line[3].split(' ')
                databaseinfo[pdb] = enaids
            for pdbkey in databaseinfo:
                if pdbkey in incross:
                    incrossdict[pdbkey] = databaseinfo[pdbkey]
                else:
                    continue
        GetNucleotideFiles(incrossdict,nucFolder)
        GetAminoFiles(allfilepdbs,'./AminoPDBlist')
        performBlast('./AminoPDBlist',crossdatabase,'./PDBlistAlignments',"local","./BLASTDB")
        getPDBs('./descargas',PdbFiles,'./input_fastas_calc_consensus')
        with open(crossdatabase,'a') as crossfile:
            pdb_ena_not = {}
            for notpdb in notincross:
                pdblink = 'https://www.rcsb.org/structure/{}'.format(notpdb)
                with urllib.request.urlopen(pdblink) as response:
                    unipId = re.findall("www.uniprot.org/uniprot/\w+\" target=",str(html))[0]
                    # Delete extra characters from the string to get only the Uniprot id
                    unipId = unipId.replace("www.uniprot.org/uniprot/","").replace("\" target=","")
                uniprotlink = "https://www.uniprot.org/uniprotkb/{}.xml".format(notpdb)
                with urllib.request.urlopen(uniprotlink) as response:
                    enas = []
                    html = response.read()
                    enaId = re.findall("<dbReference type=\"EMBL\" id=\"\w*",str(html))
                    # If no ENA ids were found search with another pdb id
                    if len(enaId) == 0:
                        continue
                    # Delete the part of the string that we are not going to use
                    # and append the ENA id according to their pdb
                    for index,ena in enumerate(enaId):
                        enaId[index] = ena.replace("<dbReference type=\"EMBL\" id=\"","")
                        enas.append(enaId[index])
                    enas = set(enas)
                    pdb_ena_not[not_pdb] = enas
                    # Save the ids in a string separated by spaces and saved it in a tabular file
                    ENAs = " ".join(enas)
                text = "\t{}\t{}\t{}\t\n".format(notpdb,unipId,ENAs)
                crossfile.write(text)
            GetNucleotideFiles(pdb_ena_not,nucFolder)
            performBlast('./AminoPDBlist',crossdatabase,'./PDBlistAlignments',"local","./BLASTDB")
            GetBestLocalAlignments('./PDBlistAlignments',nucFolder,'./descargas')
            getPDBs('./descargas',PdbFiles,'./input_fastas_calc_consensus')
