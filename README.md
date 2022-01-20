
<p align="center">
  <img src="https://scontent.fcvj5-1.fna.fbcdn.net/v/t1.6435-9/126787074_126739485894201_6763364814445598779_n.png?_nc_cat=104&ccb=1-5&_nc_sid=09cbfe&_nc_ohc=jDA0jgQszkEAX-tym7C&_nc_ht=scontent.fcvj5-1.fna&oh=00_AT_JB8K0B1N_lGqaH8ondV34lJdfUJw21B_JESBESYrL5w&oe=621065E8" width="200px" height="200px"/></p>

# FariasLab

This repository contains all the scripts and files necessary to find cotranslational sites that allow folding intermediates to form .

## Order of execution of the programs

1. ProtList.py
2. aa_to_nt.py
3. GetPDBs.py
4. MinMax.py
5. AutomaticSaknovich.sh
6. CoTransFiles.sh
7. CoTransSites.py

# Folders

- ## Extraction

  This folder aims to extract from the SCOPe database, the proteins that we will use to run the analysis.

  ### **Files:**
  
---------------------

  - [**ScopDatabaseFile.txt**](https://github.com/MichMourra/FariasLab/blob/main/Extraction/ScopDatabaseFile.txt) :page_facing_up: : This file contains all the proteins information like: Protein identifier, family number and polypeptide sequence.

<p align="center">
  <img src="https://scop.mrc-lmb.cam.ac.uk/static/media/scop_logo.358e74fa.png"/></p>


---------------------

  - [**ProtList.py**](https://github.com/MichMourra/FariasLab/blob/main/Extraction/ProtList.py) :white_check_mark: : This script has the objective to generate a list file that will contains the protein id and family number taking as arguments a fasta file with the protein data from SCOPe.

    **Arguments**

    - **inputfile**: The file path to the SCOPe database file
    - **outputfile**: The path to the new list file

    **How to run**
    ```diff
    + python3 ProtList.py -i ./ScopDatabaseFile.txt -o ScopeNewList.txt
    ```
---------------------   
 
  - [**ScopeNewList.txt**](https://github.com/MichMourra/FariasLab/blob/main/Extraction/ScopeNewList.txt) :page_facing_up: : Output file from ProtList.py , this file is a list that contains the protein id and the family id. This file serves as the input file to the aa_to_nt.py script

---------------------

- ## Preparation

  The files in this folder aims to download the fasta and pdb files necessary for execute the programs to quantify the translational properties for each protein : **MinMax.py** and **AutomaticSaknovich.sh**

  ### **Files**
---------------------

  - [**aa_to_nt.py**](https://github.com/MichMourra/FariasLab/blob/main/Preparation/aa_to_nt.py) :exclamation: : This script has as purpose to download the entire sequence of peptides of interest proteins in fasta format, taking as input file a list that contains the protein id and the family number. This data was previously captured from SCOPe database by the SCOPeList.py script.

   ```diff
   - This program has been officially deprecated due to changes in the html base of the scope database portal. 
   - However we keep the script for a possible future update
   ```
   

   **Arguments**

   - **inputPath**: Path to read the input file (list)
   - **outputPath**: Path to place the folders with all the downloaded data and error lists
   - **listNameFile**: Name of the file with the name of all the data to download
   - **cut**: If it is equal to "TRUE", only the part of the ENA sequence that corresponds to the Scope sequence is taken, otherwise (FALSE), the complete sequence provided by it will be taken.

   **How to run** :exclamation:
   
   ```diff
   - python3 aa_to_nt.py --inputPath ../Extraction --outputPath . --listNameFile ScopeNewList.txt --cut FALSE
   ```
    
---------------------
   
  - [**DownloadFastas.py**](https://github.com/MichMourra/FariasLab/blob/main/Preparation/DownloadFastas.py) :white_check_mark: : This script is able to download the fasta files using the EMBL ids extracted from the ENA file. With this ids we can get the nucleotide sequence using the entrez and then we save this data into a fasta file of the protein.

**Arguments**

 - **InputFilePath**: Path to the file that contains the ENA proteins data.
 - **OutputFolder**: Path to the folder that will contains the resulting files.
 - **Email**: Path to the folder that contains the elongation profile files.

**How to run**
```diff
+ ./DownoloadFastas.py -I ./ENA.txt -O "./SCOPeFastas/" -E cmourra@lcg.unam.mx
```
<p align="center">
  <img src="https://scop.mrc-lmb.cam.ac.uk/static/media/scop_logo.358e74fa.png" width="250px" height="100px"/></p>
  
<p align="center">
  <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/e2/Uniprot.gif/245px-Uniprot.gif" width="250px" height="100px"/></p>
  
<p align="center">
  <img src="https://upload.wikimedia.org/wikipedia/en/thumb/b/b1/EMBL_logo.svg/240px-EMBL_logo.svg.png" width="250px" height="100px"/></p>

---------------------

  - [**GetPDBs.py**](https://github.com/MichMourra/FariasLab/blob/main/Preparation/GetPDBs.py) :white_check_mark: : This script is able to download the PDB files from uniprot database by establishing a programmatic access. Takes as argument the folder path that contains all the fasta files resulting from aa-to-nt.py.

    **Arguments**

    - **InputFolder**: Path to the folder that contains the nucleotide fastas
    - **OutputPath**: Path to the folder that will contains the resulting files
    - **FastasPath**: Path to the folder that will contains the fasta files to necesarry to run AutomaticSaknovich.sh

    **How to run**
    
    ```diff
    + python3 ./GetPDBs.py -i ./descargas -o ../cg_cotrans/input_pdbs_calc_consensus/ -f ../cg_cotrans/input_fastas_calc_consensus/
    ```

    **Output**

    - **descargas**: This folder will contain the nucleotide fastas

<p align="center">
  <img src="https://cdn.rcsb.org/rcsb-pdb/v2/common/images/rcsb_logo.png" width="250px" height="100px"/></p>

---------------------

  - [**chromedriver**](https://github.com/MichMourra/FariasLab/blob/main/Preparation/chromedriver) :x: : This is the google chrome driver used by selenium to do the fastas search. If there is an error with the driver this could be caused by two things:

    1. This chromedriver is for macOS system and you're trying to run the script in linux or windows
    2. Selenium needs a new version of the chromedriver.

    In either case, you can download the chromedriver you need from the following URL:

    - https://chromedriver.chromium.org/downloads


   ```diff
   - This driver is used by the aa-to-nt.py script. 
   - But since that program is currently no longer used, the driver is unnecessary in the new DownloadFastas.py program. 
   ```


<p align="center">
  <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/a/a5/Google_Chrome_icon_%28September_2014%29.svg/245px-Google_Chrome_icon_%28September_2014%29.svg.png" width="150px" height="150px"/></p>

---------------------

  - [**descargas**](https://github.com/MichMourra/FariasLab/tree/main/Preparation/descargas) :file_folder: : 

    This folder contains the fasta files resulting from aa_to_nt.py program
    
---------------------

- ## Execution

  The files in this folder aims to execute the MinMax.py and AutomaticSaknovich.sh this in order to obtain the translational properties of each protein, properties like: rare codon frequency, native contact in the proteins, stabilizing energies of the protein and elongation rate of the protein.

  ### Files
  
---------------------

  - [**MinMax.py**](https://github.com/MichMourra/FariasLab/blob/main/Execution/MinMax.py) :white_check_mark: : This script has the task of evaluating the relative usage frequencies of the synonymous codons used to encode a protein sequence of interest and compares these results to a rigorous null model. MinMax.py takes as argument the folder path that contains all the fasta files resulting from aa-to-nt.py. 

    **Arguments**

    - **inputFolderPath**: Path to read the folder with files to evaluate
    - **outputPath**: Path to place the folder with all the generated figures.

    **How to run**
    ```diff
    + python3 MinMax.py --inputFolderPath ../Preparation/descargas --outputPath .
    ```
    **Output**

    - **figuras**: Thi folder will contains the resulting files by execute the MinMax.py program

---------------------

  - [**figuras**](https://github.com/MichMourra/FariasLab/tree/main/Execution/figuras) :file_folder: : This folder contains the csv files that has the rare codon frequency for each site of the protein and also has the graphics to see more easily the areas enriched with rare codons.

---------------------

- ## cg_cotrans

  This folder contains all the scripts and data necessary for execute the AutomaticSaknovich.sh program.

  ### Files
  
  ---------------------

  - [**Automatic_Saknovich.sh**](https://github.com/MichMourra/FariasLab/blob/main/cg_cotrans/Automatic_Saknovich.sh) :white_check_mark: : This bash script has the task of calling another python scripts to calculate three properties including the native contacts within a protein, the stabilizing energies per codon and the elongation rate during translation in the given protein sequences. This program take as first argument a folder path that contains all the fasta files resulting from aa-to-nt.py and as second argument a folder path with all the PDB files resulting from getPDBs.py.

    **How to run**
    ```diff
    + ./AutomaticSaknovich.sh
    ```
    **Input folders**

    - [input_fastas_calc_consensus](https://github.com/MichMourra/FariasLab/tree/main/cg_cotrans/input_fastas_calc_consensus) :file_folder: : This folder contains the fastas with the nucleotide sequence for each protein, this are the same fastas that we get by using the aa_to_nt.py program.
    - [input_pdbs_calc_consensus](https://github.com/MichMourra/FariasLab/tree/main/cg_cotrans/input_pdbs_calc_consensus) :file_folder: : This folder contains the pdbs of each protein, this are the same pdbs that we get by using the GetPDBs.py program.

    **Scripts executed by AutomaticSaknovich**

    1. [**calc_consensus_contacts.py**](https://github.com/MichMourra/FariasLab/blob/main/cg_cotrans/calc_consensus_contacts.py) :white_check_mark: : A script that has the purpose of calculate the native contacts formed when the proteins exits from the ribosome.
    2. [**calc_absolute_energies.py**](https://github.com/MichMourra/FariasLab/blob/main/cg_cotrans/calc_absolute_energies.py) :white_check_mark: : A script which calculates the energies for stabilizing the protein when this exits the ribosome.
    3. [**calc_elongation.py**](https://github.com/MichMourra/FariasLab/blob/main/cg_cotrans/calc_elongation.py) :white_check_mark: : This python script has the purpose of calculate the elongation rate of the protein during translation.

    **Output folder**

    - [**output_input_ALL**](https://github.com/MichMourra/FariasLab/tree/main/cg_cotrans/output_input_ALL) :file_folder: : This folder contains all the files produced by calculations made by AutomaticSaknovich.sh program.

---------------------

- ## Discovery

  This folder contains the files necessary to find the cotranslational site of each protein by using the data produced by the previous programs.

  ### Files
  
---------------------

  - [**CoTransFiles.sh**](https://github.com/MichMourra/FariasLab/blob/main/Discovery/CoTransFiles.sh) :white_check_mark: : This bash script has of purpose that the elongation rate and MinMax files are of the same proteins and are in the same proportion. This program also create a folder where the reviewed files will be located.

    **Arguments**

    - **ElongationDir**: The directory that contains all the elongation files generated by the AutomaticSaknovich.sh program.
    - **MinMaxDir**: The directory that contains all the MinMax files generated by the MinMax.py program.
    - **OutDir**: Path to the new directory that will contains the new directories for the input files and the resulting files
    - **PythonScript**: The path where the CoTransSites.py script is located

    **How to run**
    ```diff
    + ./CoTransFiles.sh ./output_input_ALL ./figuras ./CoTranslate ./CoTransSites.py
    ```
---------------------

  - [**CoTransSites.py**](https://github.com/MichMourra/FariasLab/blob/main/Discovery/CoTransSites.py) :white_check_mark: : This python script has the purpose of finding the cotranslational sites by using the elongation and MinMax files obtained by running the MinMax.py and AutomaticSaknovich.sh programs.

    **Arguments**

    - **ElongationFolder**: Path to the folder that contains the elongation profile files
    - **MinMaxFolder**: Path to the folder that contains the elongation profile files
    - **OutputFolder**: Path to the folder that will contains the resulting files

    **How to run**
    ```diff
    + python3 ./CoTransSites.py -E ./output_input_ALL -M ./figuras -O ./CoTrans
    ```
    **Output**

    - [**CoTranslate**](https://github.com/MichMourra/FariasLab/tree/main/Discovery/CoTranslate) :file_folder: : This folder contains the folders with the reviewed input files, the text file with each cotranslational site found it by the program and a durectory for each protein within graphics of ech cotranslational site.
