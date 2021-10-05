#!/bin/bash


####################################################################################################
####################################################################################################
###																																					                     ###
###														           How to execute																           ###
###																																					                     ###
###       ./CoTransFiles.sh ./output_input_ALL ./figuras ./CoTranslate ./CoTransSites.py				 ###
###																																					                     ###
####################################################################################################
####################################################################################################


# Extract the directories where the information is stored
ElongationDir=$1
MinMaxDir=$2
OutDir=$3
PythonScript=$4

# Get this script path
Path=$(pwd)

# Eliminate spaces in the MinMax folder file names
cd $MinMaxDir
for file in *; do mv "$file" `echo $file | tr ' ' '_'`; done
cd $Path

# Eliminate spaces in the Elongation folder file names
cd $ElongationDir
for file in *; do mv "$file" `echo $file | tr ' ' '_'`; done
cd $Path

# Get file paths
ElongationFiles=$ElongationDir/*elongation_profile.dat
MinMaxFiles=$MinMaxDir/*win17.csv

#Generate files directories
ElongationFolder=$OutDir/ElongationFolder
MinMaxFolder=$OutDir/MinMaxFolder
mkdir $ElongationFolder
mkdir $MinMaxFolder

# Get Elongation and MinMax file proteins into a list
ElongationList=$(ls $ElongationFiles | grep -o '/.*__' | sed 's/\/.*\///g' | sed 's/_.*$//g' | uniq)
MinMaxList=$(ls $MinMaxFiles | grep -o 'MinMax_.*_' | sed 's/MinMax_//g' | sed 's/_.*_$//g' | uniq)

# Loop through the list elements and
for EloProt in $ElongationList
do
	for MinMaxProt in $MinMaxList
	do
		# Get which proteins are the same in both sets
		if [[ $EloProt == $MinMaxProt ]];
		then
			for MinMax in $(ls $MinMaxFiles | sed 's/\.\/.*\///g')
			do
				# Get the MinMax protein file
				if [[ $MinMax == *"$MinMaxProt"* ]];
				then
					# Copy the file into the new directory
					cp $MinMaxDir/$MinMax $MinMaxFolder
				fi
			done

			for ElongationData in $(ls $ElongationFiles | sed 's/\.\/.*\///g')
			do
				# Get the Elongation profile file
				if [[ $ElongationData == *"$EloProt"* ]];
				then
					# Copy the file into the new directory
					cp $ElongationDir/$ElongationData $ElongationFolder
				fi
			done
		fi
	done
done


python3 $PythonScript -E $ElongationDir -M $MinMaxDir -O $OutDir
