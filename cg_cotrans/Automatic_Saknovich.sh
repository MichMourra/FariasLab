#!/bin/bash
# El siguiente comando puede ser utilizado para correr este programa y que no si interrumpa si se cierra sesion:
# nohup ./Automatic_Saknovich.sh &
# nohup Automatic_Saknovich.sh &


# Este programa trabaja con varios directorios. Uno que contiene todos los fastas a calcular (input_fastas_calc_consensus), otro que contiene todos los pdbs correspondientes (input_pdbs_calc_consensus) y un tercer directorio en donde estan todos los outputs de los tres programas (output_input_ALL).

# Para hacer los calculos de Saknovich se deben usar tres programas en este orden:
    #calc_consensus_contacts.py --> recibe dos entradas, un archivo fasta de nucleotidos y su archivo pdb correspondiente.
    #calc_absolute_energies.py --> recibe el output de calc_consensus_contacts.py
    #calc_elongation.py --> recibe el output de calc_absolute_energies.py

echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Ejecucion de todas las secuencias del directorio input_fastas_calc_consensus"
echo "Programas: calc_consensus_contacts.py -> input_fastas_calc_consensus.py -> calc_elongation.py"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "\n"


# Recorrer todos los fastas dentro del directorio input de fastas
for i in $(cd input_fastas_calc_consensus && ls):
do

    # Extraccion del nombre del FASTA. Eg. d1a3wa2.fasta ----> d1a3wa2
    filename=$(basename -- "$i")
    fasta="${filename%.*}"
    # Extraccion del nombre del archivo pdb. Eg d1a3wa2 ----> 1a3w. AQUI SE ASUME QUE EL NOMBRE DEL FASTA Y DEL PDB SI CORRESPONDEN (EJEMPLO: d1a3wa2 y 1a3w).
    pdb=$(echo $fasta| cut -c2-5)

    # Estos dos datos extraidos se repiten todo el tiempo en las lineas de comandos de los tres programas de Saknovich, por lo que
    # solo se sustituyen.

    echo "Archivo: $fasta"
    echo "PDB: $pdb"
    echo "\n"

    # Esta es la linea que hace todo lo de Saknovich automatica, asi que solo habria que asegurarse de que existan los tres directorios con ese nombre y con la informacion que necesita cada uno dentro para que este pequenio script funcione.
    (python3 calc_consensus_contacts.py --fasta ./input_fastas_calc_consensus/$fasta.fasta --pdbids $pdb --structure-path ./input_pdbs_calc_consensus/ --output ./output_input_ALL/$fasta.dat && python3 calc_absolute_energies.py --polymer ./output_input_ALL/$fasta.dat && python3 calc_elongation.py --gene $fasta --polymer ./output_input_ALL/$fasta._abs_energies.dat --output ./output_input_ALL) &

done
##
wait
