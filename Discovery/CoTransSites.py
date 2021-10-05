

#··················································································#
#··················································································#
#                                  How to execute                                  #
#    python3 ./CoTransSites.py -E ./output_input_ALL -M ./figuras -O ./CoTrans     #
#··················································································#
#··················································································#


#··················································································#
#··················································································#
#                                      Modules                                     #
#··················································································#
#··················································································#


import os
import matplotlib.pyplot as plt
import argparse

#··················································································#
#··················································································#
#                                     Functions                                    #
#··················································································#
#··················································································#


def obtener_pendiente(x,y):
    limit = len(x)
    valores = []
    i1 = 0
    i2 = 1
    while(i2 < limit):
        try:
            Pend = (y[i2] - y[i1])/(x[i2] - x[i1])
            rango = "{}-{}".format(x[i1],x[i2])
            values = (rango,Pend)
            valores.append(values)
        except:
            continue
        i1 = i1 + 1
        i2 = i2 + 1
    return valores

def getChange(elements):
    Positions = []
    Starts = []
    for values in elements:
        position = values[0].split("-")
        if (abs(float(values[1])) > 1):
            start = float(position[0]) - (float(position[0]) % 20)

            if start in Starts:
                continue

            Starts.append(start)
            end = start + 20

            Positions.append((start,end))
        else:
            continue
    return Positions


#··················································································#
#··················································································#
#                                     Arguments                                    #
#··················································································#
#··················································································#

# Arguments for request the minmax file paths, elongation profile paths and output path
parser = argparse.ArgumentParser()
parser.add_argument("-E", "--ElongationFolder", help="Path to the folder that contains the elongation profile files")
parser.add_argument("-M", "--MinMaxFolder", help="Path to the folder that contains the elongation profile files")
parser.add_argument("-O", "--OutputFolder", help="Path to the folder that will contains the resulting files")
args = parser.parse_args()

# Get the content of the arguments
PathOutput = args.ElongationFolder
PathMinMax = args.MinMaxFolder
PathCoT = args.OutputFolder
PathCoTrans = "{}/CoTransSites.txt".format(PathCoT)



#··················································································#
#··················································································#
#                                     Main Code                                    #
#··················································································#
#··················································································#

# Get the filepaths in the directories
Output = os.listdir(PathOutput)
MinMaxO = os.listdir(PathMinMax)

# List that contains the files paths
Elongation = []
MinMax = []

# Get the list of only the elongation profile paths
for FileNameOutput in Output:
    if "_elongation_profile.dat" in FileNameOutput:
        Elongation.append(FileNameOutput)
        continue
    else:
        continue

# Get the list with the minmax files with window size 17
for FileNameMinMax in MinMaxO:
    if "win17.csv" in FileNameMinMax:
        MinMax.append(FileNameMinMax)
        continue
    else:
        continue

# Open the file that will contain all the cotranslational sites
CoTransFile = open(PathCoTrans,"w")

Count = 0
Sites = 0
for FileElo in Elongation:
    for FileMin in MinMax:

        ElongationProtein = FileElo.replace("_elongation_profile.dat","")
        MinMaxProtein = FileMin.replace("MinMax_","")
        MinMaxProtein = MinMaxProtein[0:7]


        if (ElongationProtein == MinMaxProtein):

            print("****************")
            print("Elongation profile protein: {}".format(ElongationProtein))
            print("****************")
            print("++++++++++++++++")
            print("MinMax protein: {}".format(MinMaxProtein))
            print("++++++++++++++++")

            ElongationPath = "{}/{}".format(PathOutput, FileElo)
            MinMaxPath = "{}/{}".format(PathMinMax,FileMin)

            # Open files
            MinMaxFile = open(MinMaxPath,"r")
            ElongationFile = open(ElongationPath,"r")

            # Elongation data list
            xElongation = []
            yElongation = []

            xMinMax = []
            yMinMax = []

            # MinMax file data list
            MinMaxList = []
            RRT = []

            # Get column values of the elongation file
            i = 0
            for line in ElongationFile:
                if i == 0:
                    i = i + 1
                    continue
                line = line.split("\t")
                xElongation.append(float(line[0]))
                yElongation.append(float(line[3]))
                i = i + 1

            # Get pending values of the elongation profile .dat file
            Pending = obtener_pendiente(xElongation,yElongation)
            # Get pending values higher than 1
            Change = getChange(Pending)

            # Get MinMax data
            flag = 0
            for line in MinMaxFile:
                if flag == 0:
                    flag = 1
                    continue
                line = line.split(",")
                # Get minmax values
                try:
                    value = (float(line[0]),float(line[3]))
                    xMinMax.append(value[0])
                    yMinMax.append(value[1])
                    MinMaxList.append(value)
                    RRT.append(float(line[4]))
                except:
                    continue

            MinusValue = []
            for element in MinMaxList:
                if element[1] < 0:
                    minus = element[0] % 10
                    start = element[0] - minus - 10
                    end = start + 20
                    Range = (start,end)
                    MinusValue.append(Range)

            MinusValue = set(MinusValue)
            if ((len(MinusValue) > 0) and (len(Change) > 0)):
                CoTransFile.write(">" + FileElo + "\n")
                ProteinFolderPath = "{}/{}".format(PathCoT, ElongationProtein)

                try:
                    os.mkdir(ProteinFolderPath)
                except:
                    print(">Protein Folder already exist!")

                for MinValue in MinusValue:
                    for Elongated in Change:
                        if ((MinValue[0] >= Elongated[0]) and (MinValue[0] <= Elongated[1])):
                            Sites += 1
                            print("{} cotranslational sites found".format(Sites))
                            print("Site: {} , {}".format(Elongated[0], MinValue[1]))
                            NewData = "{},{}\n".format(Elongated[0],MinValue[1])
                            CoTransFile.write(NewData)
                            plot1 = plt.figure(1,figsize=(15,8))
                            plt.plot(xElongation,yElongation,color="Blue")
                            plt.axvline(x=Elongated[0], color='green', linestyle ="--")
                            plt.axvline(x=MinValue[1], color='green', linestyle ="--")
                            plt.title("Elongation Profile " + ElongationProtein)
                            plt.xlabel("L (residue)")
                            plt.ylabel("F_L")
                            plt.grid()
                            ImgPath1 = "{}/ElongationProfile_{}Site_{}.png".format(ProteinFolderPath, ElongationProtein, Sites)
                            plot1.savefig(ImgPath1, dpi=plot1.dpi)
                            plt.clf()


                            plot2 = plt.figure(2,figsize=(15,8))
                            plt.plot(xMinMax,yMinMax,color="Blue")
                            plt.axhline(y=0, color='r', linestyle ="--")
                            plt.axvline(x=Elongated[0], color='green', linestyle ="--")
                            plt.axvline(x=MinValue[1], color='green', linestyle ="--")
                            plt.title("MinMax " + MinMaxProtein)
                            plt.xlabel("Center of codon window")
                            plt.ylabel("%MinMax")
                            plt.grid()
                            ImgPath2 = "{}/MinMax_{}Site_{}.png".format(ProteinFolderPath, MinMaxProtein, Sites)
                            plot2.savefig(ImgPath2, dpi=plot2.dpi)
                            plt.clf()

                            continue

                        if ((MinValue[1] >= Elongated[0]) and (MinValue[1] <= Elongated[1])):
                            Sites += 1
                            print("{} cotranslational sites found".format(Sites))
                            NewData = "{},{}\n".format(MinValue[0],Elongated[1])
                            print("Site: {} , {}".format(MinValue[0], Elongated[1]))
                            CoTransFile.write(NewData)
                            plot1 = plt.figure(1,figsize=(15,8))
                            plt.plot(xElongation,yElongation,color="Blue")
                            plt.axvline(x=MinValue[0], color='purple', linestyle ="--")
                            plt.axvline(x=Elongated[1], color='purple', linestyle ="--")
                            plt.title("Elongation Profile " + ElongationProtein)
                            plt.xlabel("L (residue)")
                            plt.ylabel("F_L")
                            plt.grid()
                            ImgPath1 = "{}/ElongationProfile_{}Site_{}.png".format(ProteinFolderPath, ElongationProtein, Sites)
                            plot1.savefig(ImgPath1, dpi=plot1.dpi)
                            plt.clf()


                            plot2 = plt.figure(2,figsize=(15,8))
                            plt.plot(xMinMax,yMinMax,color="Blue")
                            plt.axhline(y=0, color='r', linestyle ="--")
                            plt.axvline(x=MinValue[0], color='purple', linestyle ="--")
                            plt.axvline(x=Elongated[1], color='purple', linestyle ="--")
                            plt.title("MinMax " + MinMaxProtein)
                            plt.xlabel("Center of codon window")
                            plt.ylabel("%MinMax")
                            plt.grid()
                            ImgPath2 = "{}/MinMax_{}Site_{}.png".format(ProteinFolderPath, MinMaxProtein, Sites)
                            plot2.savefig(ImgPath2, dpi=plot2.dpi)
                            plt.clf()

                            continue
                        else:
                            print(">No site found")
                            CoTransFile.write("#Not found\n")
                            continue


            # Close files
            MinMaxFile.close()
            ElongationFile.close()

        # Increment the counter
        Count = Count + 1

CoTransFile.close()
