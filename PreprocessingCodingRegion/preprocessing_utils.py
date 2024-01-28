from Bio import SeqIO
import subprocess
import os
import pandas as pd
import io
import multiprocessing as mp

def combineAllMSAGenes(outputFolder, arrayAccessionsLocation, numberOfGenes):
    accessionDict = {}
    accessionDeleteList = list()
    for filename in sorted(os.listdir(outputFolder)):
        if filename.__contains__("aligned"):
            latestLength = 0
            accessions = pd.read_csv(arrayAccessionsLocation, sep="\t")["Accession"].values.tolist()
            for record in SeqIO.parse(os.path.join(outputFolder,filename), "fasta"):
                latestLength = record.__len__()
                print(getAccessionFromId(record.id))
                accessionNumber = getAccessionFromId(record.id)
                if accessionNumber in accessions:
                    accessions.remove(accessionNumber)
                if accessionNumber not in accessionDict:
                    accessionDict[accessionNumber] = ""
                if accessionDict[accessionNumber] == "":
                    accessionDict[accessionNumber] += record.seq
                else:
                    accessionDict[accessionNumber] += "8" + record.seq

            for accession in accessions:
                accessionDeleteList.append(accession)
                if accession not in accessionDict:
                    accessionDict[accession] = ""
                if accessionDict[accession] == "":
                    accessionDict[accession] += "-" * latestLength
                else:
                    accessionDict[accession] += "8" + "-" * latestLength

    array = ""
    for name, sequence in accessionDict.items():
        array += (f">{name}\n{sequence}\n")


    output_file = f"{outputFolder}/sequences_aligned_{numberOfGenes}.fasta"
    with open(output_file, "w") as file:
        file.write(array)
    output_file1 = f"{outputFolder}/accessionDeleteList.txt"
    with open(output_file1, "w") as file:
        file.write(str(accessionDeleteList))
def runMash(outputFolder, filename):
    #todo error checking
    print("started mash on ", filename)
    result = subprocess.run(["mash", "triangle", os.path.join(os.path.dirname(outputFolder), filename), "-E"], capture_output=True, text=True)
    buffer = io.StringIO(result.stdout)
    output_df = pd.read_csv(buffer, sep="	", header=None)
    total_difference = output_df.loc[:, 2].sum()
    #number of sequences in coding region
    number_of_sequences = len([1 for sequence in open(os.path.join(os.path.dirname(outputFolder), filename)) if sequence.startswith(">")])
    return (filename, total_difference / number_of_sequences)

def checkAllDistances(outputFolder, numberOfGenes):
    items = []
    for filename in os.listdir(outputFolder):
        if filename.__contains__(".fasta"):
            items.append((outputFolder, filename))

    dict = {}
    with mp.Pool() as pool:
        for result in pool.starmap(runMash, items):
            dict[result[0]] = result[1]

    sorted_dict = sorted(dict.items(), key=lambda x:x[1], reverse=True)
    print(sorted_dict)
    output_file = f"{outputFolder}/difference.txt"
    with open(output_file, "w") as file:
        file.write(str(sorted_dict))
    topX = sorted_dict[:numberOfGenes]
    return topX

def runMSA(outputFolder, filename):
    #todo error checking
    result = subprocess.run([os.path.join(os.path.dirname(__file__), './make_MSA.sh'), outputFolder, filename.split(".")[0]], capture_output=True, text=True)
    print("done msa on", filename)
    return filename

def runMSAFor(outputFolder, codingRegions):
    items = []
    for item in codingRegions:
        filename = item[0]
        items.append((outputFolder, filename))

    #todo error checking
    dict = {}
    with mp.Pool() as pool:
        for result in pool.starmap(runMSA, items):
            print("MSA ", result, " finished")

def runMSAForEveryGene(outputFolder):
    for filename in os.listdir(outputFolder):
        result = subprocess.run([os.path.join(os.path.dirname(__file__), './make_MSA.sh'), outputFolder, filename.split(".")[0]], capture_output=True, text=True)
        print("done msa on", filename)

def getAccessionFromId(recordId):
    if(recordId.__contains__("join")):
        return recordId.replace("join(","").replace(")", "").split(",")[0].split(":")[0]
    else:
        return recordId.split("|")[0].split(":")[0]
