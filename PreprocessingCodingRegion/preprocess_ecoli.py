from Bio import SeqIO
import os
import argparse
import time
import pandas as pd
from collections import Counter

from preprocessing_utils import getAccessionFromId, checkAllDistances, runMSAFor, combineAllMSAGenes

def get_gene_list(coding_region_location):
    path = "./data/"
    first_gene_list = list()
    for index, directory in enumerate(os.listdir(path)):
        if os.path.isdir(os.path.join(path, directory)):
            temp_gene_list = list()
            for record in SeqIO.parse(path +"/" + directory +"/" + coding_region_location, "fasta"):
                if len(record.description.split("gene=")) > 1:
                    gene_name = record.description.split("gene=")[1].split("]")[0]
                    gene_name = gene_name.replace("/", "--")
                    if gene_name.__contains__("truncated"):
                        gene_name = gene_name.split("truncated ")[1]
                    gene_name = gene_name.strip()
                    temp_gene_list.append(gene_name)
            c = Counter(temp_gene_list)
            temp_gene_list_counted = [x for x, v in c.items() if v == 1]
            if len(first_gene_list) == 0:
                print("list length", len(temp_gene_list_counted))
                first_gene_list = temp_gene_list_counted
            else:
                print("list length", len(temp_gene_list_counted))
                first_gene_list = set(first_gene_list).intersection(set(temp_gene_list_counted))

    print(first_gene_list)
    print(len(first_gene_list))
    return first_gene_list

def createForEveryRecord(codingRegionLocation, outputFolder, geneList):
    path = "./data/"
    gene_dict = {}
    for directory in os.listdir(path):
        if os.path.isdir(os.path.join(path, directory)):
            for record in SeqIO.parse(path+"/"+directory+"/"+codingRegionLocation, "fasta"):
                if len(record.description.split("gene=")) > 1:
                    gene_name = record.description.split("gene=")[1].split("]")[0]
                    if gene_name in geneList:
                        gene_name = gene_name.replace("/", "--")
                        if gene_name.__contains__("truncated"):
                            gene_name = gene_name.split("truncated ")[1]
                        if gene_name.__contains__("ORF1a "):
                            continue
                        gene_name = gene_name.strip()
                        if gene_name.__contains__("ORF"):
                            gene_name = gene_name.split(" ")[0]
                        if gene_name not in gene_dict:
                            gene_dict[gene_name] = ""
                        gene_dict[gene_name] += f">{directory}|\n{record.seq}\n"

    for gene_name, sequences in gene_dict.items():
        output_file = f"{outputFolder}/{gene_name}.fasta"
        with open(output_file, "w") as file:
            file.write(sequences)
        print(f"gene {gene_name}")



def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description='Run AmpliDiff to find discriminatory amplicons and corresponding primers in input sequences.')

    parser.add_argument('codingRegions', type=str, help='Path to aligned sequences fasta file')
    parser.add_argument('metadata', type=str, help='Path to the metadata for the sequences')
    parser.add_argument('-o', '--output', type=str, default='.', help='Path to the folder where output will stored, default is current folder')
    parser.add_argument('-cr', '--coding_regions', type=int, default=5, help='Number of coding regions that are selected, default is 5')
    #Parse arguments
    args = parser.parse_args()

    #make sure output directory is empty
    if len(os.listdir(args.output)) > 0:
        raise ValueError("Output directory must be empty")

    gene_list = get_gene_list(args.codingRegions)
    createForEveryRecord(args.codingRegions, args.output, gene_list)
    top_differentiating_regions = checkAllDistances(args.output, args.coding_regions)
    runMSAFor(args.output, top_differentiating_regions)
    combineAllMSAGenes(args.output, args.metadata, args.coding_regions)

    end_time = time.time()
    duration = end_time - start_time
    print(duration)


if __name__ == '__main__':
    main()