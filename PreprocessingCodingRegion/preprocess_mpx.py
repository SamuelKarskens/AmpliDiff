from Bio import SeqIO
import os
import argparse
import time
import pandas as pd
from collections import Counter

from preprocessing_utils import getAccessionFromId, checkAllDistances, runMSAFor, combineAllMSAGenes

def getGeneList(codingRegionLocation):
    dict = {}
    for record in SeqIO.parse(codingRegionLocation, "fasta"):
        accession = getAccessionFromId(record.id)
        gene_name = record.description.split("[")[0].split(record.id)[1]
        gene_name = gene_name.strip()
        if accession not in dict:
            dict[accession] = list()
        dict[accession].append(gene_name)

    total_list = list()
    for item in dict:
        c = Counter(dict[item])
        list_only_singular_genes = [x for x, v in c.items() if v == 1]
        total_list.append(set(list_only_singular_genes))
    intersected = set.intersection(*total_list)
    print(intersected)
    print(len(intersected))
    return intersected
def createForEveryRecord(codingRegionLocation, outputFolder, geneList):
    gene_dict = {}
    for record in SeqIO.parse(codingRegionLocation, "fasta"):
        gene_name = record.description.split("[")[0].split(record.id)[1]
        gene_name = gene_name.strip()
        if gene_name in geneList:
            gene_name = gene_name.replace("/", "--")
            if gene_name not in gene_dict:
                gene_dict[gene_name] = ""
            gene_dict[gene_name] += f">{record.id}|\n{record.seq}\n"

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

    geneList = getGeneList(args.codingRegions)
    createForEveryRecord(args.codingRegions, args.output, geneList)
    top_differentiating_regions = checkAllDistances(args.output, args.coding_regions)
    runMSAFor(args.output, top_differentiating_regions)
    combineAllMSAGenes(args.output, args.metadata, args.coding_regions)

    end_time = time.time()
    duration = end_time - start_time
    print(duration)


if __name__ == '__main__':
    main()