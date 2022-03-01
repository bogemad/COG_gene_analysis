#!/usr/bin/env python3

import sys, os
from Bio import SeqIO

import sys, csv, os
from collections import defaultdict
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', action='store', dest='csv_file', help='Input csv of target orthogroups and genes')
    parser.add_argument('-p', action='store', dest='proteins', nargs='+', help='Protein fasta files')
    parser.add_argument('-a', action='store', dest='annotations', nargs='+', help='Annotation files from funannotate')
    parser.add_argument('-g', action='store', dest='cogs', nargs='+', help='COG annotations files outputted from extract_COG_annotations.py')
    parser.add_argument('-o', action='store', dest='out', help='Output directory')
    args = parser.parse_args()
    extract(args.csv_file, args.proteins, args.annotations, args.cogs, args.out)

def extract_all(gene_files, prots, anns, COGs, out):
    for gene_file in gene_files:
        extract(gene_file, prots, anns, COGs, out)


def extract(gene_file, prots, anns, COGs, out):
    outg = os.path.join(out, os.path.basename(gene_file)[:-4])
    protein_dir = os.path.join(outg, "protein_fastas")
    annotation_dir = os.path.join(outg, "annotations")
    os.makedirs(protein_dir, exist_ok = True)
    os.makedirs(annotation_dir, exist_ok = True)
    og_d = input_list(gene_file)
    all_proteins = input_fastas(prots)
    anno_d = consolidate_dictionaries([ input_list(a, delim='\t', header=True) for a in anns ])
    with open(anns[0]) as f:
        hl = f.readline().strip('\n').split('\t')
        header = '\t'.join(hl[:5] + hl[6:7] + hl[15:16] + ['COG\n'])
    cog_d = consolidate_dictionaries([ input_list(a) for a in COGs ])
    output_desired_protein_fastas(og_d, all_proteins, protein_dir)
    output_desired_protein_annotations(og_d, anno_d, annotation_dir, header, cog_d)
    output_desired_protein_annotations_in_one_file(annotation_dir, outg, header)


def input_list(csv_file, delim=',', header=False):
    og_d = {}
    with open(csv_file) as infile:
        cr = csv.reader(infile, delimiter=delim)
        if header == True:
            null = next(cr, None)
        for row in cr:
            og_d[row[0]] = row[1:]
    return og_d

def consolidate_dictionaries(d_list):
    for dict in d_list[1:]:
        d_list[0].update(dict)
    return d_list[0]

def input_fastas(proteins):
    dicts = []
    for fasta in proteins:
        d = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
        dicts.append(d)
    return consolidate_dictionaries(dicts)

def output_desired_protein_fastas(og_d, all_proteins, outdir):
    for og, og_prots in og_d.items():
        og_l = []
        for prot in og_prots:
            og_l.append(all_proteins[prot])
        outfile = os.path.join(outdir, "{}.fasta".format(og))
        SeqIO.write(og_l, outfile, "fasta")

def output_desired_protein_annotations(og_d, anno_d, outdir, header, cog_d):
    for og, og_prots in og_d.items():
        og_l = []
        for prot in og_prots:
            anno = prot.split('-')[0]
            og_l.append('\t'.join([anno] + anno_d[anno][:4] + anno_d[anno][5:6] + anno_d[anno][14:15] + [cog_d[anno][0]]))
        outfile = os.path.join(outdir, "{}.annotations.tsv".format(og))
        with open(outfile, 'w') as outh:
            outh.write(header)
            for line in og_l:
                outh.write(line + '\n')

def output_desired_protein_annotations_in_one_file(indir, outdir, header):
    outfile = os.path.join(outdir, "{}.unique.annotations.tsv".format(os.path.basename(outdir)))
    outh = open(outfile, 'w')
    outh.write(header)
    for file in os.listdir(indir):
        with open(os.path.join(indir, file)) as inh:
            outh.write(file.replace('.annotations.tsv','') + '\n')
            null = next(inh)
            outh.write(inh.read())
    outh.close()


if __name__ == '__main__':
    main()