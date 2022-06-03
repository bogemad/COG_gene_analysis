#!/usr/bin/env python

from Bio import SeqIO
import sys, csv, os, argparse
from collections import defaultdict
from pprint import pprint


def output_note_data(qual):
    for k,v in qual.items():
        if k == 'note':
            for item in v[0].split('; '):
                if item.startswith('COG:'):
                    return qual['locus_tag'] + [item[4:]]
    return qual['locus_tag'] + ['']

def extract(gbk, outfile):
    recs = SeqIO.parse(gbk, "genbank")
    out = open(outfile, 'w')
    
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'CDS':
                output = output_note_data(feat.qualifiers)
                out.write("{},{}\n".format(output[0], output[1]))
    out.close()

def import_data(infiles):
    COG_d = defaultdict(list)
    for infile in infiles:
        with open(infile) as inhandle:
            csv_in = csv.reader(inhandle)
            for row in csv_in:
                COG_d[infile].append(row)
    return COG_d

def calculate_totals(cog_list):
    COG_assign_totals_d = defaultdict(int)
    for row in cog_list:
        if row[1] == '':
            COG_assign_totals_d["None"] += 1
        elif len(row) == 2:
            COG_assign_totals_d[row[1]] += 1
        elif len(row) > 2:
            COG_assign_totals_d["Multiple"] += 1
        else:
            print("Issue at row {}".format(" ".join(row)))
    return COG_assign_totals_d

def generate_big_COG_d(infiles, COG_d, COGs):
    big_COG_totals_d = defaultdict(list)
    for COG in COGs:
        big_COG_totals_d[COG].append(COG)
    for infile in infiles:
        COG_assign_totals_d = calculate_totals(COG_d[infile])
        for COG in COGs:
            big_COG_totals_d[COG].append(COG_assign_totals_d[COG])
    return big_COG_totals_d

def basename_l(l):
    return [os.path.basename(v) for v in l]

def output_totals(big_COG_totals_d, infiles, COGs, outfile):
    with open(outfile, 'w') as out:
        csv_out = csv.writer(out)
        csv_out.writerow(["COG"] + basename_l(infiles))
        for COG in COGs:
            csv_out.writerow(big_COG_totals_d[COG])

def generate_summary(infiles, outfile):
    COGs = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','S','T','U','V','W','Y','Z','Multiple','None']
    COG_d = import_data(infiles)
    big_COG_totals_d = generate_big_COG_d(infiles, COG_d, COGs)
    output_totals(big_COG_totals_d, infiles, COGs, outfile)

def import_COGs(COGs):
    COG_d = {}
    with open(COGs) as cogs:
        cogs_csv = csv.reader(cogs)
        for row in cogs_csv:
            COG_d[row[0]] = row[1:]
    return COG_d

def output_gene_COG(COG_d, genes, outfile):
    gene_list = []
    with open(genes) as Genes:
        genes_csv = csv.reader(Genes)
        for row in genes_csv:
            gene_list += [ gene.split('-')[0] for gene in row[1:]]
    with open(outfile, 'w') as out:
        csv_out = csv.writer(out)
        for gene in gene_list:
            csv_out.writerow([gene] + COG_d[gene])

def consolidate_dictionaries(d_list):
    for dict in d_list[1:]:
        d_list[0].update(dict)
    return d_list[0]

def interest_genes(cog_annotation_files, gene_files, out):
    ig_out = os.path.join(out, 'interest_groups_COGs')
    outfiles = []
    os.makedirs(ig_out, exist_ok = True)
    COGds = []
    for COGs in cog_annotation_files:
        COGds.append(import_COGs(COGs))
    COG_d = consolidate_dictionaries(COGds)
    for gene_file in gene_files:
        outfile = os.path.join(ig_out, os.path.splitext(os.path.basename(gene_file))[0] + '_COGs.csv')
        output_gene_COG(COG_d, gene_file, outfile)
        outfiles.append(outfile)
    return outfiles

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', action='store', dest='gbk', nargs='+', help='Genbank files to extract COG annotations from')
    parser.add_argument('-o', action='store', dest='out', help='Output directory')
    parser.add_argument('-e', action='store', dest='genes', nargs='+', help='CSV files showing genes unique to group of interest - from venn_diagram_output')
    args = parser.parse_args()
    cog_ann = []
    for gbk in args.gbk:
        outfile = os.path.join(args.out, os.path.splitext(os.path.basename(gbk))[0] + '.cog_annotations.csv')
        cog_ann.append(outfile)
        extract(gbk, outfile)
    generate_summary(cog_ann, os.path.join(args.out, "COG_summary.csv"))
    interest_genes(cog_ann, args.genes, args.out)



