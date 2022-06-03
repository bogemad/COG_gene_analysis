#!/usr/bin/env python

import sys, csv, os, itertools
from collections import defaultdict
import argparse
import pandas as pd

def write_gene_names_for_interest_groups(row):
    list = [row[0]]
    for item in row[1:]:
        if pd.isnull(item):
            continue
        genes = item.split(', ')
        list += genes
    return list


def generate(assigned, unassigned, baseout):
    iso_d = defaultdict(int)
    out = os.path.join(baseout, "venn_output")
    os.makedirs(out, exist_ok = True)
    
    df = pd.read_table(assigned, header=0)

    empty_cols = [col for col in df.columns if df[col].isnull().all()]
    df.drop(empty_cols, axis=1, inplace=True)

    samples = list(df.columns)[3:]

    interest_groups = []
    for i in range(1, len(samples) + 1):
        interest_groups += itertools.combinations(samples, i)

    file_l = []
    csv_d = {}
    for i, group in enumerate(interest_groups):
        file_l.append(os.path.join(out, "{}.csv".format(".".join(group))))
        csv_d[interest_groups[i]] = csv.writer(open(file_l[i], 'w'))

    df = df.loc[:, ['HOG'] + samples]
    df.dropna(how='all', subset=samples, inplace=True)

    for i, row in df.iterrows():
        iso_list = []
        for col_name in samples:
            if pd.notnull(row[col_name]):
                iso_list.append(col_name)
        outrow = write_gene_names_for_interest_groups(row)
        csv_d[tuple(iso_list)].writerow(outrow)
        iso_d["&".join(iso_list)] += len(write_gene_names_for_interest_groups(row)) - 1
    
    with open(unassigned) as infile:
        udf = pd.read_table(infile, header=0)
    
    udf.dropna(how='all', subset=samples, inplace=True)
    empty_cols = [col for col in udf.columns if udf[col].isnull().all()]
    udf.drop(empty_cols, axis=1, inplace=True)

    for i, row in udf.iterrows():
        iso_list = []
        for col_name in samples:
            if pd.notnull(row[col_name]):
                iso_list.append(col_name)
        outrow = write_gene_names_for_interest_groups(row)
        csv_d[tuple(iso_list)].writerow(outrow)
        iso_d["&".join(iso_list)] += len(write_gene_names_for_interest_groups(row)) - 1

    with open(os.path.join(out, 'theileria.combinations.txt'), 'w') as outfile:
        for key in sorted(iso_d.keys()):
            outfile.write('{}\t{}\n'.format(key, iso_d[key]))
    return file_l

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', action='store', dest='orthogroups', help='Orthogroups.tsv file from orthofinder')
    parser.add_argument('-u', action='store', dest='unassigned', help='Orthogroups_UnassignedGenes.tsv file from orthofinder')
    parser.add_argument('-o', action='store', dest='out', help='Output directory')
    args = parser.parse_args()
    generate(args.orthogroups, args.unassigned, args.out)
