import pandas as pd
import requests, sys, json, os, argparse, subprocess, tqdm, tarfile, shutil
import numpy as np



class Intersector:
    def __init__(self, accession) -> None:
        self.p2rank_file = ""
        self.p2_df = pd.read_csv(self.p2rank_file, 
                                   sep='\s*,\s*',
                                   usecols=["name", "residue_ids"],
                                   skipinitialspace=True, engine='python')
        self.interpro_scan_file = ""
        self.interpro_df = pd.read_csv(interpro_scan_file,
                                    sep='\t', 
                                    header=None, 
                                    usecols=[3, 4, 6, 7], 
                                    names=["family", "identifier", "range_inf", "range_sup"])
        self.df = None
        
    def CalcIntersections(self):
        colnames = ["Family", "Identifier"]
        for i in range(0, df.shape[1]):
            colnames.append(f"Points of intersection with Pocket {i+1}")
        self.df = pd.DataFrame(columns=colnames)
        for ind, row in self.interpro_df.iterrows():
            cur_fam = row["family"]
            cur_identifier = row["identifier"]
            range_inf = row["range_inf"]
            range_sup = row["range_sup"]
            intersections = dict()
            for j, p2row in self.p2_df.iterrows():
                adjacent_residues = list(p2row["residue_ids"].split(' '))
                adjacent_residues = np.asarray(int(x[2:]) for x in adjacent_residues)
                intersections[p2row["name"]] = np.where(adjacent_residues >= range_inf and adjacent_residues <= range_sup, adjacent_residues)
                for res in adjacent_residues:
                    #if res >= range_inf and res <= range_sup:
                    #    intersections[p2row["name"]].append(str(res))
            values = [cur_name, cur_identifier]
            for k in intersections.keys():
                values.append((', ').join(intersections[k]))
            df.loc[len(df.index)] = values
###################################################
p2rank_file = "/home/rafa_br/Documents/GitHub/targetpathogen/TP/output/P9WPB2/P9WPB2_p2rank/P9WPB2_AF.pdb_predictions.csv"

p2rank_predictions = pd.read_csv(p2rank_file, sep='\s*,\s*', usecols=["name", "residue_ids"], skipinitialspace=True, engine='python')

print("P2Rank predictions")
print(p2rank_predictions)


interpro_scan_file = "~/Downloads/iprscan5-R20230713-163003-0060-38200671-p1m.tsv"
interpro_scan = pd.read_csv(interpro_scan_file, sep='\t', header=None, usecols=[3, 4, 6, 7], names=["name", "identifier", "range_inf", "range_sup"])
print("InterproScan")
print(interpro_scan)
lnames = ["Family", "Identifier"]
lnames.extend("Points of intersection on " + p2rank_predictions["name"])
df = pd.DataFrame(columns=lnames)



print("-----------------------------------")
print(df)
df.to_csv("teste.csv", sep=',', index=False)
