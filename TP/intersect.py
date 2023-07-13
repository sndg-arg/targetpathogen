import pandas as pd
import requests, sys, json, os, argparse, subprocess, tqdm, tarfile, shutil
import numpy as np



class Intersector:
    def __init__(self, accession) -> None:
        self.p2rank_file = "/home/rafa_br/Documents/GitHub/targetpathogen/TP/output/P9WPB2/P9WPB2_p2rank/P9WPB2_AF.pdb_predictions.csv"
        self.p2_df = pd.read_csv(self.p2rank_file, 
                                   sep='\s*,\s*',
                                   usecols=["name", "residue_ids"],
                                   skipinitialspace=True, engine='python')
        self.interpro_scan_file = "~/Downloads/iprscan5-R20230713-163003-0060-38200671-p1m.tsv"
        self.interpro_df = pd.read_csv(self.interpro_scan_file,
                                    sep='\t', 
                                    header=None, 
                                    usecols=[3, 4, 6, 7], 
                                    names=["family", "identifier", "range_inf", "range_sup"])
        self.df = None
        
    def CalcIntersections(self):
        colnames = ["Family", "Identifier"]
        for i in range(0, self.p2_df.shape[0]):
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
                adjacent_residues = list(int(x[2:]) for x in adjacent_residues)
                intersections[p2row["name"]] = list()
                for res in adjacent_residues:
                    if res >= range_inf and res <= range_sup:
                        intersections[p2row["name"]].append(str(res))
            values = [cur_fam, cur_identifier]
            for k in intersections.keys():
                values.append((', ').join(intersections[k]))
            self.df.loc[len(self.df.index)] = values
###################################################




if __name__ == "__main__":
    obj = Intersector("a")
    obj.CalcIntersections()
    obj.df.to_csv("teste.csv", sep=',', index=False)