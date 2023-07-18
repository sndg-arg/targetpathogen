import pandas as pd
import requests, sys, json, os, argparse, subprocess, tqdm, tarfile, shutil
import numpy as np


"""
    pathwaytools

    1 levantar (docker)
    Sever -- DB
    2  correr el cliente
    pathologi genebak
   3 API obtenemos archivo que necesitamos

   -----------------------

   tax SA -- genomas (Refseq ; complete/contig) -- armar pangenoma (roary)

"""
class Intersector:
    """Class to calculate points of intersection between a protein and its predicted pockets
    """
    def __init__(self, accession, p2rank_file, interpro_scan_file, working_dir = os.getcwd()) -> None:
        self.working_dir = working_dir
        self.accession = accession
        self.result_dir = os.path.join(working_dir, os.path.join("output", self.accession))
        if not os.path.exists(self.result_dir):
                os.mkdir(self.result_dir)
        self.out_dir = os.path.join(self.working_dir, "output")
        self.p2rank_file = p2rank_file
        self.interpro_scan_file = interpro_scan_file
        pass        
        
    def CalcIntersectionsP2(self, save_to_file = False):
        """Calculates the points of intersection between a protein and its predicted pockets from P2Rank

        Args:
            save_to_file (bool, optional): Defaults to False.

        Returns:
            pandas.DataFrame: a pandas dataframe containing the points of intersection
        """
        p2_df = pd.read_csv(self.p2rank_file, 
                                   sep='\s*,\s*',
                                   usecols=["name", "residue_ids"],
                                   skipinitialspace=True, engine='python')
        interpro_df = pd.read_csv(self.interpro_scan_file,
                                    sep='\t', 
                                    header=None, 
                                    usecols=[3, 4, 6, 7], 
                                    names=["family", "identifier", "range_inf", "range_sup"])
        colnames = ["Family", "Identifier"]
        for i in range(0, p2_df.shape[0]):
            colnames.append(f"Points of intersection with Pocket {i+1}")
        df = pd.DataFrame(columns=colnames)
        for ind, row in interpro_df.iterrows():
            cur_fam = row["family"]
            cur_identifier = row["identifier"]
            range_inf = row["range_inf"]
            range_sup = row["range_sup"]
            intersections = dict()
            for j, p2row in p2_df.iterrows():
                adjacent_residues = list(p2row["residue_ids"].split(' '))
                adjacent_residues = list(int(x[2:]) for x in adjacent_residues)
                intersections[p2row["name"]] = list()
                for res in adjacent_residues:
                    if res >= range_inf and res <= range_sup:
                        intersections[p2row["name"]].append(str(res))
            values = [cur_fam, cur_identifier]
            for k in intersections.keys():
                values.append((' ').join(intersections[k]))
            df.loc[len(df.index)] = values

        if save_to_file:
            df.to_csv(os.path.join(self.result_dir, self.accession) + "_p2ranks_intersections.csv", sep=',', index=False)
        return df
    
    def CalcIntersectionsFPocket():
        """To do!
        """
        pass



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #p2rank_file = "/home/rafa_br/Documents/GitHub/targetpathogen/TP/output/P9WPB2/P9WPB2_p2rank/P9WPB2_AF.pdb_predictions.csv"
    #interpro_scan_file = "~/Downloads/iprscan5-R20230713-163003-0060-38200671-p1m.tsv"
    parser.add_argument("accession", help="protein's PDB accession number", type=str)
    parser.add_argument("interprotscan", help="File containing the protein's InterprotScan metadata", type=str)
    parser.add_argument("p2rankfile", help="File containing the protein's predicted pockets using P2Rank", type=str)
    parser.add_argument('-o', '--working_dir', help="path of the working_dir", type=str, required=False, default=os.getcwd())
    args = parser.parse_args()

    obj = Intersector(args.accession, args.p2rankfile, args.interprotscan, args.working_dir)
    df = obj.CalcIntersectionsP2(save_to_file=True)
    print(df)
