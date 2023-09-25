import pandas as pd
import sys
import os
import argparse
import glob
import re
import io
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

    def __init__(self, accession, interpro_scan_file, working_dir=os.getcwd()) -> None:
        self.working_dir = working_dir
        self.accession = accession
        self.result_dir = os.path.join(
            working_dir, os.path.join("output", self.accession))
        if not os.path.exists(self.result_dir):
            os.mkdir(self.result_dir)
        self.out_dir = os.path.join(self.working_dir, "output")
        self.interpro_scan_file = interpro_scan_file
        pass

    def CalcIntersectionsP2(self, p2rank_file, save_to_file=False):
        """Calculates the points of intersection between a protein and its predicted pockets from P2Rank

        Args:
            save_to_file (bool, optional): Defaults to False.

        Returns:
            pandas.DataFrame: a pandas dataframe containing the points of intersection
        """
        p2_df = pd.read_csv(p2rank_file,
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
            df.to_csv(os.path.join(self.result_dir, self.accession) +
                      "_p2ranks_intersections.csv", sep=',', index=False)
        return df

    def CalcIntersectionsFPocket(self, pockets_folder, save_to_file=False):
        """Calculates the points of intersection between a protein and its predicted pockets from FPOCKET

        Args:
            save_to_file (bool, optional): Defaults to False.

        Returns:
            pandas.DataFrame: a pandas dataframe containing the points of intersection
        """
        interpro_df = pd.read_csv(self.interpro_scan_file,
                                  sep='\t',
                                  header=None,
                                  usecols=[3, 4, 6, 7],
                                  names=["family", "identifier", "range_inf", "range_sup"])

        # Read all the pdb files of each pocket
        fpocket_dir = os.path.join(self.result_dir, f"{self.accession}_AF_OUT")
        pockets_folder = os.path.join(fpocket_dir, "pockets")
        #  files contain only the atoms contacted by alpha spheres in the given pocket.
        pockets = glob.glob(os.path.join(pockets_folder, "*.pdb"))
        pockets_name_val = dict()
        for pocket in pockets:
            with open(pocket, 'r') as pocket_f:
                data = pocket_f.read()
                values = re.findall(r"\bATOM.*", data)
                df_p = pd.read_csv(io.StringIO(
                    '\n'.join(values)), delim_whitespace=True, header=None)
                pocket_name = os.path.basename(pocket).split('_')[0]
                pockets_name_val[pocket_name] = list(df_p.iloc[:, 5])
        pockets_name_val = dict(sorted(pockets_name_val.items()))
        colnames = ["Family", "Identifier"]
        for i, k in enumerate(pockets_name_val.keys()):
            colnames.append(f"Points of intersection with Pocket {i}")
        df = pd.DataFrame(columns=colnames)
        # Create a row for each family
        for ind, row in interpro_df.iterrows():
            cur_fam = row["family"]
            cur_identifier = row["identifier"]
            range_inf = row["range_inf"]
            range_sup = row["range_sup"]
            intersections = dict()
            for k in pockets_name_val.keys():
                intersections[k] = list()
                for res in pockets_name_val[k]:
                    if res >= range_inf and res <= range_sup:
                        if str(res) not in intersections[k]:
                            intersections[k].append(str(res))
            values = [cur_fam, cur_identifier]
            for k in intersections.keys():
                values.append((' ').join(intersections[k]))
            df.loc[len(df.index)] = values
        if save_to_file:
            df.to_csv(os.path.join(self.result_dir, self.accession) +
                      "_fpocket_intersections.csv", sep=',', index=False)
        return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # p2rank_file = "/home/rafa_br/Documents/GitHub/targetpathogen/TP/output/P9WPB2/P9WPB2_p2rank/P9WPB2_AF.pdb_predictions.csv"
    # interpro_scan_file = "~/Downloads/iprscan5-R20230713-163003-0060-38200671-p1m.tsv"
    parser.add_argument(
        "-a", "--accession", help="protein's PDB accession number", type=str, required=True)
    parser.add_argument("-i", "--interprotscan",
                        help="File containing the protein's InterprotScan metadata", type=str, required=True)
    parser.add_argument(
        "-p", "--p2rankfile", help="File containing the protein's predicted pockets using P2Rank", type=str, required=False)
    parser.add_argument(
        "-f", "--fpocket", help="Path containing the protein's pockets predicted pockets using FPOCKET", type=str, required=False)
    parser.add_argument("-o", '--working_dir', help="path of the working_dir",
                        type=str, required=False, default=os.getcwd())
    parser.add_argument("-s", "--save", help="Save the result to a CSV file",
                        action="store_true", default=False)
    args = parser.parse_args()
    obj = Intersector(args.accession, args.interprotscan, args.working_dir)
    if args.p2rankfile is not None:
        df = obj.CalcIntersectionsP2(args.p2rankfile, args.save)
    if args.fpocket is not None:
        obj.CalcIntersectionsFPocket(args.fpocket, args.save)
