import requests
import sys
import json
import os
import argparse
import subprocess
import tqdm
import tarfile
import shutil
import glob
import re
import io
import Bio.PDB as biopdb
import pandas as pd
import numpy as np


class AlphaFolder:
    """Class to manipulate the Pockets of an alpha fold file from its PDB accession number
    """

    def __init__(self, accession, locus_tag , working_dir=os.getcwd(), results_dir=os.getcwd() + "/output", p2rank_bin=None,
                 fpocket_bin=None,in_parsl=None, max_cpu=1) -> None:
        if p2rank_bin is not None:
            self.P2RANK_BIN = p2rank_bin
        else:
            self.P2RANK_BIN = "prank"
        if fpocket_bin is not None:
            self.FPOCKET_BIN = fpocket_bin
        else:
            CWD = os.environ.get('CWD', '/')
            if in_parsl is None:
                self.FPOCKET_BIN = f"docker run -v {CWD}:{working_dir} -w {working_dir} --user {os.getuid()}:{os.getgid()} --rm ezequieljsosa/fpocket fpocket"
            else:
                parent_dir = os.path.dirname(working_dir)
                self.FPOCKET_BIN = f"docker run -v {CWD}:{parent_dir} -w {working_dir} --user {os.getuid()}:{os.getgid()} --rm ezequieljsosa/fpocket fpocket"
        self.working_dir = working_dir
        self.accession = accession
        self.locus_tag = locus_tag
        self.uniprot_url = f'https://rest.uniprot.org/uniprotkb/{accession}.txt'
        self.af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}?key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"
        self.result_dir = os.path.join(
            results_dir, self.locus_tag)
        if not os.path.exists(self.result_dir):
            os.makedirs(self.result_dir )
        self.out_dir = self.result_dir
        self.MAX_CPU = max_cpu
        self.uniprot_text_filename = os.path.join(
            self.result_dir, f"{locus_tag}_uniprot.txt")
        self.af_pdb_filename = os.path.join(
            self.result_dir, f"{locus_tag}_af.pdb")
        self.af_plddt_filename = os.path.join(
            self.result_dir, f"{locus_tag}_pldd.json")
        
        sys.stderr.write(f"Uniprot protein accession number: {self.accession}\n" +
                         f"P2RANK binaries: {self.P2RANK_BIN}\n" +
                         f"FPOCKET binaries: {self.FPOCKET_BIN}\n" +
                         f"Working dir: {self.working_dir}\n" +
                         f"Maximum cores: {self.MAX_CPU}\n"
                         f"Output dir: {self.result_dir}\n"
                         )
        pass

    def GetUniprotFile(self):
        """Get the text file containing the features of a protein from Uniprot. The file is stored at output/{accesion}
        """
        try:
            get_request = requests.get(self.uniprot_url)
        except requests.exception.RequestException:
            sys.stderr.write(f"Protein {self.accession} not found on Uniprot DB!\n")
            return None
        if "Error" in get_request.text:
            sys.stderr.write(f"Protein {self.accession} not found on Uniprot DB!\n")
            return None
        with open(self.uniprot_text_filename, "wb") as f:
            for chunk in get_request.iter_content(chunk_size=2 ** 20):
                f.write(chunk)

    def GetAlphaFoldPrediction(self):
        """Get the PDB file containing the AlphaFold prediction from its DB. The file is stored at output/{accesion}
        """
        if os.path.isfile(self.af_pdb_filename):
            return None
        try:
            metadata = json.loads(requests.get(self.af_url).text)
            pdb_text = requests.get(metadata[0]["pdbUrl"]).text
        except (requests.exceptions.RequestException, KeyError):
            sys.stderr.write(f"AlphaFold prediction for protein {self.accession} not found in DB!\n")
            return None
        with open(self.af_pdb_filename, "w") as f:
            f.write(pdb_text)

    def RunP2rankFromFile(self):
        """Run P2RANK app from the PDB file obtained from AlphaFold DB and saves the results on output/{accesion} folder
        """
        if not os.path.exists(self.P2RANK_BIN):
            print('entro')
            # Clone P2RANK repository if it doesn't exist
            clone_path = os.path.join('../opt', 'p2rank')
            os.makedirs(clone_path, exist_ok=True)
            subprocess.run(['git', 'clone', 'https://github.com/rdk/p2rank.git', clone_path], check=True)
            # Change directory to the cloned repository
            os.chdir(clone_path)
            # Make the script executable
            subprocess.run(['chmod', '+x', 'make.sh'], check=True)
            # Run the make.sh script
            subprocess.run(['./make.sh'], check=True)
            # Update self.P2RANK_BIN to point to the cloned binary
            self.P2RANK_BIN = os.path.join(clone_path, 'distro/prank')  # Assuming the binary is located here
        dir_ = os.path.join(self.result_dir, self.locus_tag + '_p2rank')
        if os.path.isdir(dir_):
            return None
        try:
            sts = subprocess.Popen(f"{self.P2RANK_BIN} predict" +
                               f" -o {dir_} -threads {self.MAX_CPU} -c alphafold -f {self.af_pdb_filename}", shell=True, stdout=sys.stderr).wait()
        except subprocess.CalledProcessError:
            sys.stderr.write(f"Failed to run P2RANK for the {self.accession} AlphaFold.\n")             

    def RunFpocketFromFile(self):
        """Run FPOCKET app from the PDB file obtained from AlphaFold DB and saves the results on output/{accesion} folder
        """
        fpocket_out = f"{self.accession}_af_out"
        dir_ = os.path.join(self.result_dir, fpocket_out)
        if os.path.isdir(dir_):
            return None
        sts_arg = self.FPOCKET_BIN + f" -f {self.af_pdb_filename}"
        
        try:

            sts = subprocess.run(
                self.FPOCKET_BIN + f" -f {self.af_pdb_filename}", shell=True, stdout=sys.stderr)
            
        
        except subprocess.CalledProcessError:
            sys.stderr.write(f"Failed to run Fpocket for the {self.accession} AlphaFold.\n")
        path_exists = os.path.exists(
            os.path.join(self.working_dir, fpocket_out))
        if path_exists:
            os.system(
                f"mv {os.path.join(self.working_dir, fpocket_out)} {dir_}")

    def GetPlddtFromFile(self, save_to_file=True):
        """Get the PLDDT values from all the residues in the PDB file from AlphaFold's DB

        Args:
            save_to_file (bool, optional): Defaults to True.

        Returns:
            Dictionary containing all the pair residue, value
        """
        if os.path.isfile(self.af_plddt_filename):
            return None
        plddt = dict()
        parser = biopdb.PDBParser()
        structure = parser.get_structure(self.accession, self.af_pdb_filename)
        for r in structure.get_residues():
            atoms = list(r.get_atoms())
            plddt[r.get_resname()] = atoms[0].get_bfactor()
        if save_to_file:
            with open(self.af_plddt_filename, 'w') as jfile:
                json.dump(plddt, jfile)
        return plddt

    def CompressResults(self, remove_files=True):
        """Compress all the results present in the protein's accession number folder

        Args:
            remove_files (bool, optional): Defaults to True.
        """
        tar_filename = os.path.join(self.out_dir, f"{self.locus_tag}.tar.gz")
        with tarfile.open(tar_filename, mode='w:gz') as archive:
            archive.add(self.result_dir, self.locus_tag, recursive=True)
        if remove_files:
            shutil.rmtree(self.result_dir)

    def CompareResults(self, threshold=0.7):
        fpocket_dir = os.path.join(self.result_dir, f"{self.locus_tag}_af_out")
        p2rank_dir = os.path.join(self.result_dir, self.locus_tag + '_p2rank')
        p2rank_file = os.path.join(
            p2rank_dir, f"{self.locus_tag}_af.pdb_predictions.csv")
        pockets_folder = os.path.join(fpocket_dir, "pockets")
        # ---- fpocket ----
        #  files contain only the atoms contacted by alpha spheres in the given pocket.
        pockets = glob.glob(os.path.join(pockets_folder, "*.pdb"))
        pockets_name_val = dict()
        pockets_name_score = dict()
        for pocket in pockets:
            with open(pocket, 'r') as pocket_f:
                data = pocket_f.read()
                values = re.findall(r"\bATOM.*", data)
                df_p = pd.read_csv(io.StringIO(
                    '\n'.join(values)), delim_whitespace=True, header=None)
                pocket_name = os.path.basename(pocket).split('_')[0]
                pockets_name_val[pocket_name] = set(df_p.iloc[:, 5])

                score = re.findall("HEADER\s0.+", data)
                pockets_name_score[pocket_name] = float(score[0].split(":")[1])
        pockets_name_val = dict(sorted(pockets_name_val.items()))
        for k in pockets_name_val.keys():
            pockets_name_val[k] = sorted(pockets_name_val[k])

        # ---- p2rank ----
        p2_df = pd.read_csv(p2rank_file,
                            sep='\s*,\s*',
                            usecols=["name", "residue_ids", "score"],
                            skipinitialspace=True, engine='python')
        p2_residues = dict()
        for j, p2row in p2_df.iterrows():
            adjacent_residues = set(p2row["residue_ids"].split(' '))
            adjacent_residues = set(int(x[2:]) for x in adjacent_residues)
            p2_residues[p2row["name"]] = sorted(adjacent_residues)


        #---------------------------
        print(f"Threshold: {threshold}")
        print("fpocket")
        for pocket in pockets_name_val.keys():
            print(f"Size of {pocket}: {len(pockets_name_val[pocket])} -> Score {pockets_name_score[pocket]}")
        print("p2rank")
        for pocket in p2_residues.keys():
            print(f"Size of {pocket}: {len(p2_residues[pocket])} -> Score: {p2_df.set_index('name').at[pocket, 'score']}")


        #---------------------------
        sorted_p2 = dict(sorted(p2_residues.items(),
                         key=lambda x: len(x[1]), reverse=True))
        selected = dict()
        print("FPOCKET in P2")
        for p2_pocket in sorted_p2.keys():
            sorted_fp = list()
            for x in pockets_name_val.keys():
                inter = set(pockets_name_val[x]).intersection(
                    set(sorted_p2[p2_pocket]))
                sorted_fp.append((x, len(inter)))
            sorted_fp = sorted(sorted_fp, key=lambda x: x[1], reverse=True)
            for x in sorted_fp:
                if (x[1]/len(pockets_name_val[x[0]])) >= threshold:
                    if selected.get(p2_pocket) is None:
                        selected[p2_pocket] = list()
                        selected[p2_pocket].append(x[0])
                    else:
                        selected[p2_pocket].append(x[0])
                    del pockets_name_val[x[0]]
        print(selected)
        # ----------------------------------------------------------
        for s in selected.keys():
            del p2_residues[s]

        print("P2 in FPOCKET")
        sorted_pockets = dict(
            sorted(pockets_name_val.items(), key=lambda x: len(x[1]), reverse=True))
        selected2 = dict()
        for fp_pocket in sorted_pockets.keys():
            sorted_p2 = list()
            for x in p2_residues.keys():
                inter = set(p2_residues[x]).intersection(
                    set(sorted_pockets[fp_pocket]))
                sorted_p2.append((x, len(inter)))
            sorted_p2 = sorted(sorted_p2, key=lambda x: x[1], reverse=True)
            for x in sorted_p2:
                if (x[1]/len(p2_residues[x[0]])) >= threshold:
                    if selected2.get(fp_pocket) is None:
                        selected2[fp_pocket] = list()
                        selected2[fp_pocket].append(x[0])
                    else:
                        selected2[fp_pocket].append(x[0])
                    del p2_residues[x[0]]
        print(selected2)
        print(f"FPOCKET: {pockets_name_val}")
        print(f"P2: {p2_residues}")


if __name__ == "__main__":
    accessions = list()
    parser = argparse.ArgumentParser()
    parser.add_argument('accessions', help="List of protein's PDB accession numbers separated with new lines",
                        type=str,
                        nargs='*',
                        default=sys.stdin)
    parser.add_argument('-o', '--results_dir', help="path of the results",
                        type=str, required=False, default=os.getcwd())
    parser.add_argument('-w', '--working_dir', help="path of the working_dir",
                        type=str, required=False, default=os.getcwd())
    parser.add_argument('-nc', '--no_compress',
                        help="flag to not compress the results", action="store_true", default=False)
    parser.add_argument('-nf', '--no_fpocket',
                        help="flag to not run fpocket", action="store_true", default=False)
    parser.add_argument('-np', '--no_p2rank',
                        help="flag to not run p2rank", action="store_true", default=False)
    parser.add_argument('-na', '--no_alphafold',
                        help="flag to not retrieve the alphafold, PLDDT and uniprot model.\
                            Important: to run the other functions, an alphafold model needs to be in the working dir with the name ACESSION_af.pdb", action="store_true", default=False)
    parser.add_argument('-pr', '--p2rank_bin', required=False,
                        help="p2rank binary path", default=None)
    parser.add_argument('-c', '--compare', required=False,
                        help="flag to compare results", action="store_true", default=None)
    parser.add_argument('-ltag', '--locus_tag', required=False,
                        help="Locus tag to be used in the dic structure",type=str, default=None)
    parser.add_argument('-T', '--threads', help="Number of threads to be used by the supported programs", type=int, default=1)
    parser.add_argument('-parsl', '--in_parsl', required=False,
                        help="Flag to indicate that this script is being runned within a parsl pipeline 'usually in the bottom /parsl folder'",  action="store_true", default=None)
    args = parser.parse_args()

    for l in args.accessions:
        accessions.append(l.strip().upper())

    for ac in tqdm.tqdm(accessions):
        obj = AlphaFolder(ac, locus_tag= args.locus_tag,in_parsl=args.in_parsl, p2rank_bin=args.p2rank_bin, working_dir=args.working_dir, results_dir=args.results_dir, max_cpu=args.threads)
        if not args.no_alphafold:
            try:
                obj.GetUniprotFile()
            except:
                sys.stderr.write(f"Error while trying to get the uniprot file for the protein {ac}")
            try:
                obj.GetAlphaFoldPrediction()
            except:
                sys.stderr.write(f"Error while trying to get the alphafold model file for the protein {ac}!\
                                 Skipping execution!")
                pass # the alphafold model is essential to the execution of the next steps
            
            try:
                obj.GetPlddtFromFile()
            except:
                sys.stderr.write(f"Error while trying to get the PLDDT file for the protein {ac}!\
                                 Skipping execution!")
        if not args.no_p2rank:
            try:
                obj.RunP2rankFromFile()
            except:
                sys.stderr.write(f"Error while trying to run P2RANK for protein {ac}!\n")
        if not args.no_fpocket:
            try:
                obj.RunFpocketFromFile()
            except:
                sys.stderr.write(f"Error while trying to run FPOCKET for protein {ac}!\n")
        if args.compare:
            try:
                obj.CompareResults()
            except:
                sys.write.stderr(f"Error while trying to compare the results \
                                 from P2RANK and FPOCKET for\
                                 protein {ac}\n")

        if not args.no_compress:
            try:
                obj.CompressResults()
            except:
                sys.write.stderr(f"Error while trying to compress \
                                 the results for protein {ac}\n")
