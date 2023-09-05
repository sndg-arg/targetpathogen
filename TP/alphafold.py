import requests, sys, json, os, argparse, subprocess
import tqdm, tarfile, shutil, glob, re, io
import Bio.PDB as biopdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
class AlphaFolder:
    """Class to manipulate the Pockets of an alpha fold file from its PDB accession number
    """
    def __init__(self, accession, working_dir=os.getcwd(), p2rank_bin=None, fpocket_bin=None, max_cpu=1) -> None:
        if p2rank_bin is not None:
            self.P2RANK_BIN = p2rank_bin
        else:
            self.P2RANK_BIN = "/home/rafa_br/Documents/GitHub/p2rank_2.4/prank"
        if fpocket_bin is not None:
            self.FPOCKET_BIN = fpocket_bin
        else:
            self.FPOCKET_BIN = f"docker run -v {working_dir}:{working_dir} -w {working_dir} --user {os.getuid()}:{os.getgid()} --rm -it ezequieljsosa/fpocket fpocket"
        self.working_dir = working_dir
        self.accession = accession
        self.uniprot_url = f'https://rest.uniprot.org/uniprotkb/{accession}.txt'
        self.af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}?key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"
        self.result_dir = os.path.join(
            working_dir, os.path.join("output", self.accession))
        if not os.path.exists(self.result_dir):
            os.mkdir(self.result_dir)
        self.out_dir = os.path.join(self.working_dir, "output")
        self.MAX_CPU = max_cpu
        self.uniprot_text_filename = os.path.join(
            self.result_dir, f"{accession}_uniprot.txt")
        self.af_pdb_filename = os.path.join(
            self.result_dir, f"{accession}_AF.pdb")
        self.af_plddt_filename = os.path.join(
            self.result_dir, f"{accession}_pldd.json")
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
        if os.path.isfile(self.uniprot_text_filename):
            return None
        get_request = requests.get(self.uniprot_url)
        with open(self.uniprot_text_filename, "wb") as f:
            for chunk in get_request.iter_content(chunk_size=2**20):
                f.write(chunk)

    def GetAlphaFoldPrediction(self):
        """Get the PDB file containing the AlphaFold prediction from its DB. The file is stored at output/{accesion}
        """
        if os.path.isfile(self.af_pdb_filename):
            return None
        metadata = json.loads(requests.get(self.af_url).text)
        pdb_text = requests.get(metadata[0]["pdbUrl"]).text
        with open(self.af_pdb_filename, "w") as f:
            f.write(pdb_text)

    def RunP2rankFromFile(self):
        """Run P2RANK app from the PDB file obtained from AlphaFold DB and saves the results on output/{accesion} folder
        """
        dir_ = os.path.join(self.result_dir, self.accession + '_p2rank')
        if os.path.isdir(dir_):
            return None
        sts = subprocess.Popen(f"{self.P2RANK_BIN} predict" +
                               f" -o {dir_} -visualizations 0 -threads {self.MAX_CPU} -c alphafold -f {self.af_pdb_filename}", shell=True, stdout=sys.stderr).wait()

    def RunFpocketFromFile(self):
        """Run FPOCKET app from the PDB file obtained from AlphaFold DB and saves the results on output/{accesion} folder
        """
        fpocket_out = f"{self.accession}_AF_out"
        dir_ = os.path.join(self.result_dir, fpocket_out)
        if os.path.isdir(dir_):
            return None
        sts = subprocess.run(
            self.FPOCKET_BIN + f" -f {self.af_pdb_filename}", shell=True, stdout=sys.stderr)
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
        tar_filename = os.path.join(self.out_dir, f"{self.accession}.tar.gz")
        with tarfile.open(tar_filename, mode='w:gz') as archive:
            archive.add(self.result_dir, self.accession, recursive=True)
        if remove_files:
            shutil.rmtree(self.result_dir)

    def CompareResults(self, threshold = 0.2):
        fpocket_dir = os.path.join(self.result_dir, f"{self.accession}_AF_out")
        p2rank_dir = os.path.join(self.result_dir, self.accession + '_p2rank')
        p2rank_file = os.path.join(p2rank_dir, f"{self.accession}_AF.pdb_predictions.csv")
        pockets_folder = os.path.join(fpocket_dir, "pockets")
        # ---- fpocket ----
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
                pockets_name_val[pocket_name] = set(df_p.iloc[:, 5])
        pockets_name_val = dict(sorted(pockets_name_val.items()))
        for k in pockets_name_val.keys():
            pockets_name_val[k] = sorted(pockets_name_val[k])

        # ---- p2rank ----
        p2_df = pd.read_csv(p2rank_file,
                            sep='\s*,\s*',
                            usecols=["name", "residue_ids"],
                            skipinitialspace=True, engine='python')
        p2_residues = dict()
        for j, p2row in p2_df.iterrows():
            adjacent_residues = set(p2row["residue_ids"].split(' '))
            adjacent_residues = set(int(x[2:]) for x in adjacent_residues)
            p2_residues[p2row["name"]] = sorted(adjacent_residues)


        '''
        print("FPOCKET pockets:")
        for k in pockets_name_val.keys():
            print(f"-----------------{k}-----------------")
            print(f"Residues: {pockets_name_val[k]}")
        print("\nP2rank pockets:")
        for k in p2_residues.keys():
            print(f"-----------------{k}-----------------")
            print(f"Residues: {p2_residues[k]}")
        '''
        mat = np.zeros((len(p2_residues), len(pockets_name_val)))

        for i, k in enumerate(p2_residues.keys()):
            for j,l in enumerate(pockets_name_val.keys()):
                intersec = list(set(p2_residues[k]).intersection(set(pockets_name_val[l])))
                mat[i][j] = len(intersec)/len(p2_residues[k])
                if mat[i][j] >= threshold:
                    print(f"{mat[i][j]*100}% of P2RANK's {k} found on FPOCKET's {l} pocket.")

        ax = sns.heatmap(mat, linewidth=0.5)
        ax.set_title("Percentage of P2RANK in FPOCKET pockets")
        ax.set_xlabel("FPOCKET pockets")
        plt.xticks(rotation=45)
        ax.set_xticklabels(pockets_name_val.keys())
        ax.set_ylabel("P2RANK pockets")
        ax.set_yticklabels(p2_residues.keys())
        plt.tight_layout()
        plt.savefig(os.path.join(self.result_dir, "comparison_p2rank_fpocket.pdf"), dpi=300)

        t = list(pockets_name_val.keys())
        ticks = dict()
        for j,l in enumerate(t):
            for k in range(j+1, len(t)):
                ticks[f"{j} + {k}"] = set(pockets_name_val[t[j]]).union(pockets_name_val[t[k]])
        mat2 = np.zeros(((len(p2_residues), len(ticks.keys()))))
        for i, k in enumerate(p2_residues.keys()):
            for j,l in enumerate(ticks.keys()):
                intersec = list(set(p2_residues[k]).intersection(set(ticks[l])))
                mat2[i][j] = len(intersec)/len(p2_residues[k])
                if mat2[i][j] >= threshold:
                    print(f"{mat2[i][j]*100}% of P2RANK's {k} found on FPOCKET's {l} pocket union.")
        print(mat2.shape)
        print(ticks.keys())
        plt.clf()
        ax = sns.heatmap(mat2, linewidth=0.5)
        ax.set_title("Percentage of P2RANK in FPOCKET pockets")
        ax.set_xlabel("FPOCKET pockets")
        ax.set_xticks(np.arange(0, len(ticks.keys())))
        ax.set_xticklabels(ticks.keys())
        plt.xticks(rotation=90)
        ax.set_ylabel("P2RANK pockets")
        ax.set_yticklabels(p2_residues.keys())
        plt.tight_layout()
        plt.savefig(os.path.join(self.result_dir, "comparison_p2rank_fpocket_union.pdf"), dpi=300)

if __name__ == "__main__":
    accessions = list()
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help="List of protein's PDB accession numbers separated with new lines",
                        type=str,
                        nargs='*',
                        default=sys.stdin)
    parser.add_argument('-o', '--working_dir', help="path of the working_dir",
                        type=str, required=False, default=os.getcwd())
    parser.add_argument('-nc', '--no_compress',
                        help="flag to not compress the results", action="store_true", default=False)
    parser.add_argument(
        '-T', '--threads', help="Number of threads to be used by the supported programs", type=int, default=1)
    args = parser.parse_args()
    for l in args.input:
        accessions.append(l.strip().upper())
    working_dir = args.working_dir
    for ac in tqdm.tqdm(accessions):
        obj = AlphaFolder(ac, max_cpu=args.threads)
        obj.GetUniprotFile()
        obj.GetAlphaFoldPrediction()
        obj.RunP2rankFromFile()
        obj.RunFpocketFromFile()
        obj.GetPlddtFromFile()
        obj.CompareResults()
        
        if not args.no_compress:
            obj.CompressResults()
