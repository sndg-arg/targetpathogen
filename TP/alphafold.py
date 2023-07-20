import requests, sys, json, os, argparse, subprocess, tqdm, tarfile, shutil
import Bio.PDB as biopdb

class AlphaFolder:

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
        self.result_dir = os.path.join(working_dir, os.path.join("output", self.accession))
        if not os.path.exists(self.result_dir):
                os.mkdir(self.result_dir)
        self.out_dir = os.path.join(self.working_dir, "output")
        self.MAX_CPU = max_cpu
        self.uniprot_text_filename = os.path.join(self.result_dir, f"{accession}_uniprot.txt")
        self.af_pdb_filename = os.path.join(self.result_dir,f"{accession}_AF.pdb")
        self.af_plddt_filename = os.path.join(self.result_dir,f"{accession}_pldd.json")
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
        sts = subprocess.Popen(f"{self.P2RANK_BIN} predict" +  f" -o {dir_} -visualizations 0 -threads {self.MAX_CPU} -c alphafold -f {self.af_pdb_filename}", shell=True, stdout=sys.stderr).wait()

    
    def RunFpocketFromFile(self):
        """Run FPOCKET app from the PDB file obtained from AlphaFold DB and saves the results on output/{accesion} folder
        """
        fpocket_out = f"{self.accession}_AF_out"
        dir_ = os.path.join(self.result_dir, fpocket_out)
        if os.path.isdir(dir_):
            return None
        sts = subprocess.run(self.FPOCKET_BIN + f" -f {self.af_pdb_filename}", shell=True, stdout=sys.stderr)
        path_exists = os.path.exists(os.path.join(self.working_dir, fpocket_out))
        if path_exists:
            os.system(f"mv {os.path.join(self.working_dir, fpocket_out)} {dir_}")
        
    def GetPlddtFromFile(self, save_to_file = True):
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
    
    def CompressResults(self, remove_files = True):
        """Compress all the results present in the protein's accession number folder

        Args:
            remove_files (bool, optional): Defaults to True.
        """
        tar_filename = os.path.join(self.out_dir, f"{self.accession}.tar.gz")
        with tarfile.open(tar_filename, mode='w:gz') as archive:
            archive.add(self.result_dir, self.accession, recursive=True)
        if remove_files:
            shutil.rmtree(self.result_dir)


if __name__ == "__main__":
    accessions = list()
    parser = argparse.ArgumentParser()
    if not sys.stdin.isatty():
        # there is something in the stdin
        parser.add_argument('stdinput', nargs='?', type=argparse.FileType('r'),
                             default=(None if sys.stdin.isatty() else sys.stdin),
                             help="List of protein's PDB accession numbers separated with new lines")
    else:
        parser.add_argument('-i', '--input', help="protein's PDB accession number", type=str, required=True)
    parser.add_argument('-o', '--working_dir', help="path of the working_dir", type=str, required=False, default=os.getcwd())
    parser.add_argument('-nc', '--no_compress', help="flag to not compress the results",action="store_true", default=False)
    args = parser.parse_args()
    if not sys.stdin.isatty():
        for l in args.stdinput.readlines():
            accessions.append(l.strip().upper())
    else:
        accessions.append(args.input)
    working_dir = args.working_dir
    for ac in tqdm.tqdm(accessions):
        obj = AlphaFolder(ac)
        obj.GetUniprotFile()
        obj.GetAlphaFoldPrediction()
        obj.RunP2rankFromFile()
        obj.RunFpocketFromFile()
        obj.GetPlddtFromFile()
        if not args.no_compress:
            obj.CompressResults()