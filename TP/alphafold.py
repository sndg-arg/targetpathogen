import requests, sys, json, os, argparse, subprocess

"""

- FPOCKET_BIN with docker line
- Pasar a Clase
- P2RANK_BIN, FPOCKET_BIN, MAX_CPU default o que se puedan poner en el contructor
- template de la url, tambien variables de la clase
- usar para contro sys.stderr
- pLDDT descargar (scrore por residuo)
- .gz

Procesar una lista

"""
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
        self.result_dir = os.path.join(working_dir, "output")
        self.MAX_CPU = max_cpu
        self.uniprot_text_filename = os.path.join(self.working_dir, f"{accession}_uniprot.txt")
        self.af_pdb_filename = os.path.join(self.working_dir,f"{accession}_AF.pdb")
        sys.stderr.write(f"Uniprot protein accession number: {self.accession}\n" +
                         f"P2RANK binaries: {self.P2RANK_BIN}\n" +
                         f"FPOCKET binaries: {self.FPOCKET_BIN}\n" + 
                         f"Working dir: {self.working_dir}\n" +
                         f"Maximum cores: {self.MAX_CPU}\n"
                         f"Output dir: {self.result_dir}\n"
                         )
        pass

    def GetUniprotFile(self):
        get_request = requests.get(self.uniprot_url)
        with open(self.uniprot_text_filename, "wb") as f:
            for chunk in get_request.iter_content(chunk_size=2**20):
                f.write(chunk)
    
    def GetAlphaFoldPrediction(self):
        metadata = json.loads(requests.get(self.af_url).text)
        pdb_text = requests.get(metadata[0]["pdbUrl"]).text
        with open(self.af_pdb_filename, "w") as f:
            f.write(pdb_text)
    
    def RunP2rankFromFile(self):
        sts = subprocess.Popen(f"{self.P2RANK_BIN} predict" +  f" -o {os.path.join(self.result_dir, f"{self.accession} + ")} -visualizations 0 -threads {self.MAX_CPU} -c alphafold -f {self.accession}_AF.pdb", shell=True, stdout=sys.stderr).wait()
    
    def RunFpocketFromFile(self):
        sts = subprocess.run(self.FPOCKET_BIN + f" -f {self.af_pdb_filename}", shell=True, stdout=sys.stderr)
        fpocket_out = f"{self.accession}_AF_out"
        path_exists = os.path.exists(os.path.join(self.working_dir, fpocket_out))
        if path_exists:
            if not os.path.exists(self.result_dir):
                os.mkdir(self.result_dir)
            os.system(f"mv {os.path.join(self.working_dir, fpocket_out)} {os.path.join(self.result_dir, fpocket_out)}")
        
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="protein's PDB accession number", type=str, required=True)
    parser.add_argument('-o', '--working_dir', help="path of the working_dir", type=str, required=False, default=os.getcwd())
    args = parser.parse_args()
    working_dir = args.working_dir
    obj = AlphaFolder(args.input)
    obj.GetAlphaFoldPrediction()
    obj.RunP2rankFromFile()
    obj.RunFpocketFromFile()