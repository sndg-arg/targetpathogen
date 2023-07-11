import requests, re, json, os, argparse, subprocess

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

P2RANK_BIN = "/home/rafa_br/Documents/GitHub/p2rank_2.4/prank"
FPOCKET_BIN = "/home/rafa_br/Documents/GitHub/fpocket/bin/fpocket"
MAX_CPU = 4
def GetAFFromUniprot(accession):
    url = f'https://rest.uniprot.org/uniprotkb/{accession}.txt'
    get_request = requests.get(url)
    filename = f"{accession}_uniprot.pdb"
    with open(filename, "wb") as f:
        for chunk in get_request.iter_content(chunk_size=2**20):
            f.write(chunk)
    text = get_request.text
    pattern = r"AlphaFoldDB;\s\S+;"
    af = (re.search(pattern, text).group()).split(';')
    return af[1].strip()

def GetAlphaFoldPrediction(accession):
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}?key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"
    metadata = json.loads(requests.get(url).text)
    pdb_text = requests.get(metadata[0]["pdbUrl"]).text
    with open(f"{accession}_AF.pdb", "w") as f:
        f.write(pdb_text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="protein's PDB accession number", type=str, required=True)
    parser.add_argument('-o', '--output', help="path of the output", type=str, required=False, default=os.getcwd())
    args = parser.parse_args()
    working_dir = args.output
    sys.stderr.write(f"Uniprot protein accession number: {args.input}\n")
    accession = GetAFFromUniprot(args.input)
    GetAlphaFoldPrediction(accession)
    print(f"Uniprot Alpha Fold DB entry: {accession}")
    out_dir = os.path.join(working_dir, "output")
    sts = subprocess.Popen(f"{P2RANK_BIN} predict" +  f" -o {os.path.join(out_dir, accession)} -visualizations 0 -threads {MAX_CPU} -c alphafold -f {accession}_AF.pdb", shell=True).wait()
    # sts = subprocess.Popen(f"{FPOCKET_BIN}" +  f" -f {accession}_AF.pdb", shell=True).wait()
    sts = subprocess.run(f"docker run -v {working_dir}:{working_dir} -w {working_dir} --user {os.getuid()}:{os.getgid()} --rm -it ezequieljsosa/fpocket fpocket -f {os.path.join(working_dir, accession)}_AF.pdb", shell=True)
    fpocket_out = f"{accession}_AF_out"
    path_exists = os.path.exists(os.path.join(working_dir, fpocket_out))
    if path_exists:
        os.system(f"mv {os.path.join(working_dir, fpocket_out)} {os.path.join(out_dir, fpocket_out)}")