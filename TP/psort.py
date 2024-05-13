from collections import defaultdict

from bioseq.io.SeqStore import SeqStore
import pythoncyc
import pickle
import networkx as nx
from threading import Thread
from TP import execute
import gzip
import shutil
from Bio import SeqIO
from pathlib import Path
import argparse
import os
from time import sleep
import subprocess
import glob

class Psort:
    DOCKERIMAGENAME = "brinkmanlab"
    DOCKERCONTAIERNAME = "psortb_commandline:1.0.2"

    #
    def __init__(self, accession):

        self.accession = accession
        if not os.path.exists('psort'):
            print('Psort folder does not exist, creating one...')
            os.makedirs('psort')
            execute(
                f'sudo docker pull {Psort.DOCKERIMAGENAME}/{Psort.DOCKERCONTAIERNAME} && '
                f'wget -O psort/psortb https://raw.githubusercontent.com/brinkmanlab/psortb_commandline_docker/master/psortb && '
                f'chmod +x psort/psortb'
            )
        

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Pathologic')
    parser.add_argument('accession', 
                        help="Accession code of the genome")
    parser.add_argument('-p', '--positive', action="store_true",
                        help="Indicates that the bacteria is a Gram positive", default=False)
    parser.add_argument('-n', '--negative', action="store_true",
                        help="Indicates that the bacteria is a Gram negative", default=False)
    args = parser.parse_args()

    ps = Psort(args.accession)
    seqstore = SeqStore("./data")
    genome_dir = seqstore.db_dir(args.accession)
    faa_gz_file = seqstore.faa(args.accession)

    unzip_command = f'gzip -dk {faa_gz_file}' 
    subprocess.run(unzip_command, shell=True)
    faa_file = seqstore.faa_decompress(args.accession)

    if os.path.exists('./tmp'):
        shutil.rmtree('tmp')
    os.makedirs('tmp')

    if args.negative:
        command = f'psort/psortb -n -o terse --seq {faa_file} --outdir tmp'
    if args.positive:
        command = f'psort/psortb -p -o terse --seq {faa_file} --outdir tmp'
    print("Starting command execution...")
    process = subprocess.Popen(command, shell=True, text=True)
    print("Command execution completed.")
    
    try:
        output, errors = process.communicate(timeout=1200)  # Adjust the timeout as needed
        print(output)
    except subprocess.TimeoutExpired:
        print("The command timed out.")
        
    # Use glob to find all.txt files in the directory
    psort_list = glob.glob(f'./tmp/*.txt')
    psort_out = psort_list[0]
    destination_file_path = os.path.join(genome_dir, 'psort_res')
    print(psort_out)
    print(destination_file_path)

    shutil.move(psort_out, destination_file_path)               


