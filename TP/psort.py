from collections import defaultdict
import os
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
import pandas as pd

class Psort:
    DOCKERIMAGENAME = "brinkmanlab"
    DOCKERCONTAIERNAME = "psortb_commandline:1.0.2"

    # If Psort does not exists get the docker image and the psortb wrapper
    def __init__(self, accession):
        self.accession = accession
        if not os.path.exists('psort'):
            print('Psort folder does not exist, creating one...')
            os.makedirs('psort')
            execute(
                f'docker pull {Psort.DOCKERIMAGENAME}/{Psort.DOCKERCONTAIERNAME} && '
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
    parser.add_argument('--tpwebdir', default=os.environ.get("BIOSEQDATADIR", "."))
    parser.add_argument('--timeout', type=int, default=1200,
                        help="Timeout in seconds for the psortb command")
    args = parser.parse_args()

    # Path handling
    #ps = Psort(args.accession)

    tpwebdir = args.tpwebdir
    seqstore = SeqStore(f"{tpwebdir}/data")
    genome_dir = seqstore.db_dir(args.accession)

    faa_decomp = seqstore.faa_decompress(args.accession)

    if not os.path.exists(faa_decomp):
        faa_gz_file = seqstore.faa(args.accession)
        unzip_command = f'gzip -dk {faa_gz_file}' 
        subprocess.run(unzip_command, shell=True)

    faa_file = seqstore.faa_decompress(args.accession)

    tmp = f'{tpwebdir}/tmp'

    # Tmp folder is needed to catch the output file
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)

    # Based on the argument recived choose a command.
    if args.negative:

        command = f'{tpwebdir}/psort/psortb -n -o terse --seq {faa_file} --outdir {tmp}'
    if args.positive:

        command = f'{tpwebdir}/psort/psortb -p -o terse --seq {faa_file} --outdir {tmp}'


    # Run the process and timeout at 20' to avoid endless execution.
    process = subprocess.run(command, shell=True, timeout=1200)
    # try:
    #     timeout = args.timeout
    #     output, errors = process.communicate(timeout=timeout)
    # except subprocess.TimeoutExpired:
    #     print("The command timed out.")
        
    # Use glob to find all.txt files in the directory and move the results to its corresponding directory.
    

    psort_list = glob.glob(f"{tpwebdir}/tmp/*.txt")
    psort_out = psort_list[0]
    destination_file_path = os.path.join(genome_dir, 'psort_res')
    shutil.move(psort_out, destination_file_path)

    psort_res = seqstore.psort(args.accession)
    print(psort_res)
    df = pd.read_csv(psort_res, sep='\t')
    # Modify the SeqID column to only store the first word
    df['SeqID'] = df['SeqID'].apply(lambda x: x.split()[0])
    # Drop the Score column
    df.drop('Score', axis=1, inplace=True)
    # Rename columns
    df.rename(columns={'SeqID': 'gene', 'Localization': 'Celular_localization'}, inplace=True)
    db_dir = seqstore.db_dir(args.accession)
    csv_filename = 'psort.tsv' 
    csv_path = os.path.join(db_dir,csv_filename)
    df.to_csv(csv_path, sep='\t', index=False)
    print(df)



