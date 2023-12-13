import time

import pandas as pd
import requests
import sys
import json
import os
import argparse
import subprocess as sp
import tqdm
import tarfile
import shutil
import Bio.PDB as biopdb
import tempfile



class FoldSeek:

    DEFAULT_DATABASES = """-F 'database[]=pdb100'"""

    def __init__(self):
        self.cmd_search = "curl -X POST -F q=@{filepath} -F 'mode=3diaa' {outputfmt} {databases} https://search.foldseek.com/api/ticket"

    def submitJob(self, filepath, extra=None):
        if not extra:
            extra = {"outputfmt": "", "databases": FoldSeek.DEFAULT_DATABASES}
        extra["filepath"] = filepath
        cmd = self.cmd_search.format(**extra)
        sys.stderr.write(cmd)
        res = json.loads(sp.check_output(cmd, shell=True).decode("utf-8"))
        # {"id":"PvH__CrMjPXgueMorGjjuElX2KI2IHrAcXkC2w","status":"PENDING"}
        return FoldSeekJob(res["id"])


class FoldSeekJob:
    def __init__(self, ticket, outfile=None, sleeptime=3):
        self.ticket = ticket
        self.sleeptime = sleeptime
        self.outfile = outfile
        if not self.outfile:
            self.outfile = tempfile.NamedTemporaryFile(delete=False)
            self.outfile.close()
            self.outfile = self.outfile.name
        self.df = None
        self.cmd_status = "curl -X GET https://search.foldseek.com/api/ticket/{ticket}"
        self.cmd_download = "curl -X GET https://search.foldseek.com/api/result/download/{ticket} > {outfile}"

    def wait_complete(self):
        res = json.loads(sp.check_output(self.cmd_status.format(ticket=self.ticket), shell=True).decode("utf-8"))

        # {"id":"UFMNJNBWBpGD26rdNpCfr6RorMkbzF-0GZFy1w","status":"COMPLETE"}
        if res["status"] == "ERROR":
            return Exception("error procesing FoldSeek")
        elif res["status"] == "PENDING":
            time.sleep(self.sleeptime)
            self.checkStatus()
        elif res["status"] == "COMPLETE":
            return
        else:
            return Exception("unknown FoldSeek server status:" + res["status"])

    def download_result(self):
        sp.call(self.cmd_download.format(ticket=self.ticket, outfile=self.outfile), shell=True)

    def process_result(self):
        cols = "query pdb_desc fident alnlen mismatch gapopen	qstart qend tstart tend x10	x11	x12	x13	x14	seq1 seq2 nums seq3 tax taxname".split()

        time.sleep(1) # needs time to close the file before reading
        self.df = pd.read_csv(self.outfile, sep="\t", names=cols, index_col=False, compression='gzip')
        self.df = self.df.dropna(subset=["pdb_desc"])
        self.df["pdb"] = [x.split("_")[0] for x in self.df.pdb_desc]
        self.df["chain"] = [x.split("_")[1].split()[0] for x in self.df.pdb_desc]

    def save(self, newfile):
        shutil.copy(self.outfile, newfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_file', type=str)
    parser.add_argument('--wait', type=int, default=3)

    parser.add_argument('--databases', default="",
                        help='''list separateby coma. Default=pdb100. Options:
                            afdb50, afdb-swissprot, afdb-proteome, cath50, mgnify_esm30, pdb100, gmgcl_id''')

    args = parser.parse_args()

    assert os.path.exists(args.pdb_file), f'"{args.pdb_file}" does not exist'

    fs = FoldSeek()
    extra = {}

    #    extra["outputfmt"] = "F='--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage'"
    if args.databases:
        extra["databases"] = FoldSeek.DEFAULT_DATABASES

    job = fs.submitJob(args.pdb_file, extra=extra)


    job.wait_complete()
    job.download_result()
    job.process_result()
    sys.stderr.write(f"{job.outfile}\n")

    sys.stdout.write(job.df.to_csv(sep="\t",index=False))

