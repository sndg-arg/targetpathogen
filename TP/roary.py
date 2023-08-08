import glob
import os
import subprocess
import sys
import argparse
from pathlib import Path


class Roarier:
    def __init__(self, kingdom, genus, working_dir=os.getcwd()) -> None:
        self.working_dir = working_dir
        self.kingdom = kingdom
        self.genus = genus
        self.out_dir = os.path.join(self.working_dir, "output")
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        self.gffs_folder = os.path.join(self.out_dir, "gff_files")
        if not os.path.exists(self.gffs_folder):
            os.mkdir(self.gffs_folder)
        self.PROKKA_bin = f"docker run --user {os.getuid()}:{os.getgid()} --rm -it -v {self.working_dir}:/data staphb/prokka:latest prokka"
        self.roary_bin = f"docker run --user {os.getuid()}:{os.getgid()} --rm -it -v {self.out_dir}:/data sangerpathogens/roary roary"

    def run_prokka(self) -> None:
        genomes = glob.glob(os.path.join(self.working_dir, "*.fna"))
        for g in genomes:
            locustag = Path(Path(g).stem).stem
            PROKKA_params = f"--kingdom {self.kingdom} --outdir /data/output/prokka/PROKKA_{locustag} --genus {self.genus} --locustag {locustag} {os.path.basename(g)} --force"
            subprocess.run(f"{self.PROKKA_bin} {PROKKA_params}",
                           shell=True, stdout=sys.stderr)
            PROKKA_out_folder = os.path.join(os.path.join(
                self.out_dir, "prokka"), f"PROKKA_{locustag}")
            # get only the first gff file (usually there is only one)
            gff_file = glob.glob(os.path.join(PROKKA_out_folder, "*.gff"))[0]
            os.system(f"cp {gff_file} {self.gffs_folder}")
            gff_file = os.path.join(
                self.gffs_folder, os.path.basename(gff_file))
            os.rename(gff_file, os.path.join(
                self.gffs_folder, 'PROKKA_' + locustag + '.gff'))

    def run_roary(self) -> None:
        gff_list = glob.glob(os.path.join(self.gffs_folder, "*.gff"))
        for i in range(0, len(gff_list)):
            gff_list[i] = os.path.basename(gff_list[i])
            gff_list[i] = os.path.join("/data/gff_files", gff_list[i])
        roary_out = "roary"
        roary_params = f"-f /data/{roary_out} -e -n -v {(' ').join(gff_list)}"
        subprocess.run(f"{self.roary_bin} {roary_params}",
                        shell=True, stdout=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Roary pipeline")
    parser.add_argument(
        "working_dir", help="path of a directory containing the .gff files", type=str)
    parser.add_argument("kingdom", help="kingdom of the sequences", type=str)
    parser.add_argument("genus", help="genus of the sequences", type=str)
    args = parser.parse_args()
    obj = Roarier(kingdom=args.kingdom, genus=args.genus,
                  working_dir=args.working_dir)
    obj.run_prokka()
    obj.run_roary()
