"""

https://github.com/tseemann/prokka
"""
from Bio.SeqFeature import SeqFeature, FeatureLocation

import fileinput
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.SeqIO as bpio
import sys


def lt(f):
    return f.qualifiers.get("locus_tag", [""])[0]


class GenebankUtils:
    """

    """

    def __init__(self):
        pass

    def validategb(self, gbk_handle, handle_out):
        # for contig in contigs:
        #     features = []
        #     for f in contig:
        #         if f.type in ["CDS","mRNA"]
        genome_protein_count = 0
        for contig in gbk_handle:
            contig_protein_count = 0
            for f in contig.features:
                if f.type == "mat_peptide":
                    contig_protein_count += 1
                    seq = str(f.extract(contig).seq.translate())
                    if seq.count("*") > 1:
                        handle_out.write(
                            f"mat_peptide {contig.id}:{f.location.start}-{f.location.end} has more than one stop codon\n")

                elif f.type == "CDS":
                    contig_protein_count += 1
                    if "locus_tag" not in f.qualifiers:
                        handle_out.write(
                            f"CDS {contig.id}:{f.location.start}-{f.location.end}  has no locus_tag\n")
                    if "translation" not in f.qualifiers and not f.qualifiers.get("pseudo", []):
                        handle_out.write(
                            f"CDS {contig.id}:{f.location.start}-{f.location.end}  has no translation property\n")
                    elif "translation" in f.qualifiers:
                        seq = Seq(f.qualifiers["translation"][0])
                        if seq.count("*") > 1:
                            handle_out.write(
                                f"CDS {contig.id}:{f.location.start}-{f.location.end}  has more than one stop codon\n")

            genome_protein_count += contig_protein_count
            handle_out.write(f"{contig.name} protein_count: {contig_protein_count}\n")

        if not genome_protein_count:
            handle_out.write('Genome has no proteins!!!\n')

    def fixgb(self, h_gb, hw_gb, new_lt=None):
        contigs = list(bpio.parse(h_gb, "gb"))
        remove = []
        genome_protein_count = 0
        ltnum = 1
        ltnum_matpep = 1
        for contig in contigs:
            for f in contig.features:
                plasmid = False
                if f.type == "source":
                    if "plasmid" in f.qualifiers:
                        plasmid = True
                if f.type == "mat_peptide":
                    genome_protein_count += 1
                    curr_cds.qualifiers["polyprotein"] = ["true"]
                    if "note" in f.qualifiers and (len(f.qualifiers["note"][0]) < 6):
                        f.qualifiers["gene"] = f.qualifiers["note"][0]
                    elif "product" in f.qualifiers and (len(f.qualifiers["product"][0].split()[0]) <= 6):
                        f.qualifiers["gene"] = f.qualifiers["product"][0].split()[0]
                    if not lt(f):
                        f.qualifiers["locus_tag"] = [lt(curr_cds) + "_" + str(ltnum_matpep).zfill(3)]
                        ltnum_matpep += 1
                    if "translation" not in f.qualifiers:
                        seq = str(f.extract(contig).seq.translate())
                        f.qualifiers["translation"] = [seq]
                        if seq.count("*") > 1:
                            remove.append(f)

                elif f.type == "CDS":
                    if plasmid:
                        f.qualifiers["plasmid"] = ["True"]
                    genome_protein_count += 1
                    curr_cds = f
                    if "locus_tag" not in f.qualifiers:
                        if not new_lt:
                            sys.stderr.write(
                                f"CDS {f.location.start}:{f.location.end} has more than one stop codon\n exiting...\n")
                            remove.append(f)
                            # sys.exit(1)
                        else:
                            f.qualifiers["locus_tag"] = [new_lt + "_" + str(ltnum).zfill(5)]
                            ltnum += 1
                    if "translation" not in f.qualifiers:
                        seq = str(f.extract(contig).seq.translate())
                        if seq.count("*") > 1:
                            f.qualifiers["pseudo"] = ["true"]
                        else:
                            f.qualifiers["translation"] = [seq]

        if genome_protein_count:
            sys.stderr.write(f"Genome has f{genome_protein_count} proteins\n")
            bpio.write(contigs, hw_gb, "gb")
        else:
            sys.stderr.write("Genome has no proteins!!!")
        sys.stderr.write("Finished!")

    def lt(self, feature):
        return feature.qualifiers.get("locus_tag", [""])[0]

    def region_from_lts(self, lts, contig, region_name, rtype="REGION", strand_lts=None):

        region = []
        strand = None

        for f in contig.features:
            if self.lt(f) in lts:
                region.append(f)
                if strand_lts and (self.lt(f) in strand_lts):
                    strand = f.location.strand

        locations = [x.location.start for x in region] + [x.location.end for x in region]
        if locations:
            return SeqFeature(type=rtype, location=FeatureLocation(
                start=min(locations), end=max(locations), strand=strand), qualifiers={"name": [region_name]})
        return None

    def proteins(self, sequences, otype="prot"):

        for contig in sequences:
            for protein in self.proteins_from_sequence(contig, otype):
                yield protein

    def contig4target(self, contig):
        ant_cds = None
        polyprots = []
        for feature in contig.features:
            if feature.type in ["CDS", "RNA", "mat_peptide"]:
                location = f'{contig.id}:{feature.location.start}-{feature.location.end}'

                locus_tag = feature.qualifiers.get("locus_tag", [location + "_" + feature.type])[0]
                gene = feature.qualifiers["gene"][0] if "gene" in feature.qualifiers else ""
                if (feature.type == "CDS"):
                    ant_cds = feature

                if feature.type == "mat_peptide":
                    ant_cds.qualifiers["polyprotein"] = ["true"]
                    gene = gene + "_" + feature.qualifiers["product"][0]
                    locus_tag = locus_tag + "_" + feature.qualifiers["product"][0].replace(" ", "_")
                    feature.qualifiers["gene"] = [gene]
                    feature.qualifiers["locus_tag"] = [locus_tag]
                    polyprots += ant_cds.qualifiers["locus_tag"]
                    feature.qualifiers["from_polyprotein"] = ant_cds.qualifiers["locus_tag"]
                    feature.type = "CDS"

                if ((feature.type == "CDS") and ("pseudo" not in feature.qualifiers) and (
                        "translation" not in feature.qualifiers)) or (feature.type == "mat_peptide"):

                    seq = str(feature.extract(contig.seq).translate())
                    feature.qualifiers["translation"] = [seq]
            contig.features = [f for f in contig.features if f.qualifiers.get("locus_tag",[""])[0] not in polyprots]
        return contig

    def proteins_from_sequence(self, contig, otype="prot"):
        org = (" [" + contig.annotations["organism"].replace("[", "_").replace("]", "_") + "]"
               ) if "organism" in contig.annotations else ""
        for feature in contig.features:
            if feature.type in ["CDS", "RNA", "mat_peptide"]:
                seq = None
                location = f'{contig.id}:{feature.location.start}-{feature.location.end}'
                description = location + " " + feature.qualifiers["product"][0] if "product" else (
                    feature.qualifiers["note"][0] if "note" in feature.qualifiers else "")

                locus_tag = feature.qualifiers.get("locus_tag", [location + "_" + feature.type])[0]
                gene = feature.qualifiers["gene"][0] if "gene" in feature.qualifiers else ""

                if feature.type == "mat_peptide":
                    gene = gene + "_" + feature.qualifiers["product"][0]
                    locus_tag = locus_tag + "_" + feature.qualifiers["product"][0].replace(" ", "_")
                if otype == "prot":
                    if feature.type == "CDS" and "pseudo" not in feature.qualifiers:
                        seq = Seq(feature.qualifiers["translation"][0])
                    elif feature.type == "mat_peptide":
                        seq = feature.extract(contig.seq).translate()
                else:
                    seq = feature.extract(contig.seq)
                if seq:
                    record = SeqRecord(id=locus_tag, name=gene, description=description + org, seq=seq)
                    yield record
                else:
                    assert "pseudo" in feature.qualifiers
    
    @staticmethod
    def download_assembly(assembly_accession, dst_dir, dtype="genomic.gbff.gz", force=False):
        # assembly_name, last_assembly_accession = NCBI.assembly_name_from_acc(assembly_accession)
        assembly_accession_no_ver = assembly_accession if assembly_accession[-2] != "." else assembly_accession[:-2]

        # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/435/
        url = "/".join([ftp_url
                           , assembly_accession_no_ver[0:3], assembly_accession_no_ver[4:7]
                           , assembly_accession_no_ver[7:10], assembly_accession_no_ver[10:13]
                        # , last_assembly_accession + "_" + assembly_name.replace(" ", "_").replace("#", "_")
                        # , last_assembly_accession + "_" + assembly_name.replace(" ", "_").replace("#",
                        #                                                                           "_") + "_" + dtype
                        ]) + "/"
        r = requests.get(url)
        download_url = ""
        acc =""
        if r.status_code == 200:
            accessions = [x.split("</")[0].replace("/","") for x in r.text.split(">") if x.startswith(assembly_accession_no_ver)]
            # GCF_000158435.(1)_ASM15843v1/
            accessions = sorted(accessions, key=lambda x: int(x.split("_")[1][-1]))
            if accessions:
                acc = accessions[-1]
                download_url = f'{url}{acc}/{acc}_{dtype}'

        if not download_url:
            err = f"{assembly_accession} not found at {url}"
            _log.error(err)
            raise FileNotFoundError(err)

        assert acc

        out_file = f'{dst_dir}/{"_".join(acc.split("_")[:2]) }.{dtype}'
        if force or not os.path.exists(out_file):
            download_file(download_url, out_file, ovewrite=force)
        else:
            _log.debug(f'{out_file} exists')
        # execute("gunzip -c  " + out_file + " > " +  out_file[:-3])
        return out_file


if __name__ == '__main__':
    import argparse
    import os
    from tqdm import tqdm

    parser = argparse.ArgumentParser(description='Utils over genebank file')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)
    cmd = subparsers.add_parser('fix', help='fix genebank file to be imported')
    cmd.add_argument('input_bgk')
    cmd.add_argument('output_gb', nargs='?', default=sys.stdout)
    cmd.add_argument('--new_lt', default=None, help="new locus tag. Use ONLY if no locus tag is present")

    cmd = subparsers.add_parser('genes', help='extracts the list of gene products from ')
    cmd.add_argument('input_bgk')
    cmd.add_argument('output_fna', nargs='?', default=sys.stdout)
    cmd.add_argument('-otype', help="output type", choices=["prot", "nucl"], default="prot")

    cmd = subparsers.add_parser('4target', help='extracts the list of gene products from ')
    cmd.add_argument('input_bgk')
    cmd.add_argument('output_gbk', nargs='?', default=sys.stdout)

    cmd = subparsers.add_parser('validate', help='validates genebank file')
    cmd.add_argument('input_bgk')

    cmd = subparsers.add_parser('download', help='downloads an NCBI accession')
    cmd.add_argument( 'accession')
    cmd.add_argument('-o', '--output', help="output_directory", default="./")
    cmd.add_argument('--force', action="store_true")

    args = parser.parse_args()
    utils = GenebankUtils()

    if args.command == "4target":
        if args.input_bgk != "-":
            h = open(args.input_bgk)
        else:
            h = sys.stdin
        for contig in bpio.parse(h,"gb"):
            bpio.write(utils.contig4target(contig), sys.stdout, "gb")

    if args.command == "validate":
        gbk_h = bpio.parse(fileinput.input(args.input_bgk), "gb")
        utils.validategb(gbk_h, sys.stderr)

    if args.command == "genes":
        if isinstance(args.input_bgk, str):
            if not os.path.exists(args.input_bgk):
                sys.stderr.write(f'{args.input_bgk} not found')
                sys.exit(1)
        if args.input_bgk != "-":
            gbk_h = bpio.parse(args.input_bgk, "gb")
        else:
            gbk_h = bpio.parse(fileinput.input(args.input_bgk), "gb")

        if isinstance(args.output_fna, str):
            h_faa = open(args.output_fna, "w")
        else:
            h_faa = args.output_fna
        try:
            for seq in utils.proteins(gbk_h, otype=args.otype):
                bpio.write(seq, h_faa, "fasta")
        finally:
            try:
                h_faa.close()
            except:
                pass

    if args.command == "fix":
        try:
            if args.input_bgk != "-":
                h = open(args.input_bgk)
            else:
                h = sys.stdin
            utils.fixgb(h, args.output_gb, args.new_lt)
        finally:
            h.close()

    if args.command == "download":
        downloaded_file = GenebankUtil.sdownload_assembly(args.accession, dst_dir=args.output, force=args.force)
        if not os.path.exists(downloaded_file) or (os.path.getsize(downloaded_file) < 10):
            sys.stderr.write(f"error downloading {args.accession}")
        else:
            sys.stdout.write(downloaded_file)
            sys.stdout.write("\n")
