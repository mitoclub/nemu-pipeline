import re

from Bio import SeqIO
import click


@click.command("formatter", help="reformat midori fasta headers")
@click.option("-i", "--inp", required=True, type=click.Path(True), help="path to input fasta file")
@click.option("-o", "--out", required=True, type=click.Path(writable=True), help="path to output fasta file")
def main(inp, out):
    seqs = {}
    for rec in SeqIO.parse(inp, "fasta"):
        header = rec.description
        raw_acc, taxa = header.split("###")
        acc, place = re.match("(\w+\.\d)\.(.+)", raw_acc).groups()
        taxa = taxa.removeprefix("root_1;Eukaryota_2759;")
        species = " ".join(taxa.split(";")[-1].split("_")[:-1])
        rec.id = acc
        rec.description = f"{species} {place} ###{taxa}"
        if acc in seqs:
            if len(seqs[acc]) < len(rec):
                seqs[acc] = rec
        else:
            seqs[acc] = rec

    SeqIO.write(list(seqs.values()), out, "fasta")


if __name__ == "__main__":
    main()
