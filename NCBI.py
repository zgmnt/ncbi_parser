from Bio import SeqIO
from Bio import Entrez
import os

Entrez.email =   # Always tell NCBI who you are


def get_ids(type, gene):
    IDs = []
    if type not in ["nucleotide", "protein"]:
        return IDs

    gene = "[" + gene + "]"
    handle = Entrez.esearch(db=type, term=gene)
    record = Entrez.read(handle)
    retmax_number = record["Count"]

    handle = Entrez.esearch(db=type, term=gene, retmax=retmax_number)
    record = Entrez.read(handle)
    IDs = record["IdList"]
    handle.close()

    return IDs


def do_everything(type, ids, gene):
    lenids = len(ids)
    # check gene folder
    path = "../../../Desktop/Tetracycline/" + gene
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)


    if type not in ["nucleotide", "protein"]:
        return

    # fasta in right folder
    if type == "nucleotide":
        seq_name = "/nucleotide/"
    else:
        seq_name = "/protein/"

    for id in ids:
        idsindex = ids.index(id) + 1
        print("record", idsindex, " out of ", lenids, "   ", (idsindex/lenids)*100, "%")
        handle = Entrez.efetch(db=type, id=id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        name = record.features[0].qualifiers["organism"][0]

        # FASTA
        handle = Entrez.efetch(db=type, id=id, rettype="fasta", retmode="text")
        fasta = handle.read()

        # check for organism folders
        path = "../../../Desktop/Tetracycline/" + gene + "/" + name
        if os.path.exists(path):
            pass
        else:
            os.mkdir(path)

        if type == "protein":
            # check for protein folders
            path = "../../../Desktop/Tetracycline/" + gene + "/" + name + "/protein"
            if os.path.exists(path):
                pass
            else:
                os.mkdir(path)
        else:
            # check for nuc folders
            path = "../../../Desktop/Tetracycline/" + gene + "/" + name + "/nucleotide"
            if os.path.exists(path):
                pass
            else:
                os.mkdir(path)

        # check fasta in folder

        path = "../../../Desktop/Tetracycline/" + gene + "/" + name + seq_name + id + ".fasta"

        if os.path.exists(path):
            continue
        else:
            with open(path, "w") as file:
                file.write(fasta)

        handle.close()


def search(gene):
    type_seraches = ["nucleotide", "protein"]
    for type_serach in type_seraches:
        print("type search ", type_serach, type_seraches.index(type_serach), " out of ", len(type_seraches))
        ids = get_ids(type_serach, gene)
        do_everything(type_serach, ids, gene)


gene = input("type gene...")
search(gene)
input("done")
