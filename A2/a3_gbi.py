import Bio
from Bio import SeqIO



for record in SeqIO.parse(r"material/unknown_protein.fasta", "fasta"):
    print(record.seq)


# Needleman-Wunsch

