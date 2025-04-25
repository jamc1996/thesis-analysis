# Calculation of the number of tryptic and non-specific peptides in human uniprot
# database. 2 missed cleavages are allowed for trypsin.

from Bio import SeqIO
import re
from inspire.utils import fetch_proteome


def main():
    human_prot = [
        (x.name, str(x.seq)) for x in SeqIO.parse(
            'data/uniprotkb_taxonomy_id_9605_AND_reviewed_2025_01_06.fasta', 'fasta',
        )
    ]

    ip_search_space = set()
    tryptic_search_space = set()

    for prot in human_prot:
        print(prot[0])
        prot_seq = prot[1]

        for length in range(8, 15):
            for idx in range(len(prot_seq)-(length-1)):
                ip_search_space.add(prot_seq[idx:idx+length])
        
        pattern = ['K', 'R']
        regex = re.compile(r'\b(' + '|'.join(pattern) + r')\b')
        kr_positions = (
            [-1] + sorted(
                [m.start() for m in re.finditer('R', prot_seq)] +
                [m.start() for m in re.finditer('K', prot_seq)]
            ) + [len(prot_seq)]
        )
        n_krs = len(kr_positions)

        for idx, kr_pos in enumerate(kr_positions[:-1]):
            start_pos = kr_pos + 1
            end_pos = kr_positions[idx+1] + 1
            if (end_pos - start_pos >= 7) & (end_pos - start_pos <= 30):
                tryptic_search_space.add(prot_seq[kr_pos+1:end_pos])
            if idx < n_krs - 2:
                end_pos = kr_positions[idx+2] + 1
                if (end_pos - start_pos >= 7) & (end_pos - start_pos <= 30):
                    tryptic_search_space.add(prot_seq[kr_pos+1:end_pos])
            if idx < n_krs - 3:
                end_pos = kr_positions[idx+3] + 1
                if (end_pos - start_pos >= 7) & (end_pos - start_pos <= 30):
                    tryptic_search_space.add(prot_seq[kr_pos+1:end_pos])

    print(len(tryptic_search_space))
    print(len(ip_search_space))

if __name__ == '__main__':
    main()
