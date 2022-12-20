from pymol import cmd


amino_acids = ['V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C']

print (amino_acids)

motifs = [['2','C'],['2','5'],['2','4'],['4','C'],['5','C'],['2','5','C'],['2','4','C']]

motifs = [['2','C']]

peptide_length = 9


def make_peptide(sequence:str):
    print (sequence)
    cmd.fab(sequence, sequence.lower())
    cmd.save(f'inputs/scan_peptides/with_hydrogens/{sequence.lower()}_raw.pdb')
    cmd.remove('hydro')
    cmd.save(f'inputs/scan_peptides/without_hydrogens/{sequence.lower()}_raw.pdb')
    cmd.delete('all')


base_peptide = peptide_length*'A'
print (base_peptide)

peptide_set = {}

for motif in motifs:
    print (motif)
    this_motif = []
    for position in motif:
        if position == 'C':
            position = peptide_length
        this_motif.append(int(position))

    print (this_motif)        

    i = 0

    motif_peptide = ''
    while i < peptide_length:
        i += 1
        if i in this_motif:
            motif_peptide += str(i)
        else:
            motif_peptide += '-'

    print (motif_peptide)

    possibles = {}

    j = 0
    for position in this_motif:
        this_position = str(position)
        possibles[j] = []
        if j == 0:
            for amino_acid in amino_acids:
                if amino_acid != 'A':
                    this_possible = motif_peptide.replace(this_position, amino_acid)
                    possibles[j].append(this_possible)
        else:
            this_possibles = possibles[j - 1]
            for this_possible in this_possibles:
                for amino_acid in amino_acids:
                    if amino_acid != 'A':
                        new_possible = this_possible.replace(this_position, amino_acid)
                        possibles[j].append(new_possible)
        j += 1
        print (position)

    for peptide in possibles[len(possibles) - 1]:
        peptide = peptide.replace('-','A')
        print (peptide)
        make_peptide(peptide)
    




    
