from pymol import cmd

with open('sequences/sequences.txt','r') as input_file:
    peptides = [line.split(',')[1].strip() for line in input_file.readlines() if len(line) > 0]


for peptide in peptides:
    print (peptide)
    cmd.fab(peptide, peptide.lower())
    cmd.save(f'inputs/peptides/with_hydrogens/{peptide.lower()}_raw.pdb')
    cmd.remove('hydro')
    cmd.save(f'inputs/peptides/without_hydrogens/{peptide.lower()}_raw.pdb')
    cmd.delete('all')