import os

with open('sequences/sequences.txt','r') as input_file:
    peptides = [peptide.strip() for peptide in input_file.readlines()]


structures = ['1hhk','2bst','3d18','5iue']


for peptide in peptides:
    folder_name = f'predictions/gnina_dock_test/'

    exhaustiveness = 50
    i = 0
    for structure in structures:
        print (peptide)
        print (structure)
        abd_filename = f'structures/{structure}_1_abd.pdb'
        peptide_filename = f'inputs/{peptide.lower()}_raw.pdb'
        autobox_filename = f'structures/{structure}_1_peptide.pdb'
        docked_filename = f'{folder_name}/{peptide.lower()}_{structure}_ex{exhaustiveness}.sdf'
        gnina_command = f'tools/gnina -r {abd_filename} -l {peptide_filename} --autobox_ligand {autobox_filename} -o {docked_filename} --exhaustiveness {exhaustiveness} --seed 0'
            
        os.system(gnina_command)