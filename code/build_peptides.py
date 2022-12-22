from pymol import cmd
import os
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser, Select
from Bio.PDB.PDBIO import PDBIO



parser = PDBParser()


class NotDisorderedNotOxt(Select):
    def accept_atom(self, atom):
        return_value = False
        if not atom.get_name() == 'OXT':
            if not atom.is_disordered():
                return_value = True
        if atom.get_altloc() == "A":
            return_value = True
        return return_value


def build_peptide(peptide:str):
    peptide_file_name = f'inputs/peptides/without_hydrogens/{peptide.lower()}_raw.pdb'
    if not os.path.exists(peptide_file_name):
        cmd.fab(peptide, peptide.lower())
        cmd.save(f'inputs/peptides/with_hydrogens/{peptide.lower()}_raw.pdb')
        cmd.remove('hydro')
        cmd.save(peptide_file_name)
        cmd.delete('all')
    return peptide_file_name


def cleanup_native_peptide_structure(pbd_code:str):
    print (pdb_code)
    peptide_structure_pdb_file = f'structures/{pdb_code}_1_peptide.pdb'
    revised_peptide_structure_pdb_file = f'structures/{pdb_code}_1_peptide_clean.pdb'
    print (revised_peptide_structure_pdb_file)
    structure = parser.get_structure("peptide",peptide_structure_pdb_file)
    io = PDBIO()
    io.set_structure(structure)
    io.save(revised_peptide_structure_pdb_file, select=NotDisorderedNotOxt())
    return revised_peptide_structure_pdb_file


def fetch_native_peptide_structure(pdb_code:str) -> str:
    peptide_structure_pdb_file = f'structures/{pdb_code}_1_peptide.pdb'

    if not os.path.exists(peptide_structure_pdb_file):
        peptide_structure_url = f'https://coordinates.histo.fyi/structures/downloads/class_i/without_solvent/{pdb_code}_1_peptide.pdb'
        curl_command = f'curl {peptide_structure_url} --output {peptide_structure_pdb_file}'
        os.system(curl_command)
    return peptide_structure_pdb_file


def compare_peptide_atoms(built_peptide_filename:str, native_peptide_filename:str):
    crystal = PandasPdb().read_pdb(native_peptide_filename)
    built = PandasPdb().read_pdb(built_peptide_filename)
    crystal_atom_count = len(crystal.df['ATOM'])
    built_atom_count = len(built.df['ATOM'])
    print (crystal_atom_count)
    print (built_atom_count)
    if built_atom_count != crystal_atom_count:
        print ('UNEQUAL')
    else:
        print ('EQUAL')

peptides = {}


with open('sequences/sequences.txt','r') as input_file:
    for line in input_file.readlines():
        if len(line) > 0:
            pdb_code = line.split(',')[0].strip()
            peptide = line.split(',')[1].strip()
            peptides[peptide] = {
                'pdb_code':pdb_code,
                'sequence':peptide
            }


for peptide in peptides:
    print ('---------')
    print (peptide)
    pdb_code = peptides[peptide]['pdb_code']
    print (pdb_code)
    built_peptide_filename = build_peptide(peptide)
    native_peptide_filename = fetch_native_peptide_structure(pdb_code)
    revised_native_peptide_filename = cleanup_native_peptide_structure(pdb_code)
    print (peptides[peptide])
    print(revised_native_peptide_filename)
    compare_peptide_atoms(built_peptide_filename, revised_native_peptide_filename)
    