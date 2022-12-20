from typing import Dict, List

from pymol import cmd
from biopandas.pdb import PandasPdb


import os
import subprocess
import datetime
import json



abd_pdb_code = '1hhk'
exhaustiveness = 50
folder_name = 'predictions/box-test-35-15-20'
abd_filename = f'structures/{abd_pdb_code}_1_abd.pdb'




def get_gnina_info() -> Dict:
    result = subprocess.run(['tools/gnina','--version'], stdout=subprocess.PIPE)
    version_string = result.stdout.decode('utf-8').split('  ')

    gnina_version = version_string[0].split(' ')[1].strip()
    gnina_head = version_string[0].split(' ')[2].strip()
    gnina_build_date = version_string[1].replace('\n','').strip()

    return {
        'version':gnina_version,
        'head':gnina_head,
        'build_date': gnina_build_date
    }


def fetch_native_peptide_structure(pdb_code:str, abd_pdb_code:str) -> str:
    peptide_structure_url = f'https://coordinates.histo.fyi/structures/downloads/class_i/without_solvent/{pdb_code}_1_peptide.pdb'
    peptide_structure_path = f'structures/{pdb_code}_1_peptide'

    peptide_structure_pdb_file = f'{peptide_structure_path}.pdb'
    peptide_structure_sdf_file = f'{peptide_structure_path}_min.sdf'

    if not os.path.exists(peptide_structure_pdb_file):
        curl_command = f'curl {peptide_structure_url} --output {peptide_structure_pdb_file}'
        os.system(curl_command)

    if not os.path.exists(peptide_structure_sdf_file):
        gnina_command = f'tools/gnina --minimize -r structures/{abd_pdb_code}_1_abd.pdb -l {peptide_structure_pdb_file} -o {peptide_structure_sdf_file}'
        os.system(gnina_command)

    return peptide_structure_pdb_file


def generate_statistics(frame):
    statistics = {}
    raw_statistics = frame.split('END')[1]
    for line in raw_statistics.splitlines():
        if '>' in line:
            key = line.replace('>','').replace('<','')
        elif len(line) > 0 and '$' not in line:
            value = line
            if key is not None:
                statistics[key] = value
                key = None

    return statistics


def get_structure_name(filename:str):
    if '/' in filename:
        filepath_parts = filename.split('/')
        structure_name = filepath_parts[len(filepath_parts) - 1]
    else:
        structure_name = filename
    if '.' in structure_name:
        structure_name = structure_name.split('.')[0]
    return structure_name


def align_peptide(peptide_structure_filename:str, frame_filename:str) -> Dict:

    peptide_structure = get_structure_name(peptide_structure_filename)
    frame = get_structure_name(frame_filename)

    print(peptide_structure)
    print(frame)

    cmd.load(peptide_structure_filename, quiet=1)
    cmd.load(frame_filename, quiet=1)

    align = cmd.align(peptide_structure, frame, cycles=5, cutoff=2.0, mobile_state=-1, target_state=-1)

    alignment = {
        'rmsd':align[0],
        'atoms':align[1],
        'total_atoms':align[4]
    }
    return alignment



box = {
    'centre_of_mass':{
        'x': -43.190,
        'y': 57.376,
        'z':63.679
    },
    'sides': {
        'x':35,
        'y':15,
        'z':20
    }
}


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


i = 0
for peptide in peptides:

    if i < 10:
        peptide_filename = f'inputs/peptides/without_hydrogens/{peptide.lower()}_raw.pdb'


        output_filepath = f'{folder_name}/{peptide.lower()}_{abd_pdb_code}_ex{exhaustiveness}'
        
        
        docked_filename = f'{output_filepath}.sdf'
        stats_filename = f'{output_filepath}.stats.json'
        frames_filename = f'{output_filepath}.frames.json'
        params_filename = f'{output_filepath}.params.json'



        pdb_code = peptides[peptide]['pdb_code']

        peptide_structure_filename = fetch_native_peptide_structure(pdb_code, abd_pdb_code)
        

        gnina_command = f'tools/gnina -r {abd_filename} -l {peptide_filename} --center_x {box["centre_of_mass"]["x"]} --center_y {box["centre_of_mass"]["y"]} --center_z {box["centre_of_mass"]["z"]} --size_x {box["sides"]["x"]} --size_y {box["sides"]["y"]} --size_z {box["sides"]["z"]} -o {docked_filename} --exhaustiveness {exhaustiveness} --seed 0'
                
        print (f'{peptide} docked to {abd_pdb_code}')

        params = {
            'abd_structure':abd_pdb_code,
            'peptide':peptide,
            'pdb_code':pdb_code,
            'gnina_start_time': datetime.datetime.now().isoformat(),
            'box': box,
            'gnina':get_gnina_info(),
            'command':gnina_command
        }

        stats = {}

        #os.system(gnina_command)

        print (f'Native peptide - {peptide_structure_filename}')
        print (f'Docked peptide - {docked_filename}')


        split_token = '$$$$'
        

        with open(docked_filename) as file:
            sdf_file = file.read()

        frames = sdf_file.split(f'{split_token}\n')

        j = 1

        best_rmsd = 100
        for frame in frames:
            if len(frame) > 0:
                print (f'Frame {j}')
                frame_filename = f'{output_filepath}_{j}.sdf'
                with open(frame_filename, 'w') as frame_file:
                    frame_file.write(frame)
                
                frame_pdb_filename = f'{output_filepath}_{j}.pdb'
                rmsd_filename = f'{output_filepath}_{j}_rmsd.txt'

                obabel_command = f'obabel {frame_filename} -O {frame_pdb_filename} -d'
                os.system(obabel_command)

                frame_stats = generate_statistics(frame)

                crystal = PandasPdb().read_pdb(peptide_structure_filename)
                docked = PandasPdb().read_pdb(frame_pdb_filename)

                print (len(crystal.df['ATOM']))
                print (len(docked.df['ATOM']))
                
                
#                alignment = align_peptide(peptide_structure_filename, frame_filename)

                stats[str(j)] = frame_stats
                #stats[str(j)]['alignment'] = alignment

                #if alignment['rmsd'] < best_rmsd:
                #    stats['best_pose'] = str(j) 
                #    stats['best_rmsd'] = alignment['rmsd']
                #    best_rmsd = alignment['rmsd'] 

                j += 1


        params['gnina_end_time'] = datetime.datetime.now().isoformat()

        #print (stats['best_pose'])
        #print (stats['best_rmsd'])


        with open(stats_filename,'w') as stats_file:
            json.dump(stats, stats_file)

        with open(frame_filename,'w') as frames_file:
            json.dump(frames, frames_file)

        with open(params_filename,'w') as params_file:
            json.dump(params, params_file)
    i += 1

