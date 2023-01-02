import requests
import json

alleles = [
    'hla_b_27_05',
    'hla_b_35_01',
    'hla_b_08_01',
    'hla_b_57_01',
    'hla_b_27_09',
    'hla_b_15_02',
    'hla_a_24_02',
    'hla_b_07_02',
    'hla_e_01_03',
    'hla_b_44_02',
    'hla_b_15_01'
]



for allele in alleles:
    print(allele)
    url = f'https://datasette.histo.fyi/core.json?sql=select%0D%0A%20%20pdb_code%2C%20peptide_sequence%0D%0Afrom%0D%0A%20%20core%0D%0Awhere%0D%0A%20%20peptide_length%20%3D%209%0D%0A%20%20and%20peptide_features%20%3D%20%27correct_sequence_and_length%27%0D%0A%20%20and%20complex_type%20%3D%20%27class_i_with_peptide%27%0D%0A%20%20and%20allele_slug%20%3D%20%27{allele}%27%0D%0A%20%20group%20by%20peptide_sequence'
    r = requests.get(url)
    rows = r.json()['rows']
    peptides = []
    for row in rows:
        thisrow = {'pdb_code':row[0],'peptide_sequence':row[1]}
        peptides.append(thisrow)
    with open(f'sequences/{allele}_peptides.json', 'w') as peptides_file:
        json.dump(peptides, peptides_file)






