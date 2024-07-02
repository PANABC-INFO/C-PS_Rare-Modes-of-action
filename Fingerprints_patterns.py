#!/usr/bin/env python
# coding: utf-8

# RDKit library has to be installed

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
import os

# Define the path to the sturtural patterns file
patterns_structures_file = "example/substructure_patterns_example.csv"
patterns_structures_df = pd.read_csv(patterns_structures_file)

# Input the file names or will run the example file
set_flag_input = True

if not set_flag_input:
    # Input CSV file should have columns "MoleculeID" and "SMILES"
    input_file = "input_file_here"
    output_file = "output_file_here"
    molecule_data_df = None
    if not os.path.isfile(input_file):
        print(f"ERROR: File {input_file} not found, PLEASE CHECK YOUR INPUT FILES")
    else:
        molecule_data_df = pd.read_csv(input_file)
else:   
    # Example input CSV file with the columns "Molecule Name" and "SMILES"
    example_input = "example/example_HD_molecules.csv"
    molecule_data_df = pd.read_csv(example_input)

# Convert the smiles of structural patterns to molecule objects
show_molecules = True
pattern_molecules = [Chem.MolFromSmiles(smiles) for smiles in list(patterns_structures_df["SMILES"])]

# hash the SMILES of the structural patterns as header for the output
headers = list(patterns_structures_df["SMILES"])
headers.insert(0, "Patterns Name")

# Show the strucural patterns as molecules
if show_molecules:
    Draw.MolsToGridImage(pattern_molecules, subImgSize=(250, 250), molsPerRow=5, maxMols=100)

# Convert the screening molecules smiles to molecule objects
screening_molecules = [Chem.MolFromSmiles(smiles) for smiles in list(molecule_data_df["SMILES"])]
screening_molecule_names = [name for name in list(molecule_data_df["Molecule Name"])]

# Show the screening molecules as molecules 
if show_molecules:
    Draw.MolsToGridImage(screening_molecules, subImgSize=(250, 250), molsPerRow=5, legends=screening_molecule_names)

# Perform the substructure search of the strucutral patterns and convert the results into a bit-string
fingerprints_dict = {}
for idx, row in molecule_data_df.iterrows():
    molecule_id = row["Molecule Name"]
    molecules_smiles = row["SMILES"]
    mol = Chem.AddHs(Chem.MolFromSmiles(molecules_smiles))
    fingerprints_bits = [1 if mol.HasSubstructMatch(ring) else 0 for ring in pattern_molecules]
    fingerprints_bits.insert(0, molecule_id)
    fingerprints_dict[idx] = fingerprints_bits

# Store the resulted data into file
fingerprints_dict_df = pd.DataFrame.from_dict(fingerprints_dict, orient='index', columns=headers)
display(fingerprints_dict_df)
if set_flag_input:  
    fingerprints_dict_df.to_csv("NewData_HD_Molecules_01_six.csv", index=False) 
else:
    fingerprints_dict_df.to_csv(output_file, index=False)




