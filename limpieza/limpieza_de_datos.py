#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import SaltRemover
from molvs.standardize import Standardizer
from molvs.tautomer import TautomerCanonicalizer

"""
remover átomos no deseados
"""

# delete organometalics
def check_atoms(mol):
    """
    filter compounds that contains just desired atoms
    """
    desired_elements = {"H", "B", "C", "N", "O", "F", "P", "S", "Cl"}
    elements = set([i.GetSymbol() for i in mol.GetAtoms()])
    atom_difference = len(elements - desired_elements)
    return atom_difference


def remove_salts(mol):
    remover = SaltRemover.SaltRemover()
    res = remover.StripMol(mol)
    # return Chem.MolToSmiles(res)
    return res


def standarization(mol):
    s = Standardizer()
    s_mol = s.standardize(mol)
    return s_mol


def get_tautomer(mol):
    TC = TautomerCanonicalizer()
    return TC(mol)


"""
execution
"""


def clean_database(input_file, smiles_column, output_file):
    route = "/tmpu/jlmf_g/bide_a/AI_BD_OP/clean_ChEMBL/"
    # route = "/home/babs/Documents/DIFACQUIM/AI_BD_OP/clean_ChEMBL/"
    data = pd.read_csv(route + input_file, sep=",")
    data = data.reset_index()
    data = data.drop("index", axis=1)
    data = data.drop("standard_inchi", axis=1)
    data = data.drop("standard_inchi_key", axis=1)
    ##
    smi = data[smiles_column].tolist()

    #
    cleaned_smiles = list()
    index = list()
    for i in range(len(smi)):
        # print("smiles: ", smi[i])
        try:
            mol = Chem.MolFromSmiles(smi[i])
            if mol == None:
                continue
                # print("smile not valid")
            else:
                r = check_atoms(mol)
                if r == 0:
                    ns_mol = remove_salts(mol)
                    s_mol = standarization(ns_mol)
                    index.append(i)
                    cleaned_smiles.append(Chem.MolToSmiles(s_mol))
                else:
                    continue

        except:
            return "function not valid"

    cleaned_data = data.loc[index, :]
    cleaned_data = cleaned_data.reset_index()
    cleaned_data = cleaned_data.drop("index", axis=1)
    cleaned_data["cleaned_smiles"] = cleaned_smiles
    print("cleaned data columns: ", cleaned_data.columns)
    cleaned_data.to_csv(route + output_file, sep=",", index=True)
    print(
        f'{"initial compounds "}{data.shape[0]}{", "} {"final compounds "}{cleaned_data.shape[0]}'
    )
    return cleaned_data


cleaned_data = clean_database(
    "1_million_from_ChEMBL.csv", "canonical_smiles", "cleaned_from_1M_ChEMBL_stereo.csv"
)