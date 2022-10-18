""" Compute FP """
import os
import pandas as pd
from fingerprint_calculation.comp_fp_BV import *
from rdkit import DataStructs


class CompFP:
    def __init__(self, input_file: str, fp_name: str):
        self.root = os.getcwd()
        self.data = pd.read_csv(f"{self.root}/data/{input_file}")
        self.input_file = input_file
        print(self.data.head(2))
        self.fp_name = fp_name
        self.func_fp = self.functions_name[fp_name]

    @property
    def functions_name(self):
        functions_fp = {
            "ecfp4": morgan2,
            "ecfp6": morgan3
        }
        return functions_fp

    def fp_array(self):
        smiles = self.data.SMILES.to_list()
        fp = self.func_fp(smiles)
        output = []
        for f in fp:
            arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(f, arr)
            output.append(arr)
        exp_fp = np.asarray(output)
        return exp_fp

    def write_output(self):
        data = self.fp_array()
        explicit_vect = pd.DataFrame(
            data=data,
            columns=[str(i) for i in range(data.shape[1])],
            # index=[i for i in range(self.y.shape[0])],
        )
        final = pd.concat([self.data, explicit_vect], axis=1)
        print(final.head(5))
        output_file = f"{self.input_file.replace('.csv','')}_{self.fp_name}.csv"
        print(final.shape)
        # save output
        final.to_csv(f'{self.root}/data/{output_file}', index=False)


a = CompFP("reference_database.csv", "ecfp6")
a.write_output()
