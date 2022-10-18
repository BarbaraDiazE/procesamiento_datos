""" Compute FP """
import os
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

from compute_fp.com_fp_BV import morgan2, morgan3, topological, maccskeys, atom

# root = {"root": "/home/barbara/Documents/DIFACQUIM/PPI_classifier/phase-1/Databases/",
#         "root_Info" : "/home/barbara/Documents/DIFACQUIM/PPI_classifier/phase-1/Supervised_Models/SVM/Info",
#         "root_ROC" : "/home/barbara/Documents/DIFACQUIM/PPI_classifier/phase-1/Supervised_Models/SVM/ROC"}
root = {"root": os.getcwd()}


class CompFP:
    def __init__(self, input_file):
        self.Data = pd.read_csv(
            str(root["root"]) + str(input_file), index_col="Unnamed: 0"
        )
        print(self.Data.head())

    def fp_array(self, fp_func):
        smiles = self.Data.SMILES.to_list()
        fp = fp_func(smiles)
        output = []
        for f in fp:
            arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(f, arr)
            output.append(arr)
        self.y = np.asarray(output)
        print(self.y.shape)
        # print(self.y)

    def write_output(self, output_file):
        explicit_vect = pd.DataFrame(
            data=self.y,
            columns=[str(i) for i in range(self.y.shape[1])],
            index=[i for i in range(self.y.shape[0])],
        )
        print(explicit_vect.head(3))
        print("column names", explicit_vect.columns)
        print(self.Data.shape, explicit_vect.shape)
        frames = [self.Data, explicit_vect]
        Final = pd.concat([self.Data, explicit_vect], axis=1)
        print(Final.head(5))
        # Final
        Final.to_csv(f'{"/home/babs/Desktop/Data_Phase_2/final_files/"}{output_file}')


a = CompFP("dataset_p2.csv")
a.fp_array(atom)
a.write_output("dataset_atompairs.csv")
