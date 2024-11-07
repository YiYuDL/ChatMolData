#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 22:44:17 2024

@author: user
"""

import os
from langchain.tools import BaseTool
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdMolDescriptors
from chemcrow.utils import *
from rdkit.Chem.SaltRemover import SaltRemover
from pydantic import ValidationError
import pandas as pd
from rdkit.Chem import Descriptors, QED
from rdkit.Chem.rdMolDescriptors import CalcTPSA
from rdkit.Chem.Crippen import MolLogP
import matplotlib.pyplot as plt
from rdkit.Chem import Draw
from openpyxl import Workbook
from openpyxl.drawing.image import Image
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl.styles import Alignment
from io import BytesIO
import re

def clean(a):
        try:
            mol = Chem.MolFromSmiles(a)
            Chem.RemoveHs(mol)
            remover = SaltRemover()
            mol1 = remover.StripMol(mol)
            aa = Chem.MolToSmiles(mol1)
            return aa
        except:
            return None

def read_file(file):
    try:
        
        df_smiles = pd.read_csv(file)
        return df_smiles
        
    except:
        df_smiles = pd.read_excel(file)
        
        return df_smiles
        
    else:
        return 'Error, please check if the file is csv(xlsx) or include smiles strings'

class CMD_Preprocess(BaseTool):
    name = "cmd_preprocess"
    description = (
        "Input csv or xlsx file within molecular SMILES data, returns preprocessed csv file."
    )

    def __init__(self):
        super(CMD_Preprocess, self).__init__()

    def _run(self, path_file: str) -> str:
        try:
            #print("HERE",path_file,"THERE")
            # print(path_file)
            df_smiles = read_file(path_file)
            # print(len(df_smiles))
            # print(df_smiles)
            
            aa = df_smiles.columns.tolist()
            
            name = aa[0]
            df1 = df_smiles[name].map(clean)
            df1 = df1.dropna()
            # print(df1[4:10])
            # print(len(df1))
           
            nn = len(df_smiles) - len(df1)
            
            file2 = "prerocessed_" + path_file
            print
            if df1.empty:
                return "Error, no data left."
            elif nn == 0 :
                df1.to_csv(file2,index = False)
                print("The csv has been preprocessed, none of molecular smiles is invalid")
                return f"None of molecular smiles is invalid, the updated file is {file2}."
            else:
                df1.to_csv(file2,index = False)
                print(f"The csv has been preprocessed, {nn} of molecular smiles is invalid and the file has been updated")
                return f"{nn} of molecular smiles is invalid, the updated file is {file2}."
        except:
            print(path_file)
            
    async def _arun(self, path_file: str) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError()


class CMD_prediction(BaseTool):
    name = "cmd_prediction"
    description = "Input csv or xlsx file including molecular SMILES data, returns csv with multiple properties showing in new columns."

    def __init__(
        self,
    ):
        super(CX_prediction, self).__init__()

    def _run(self, path_file: str) -> str:
        # print('now'+path_file)
        print(path_file)
        df_smiles = read_file(path_file)
        
        results = []
        
        aa = df_smiles.columns.tolist()
        
        name = aa[0]

        try:
            for i, row in df_smiles.iterrows():
                smiles = row[name]
                mol = Chem.MolFromSmiles(smiles)
                
                mw = Descriptors.MolWt(mol)  # 分子量
                alogp = MolLogP(mol)  # ALogP
                tpsa = CalcTPSA(mol)  # tPSA
                qed = QED.qed(mol)  # QED
                #sa_score = sascorer.calculateScore(mol)  # SA Score
                
                results.append({
                    #"ID": i + 1,
                    "Smiles": smiles,
                    "MW": mw,
                    "ALogP": alogp,
                    "tPSA": tpsa,
                    "QED": qed,
                    #"SA Score": sa_score
                })
            
            df_results = pd.DataFrame(results)
            
            dot_index = path_file.rfind('.')
            new_filename = "predict_" + path_file
            print()
            df_results.to_csv(new_filename, index=False)
            return f"The updated file is {new_filename}."
            
        except:
            return 'Error, please check if the file is preprocessed'

    async def _arun(self, path_file: str) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError()

class CMD_distribution(BaseTool):
    name = "cmd_distribution"
    description = "Input csv file with molecular properties, returns histogram showing distributions of different properties, respectively."

    def __init__(
        self,
    ):
        super(CX_distribution, self).__init__()

    def _run(self, path_file: str) -> str:
        
        try:
            df1 = read_file(path_file)
            # print(df1.head())
            properties = ["MW", "ALogP", "tPSA", "QED", ]
            path1 = os.path.dirname(path_file)
            # print(f"result, {path1}")
            try:
                for prop in properties:
                    
                   
                    # print(df1[prop])
                    plt.figure(figsize=(10, 6))
                    plt.hist(df1[prop], bins=30, edgecolor='k', alpha=0.7)
                    plt.title(f'Distribution of {prop}')
                    plt.xlabel(prop)
                    plt.ylabel('Frequency')
                    plt.grid(True)
                    # print(f'{path1}/{prop}_frequency_distribution.png')
                    # plt.savefig(f'{path1}/{prop}_frequency_distribution.png')
                    plt.show()
                    
                return "Figures have been shown"
                    
            except:
                
                print(f"no column named {prop}")
        
        except:
            return "please check the if the file is what you plan to analyze."
        
       
    async def _arun(self, path_file: str) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError()
        
        
class Structure2Smiles(BaseTool):
    name = "Structure2Smiles"
    description = (
        "Input cif, mol, sdf, png or jpg file representing molecular structures, returns molecular SMILES."
    )

    def __init__(self):
        super(CX_preprocess, self).__init__()

    def _run(self, path_file: str) -> str:
       
        try:
            if path_file[-3:0] == 'png' or path_file[-3:0] == 'jpg':
                smiles = ''
                
            if path_file[-3:0] == 'cif':
                mol = Chem.MolFromMolCIFfile(file_path)
                smiles = Chem.MolToSmiles(mol)
                
            elif path_file[-3:0] == 'mol':
                mol = Chem.MolFromMOLfile(file_path)
                smiles = Chem.MolToSmiles(mol)
            
            elif path_file[-3:0] == 'sdf':
                mol = Chem.MolFromSDFfile(file_path)
                smiles = Chem.MolToSmiles(mol)
                
            else:
                smiles = ''
            return smiles
        except:
            print("The file is unsupported")
            return "The task execution is completed, please check the error in this step."
        
        if df1.empty:
            return "Error, no data left."
        else:
            df1.to_csv(path_file,index = False)
            return "The csv has been preprocessed, and the file has been updated"

    async def _arun(self, path_file: str) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError()
        
