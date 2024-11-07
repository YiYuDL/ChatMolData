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
from chatmoldata.utils.utils import *
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


from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem import rdMolEnumerator
from IPython.display import display

def align_bundle_coords(bndl):
    ps = rdFMCS.MCSParameters()
    for m in bndl:
        Chem.SanitizeMol(m)
    mcs = rdFMCS.FindMCS(bndl,completeRingsOnly=True)
    q = Chem.MolFromSmarts(mcs.smartsString)
    rdDepictor.Compute2DCoords(q)
    for m in bndl:
        rdDepictor.GenerateDepictionMatching2DStructure(m,q)
        
def Substruc_Prepare(mol):
    qry = Chem.MolFromMolFile(mol)
    qry_bundle = rdMolEnumerator.Enumerate(qry)
    qry_addHs_bundle = []
    
    for molecule in qry_bundle:
        mol1 = Chem.AddHs(molecule, addCoords=True)
        qry_addHs_bundle.append(mol1)
    
    align_bundle_coords(qry_bundle)
    print("子结构列表：")
    ac = Draw.MolsToGridImage(qry_bundle, subImgSize=(200, 200), molsPerRow=len(qry_bundle))
    display(ac)
    return qry_addHs_bundle

def Substruc_Prepare2(mol):
    qry = Chem.MolFromMolBlock(mol)
    qry_bundle = rdMolEnumerator.Enumerate(qry)
    qry_addHs_bundle = []
    
    for molecule in qry_bundle:
        mol1 = Chem.AddHs(molecule, addCoords=True)
        qry_addHs_bundle.append(mol1)
    
    align_bundle_coords(qry_bundle)
    print("子结构列表：")
    ac = Draw.MolsToGridImage(qry_bundle, subImgSize=(200, 200), molsPerRow=len(qry_bundle))
    
    display(ac)
    return qry_addHs_bundle

def SubSearch(substruct_list, csv_file):
    mols_no_H = []
    mols_addHs = []
    
    # smiles_list = df
    ''' gai '''
    df1 = read_file(csv_file)
    print("df1 len: ",len(df1))
    smiles_columns = [col for col in df1.columns if re.search(r'(?i)smiles', col)][0]
    ic_columns = [col for col in df1.columns if re.search(r'(?i)IC50', col)][0]
    # df2 = df1.groupby(smiles_columns, as_index=False)[ic_columns].mean()
    df_smiles = df1[smiles_columns].tolist()
    
    for molecule in df_smiles:
        mol1 = Chem.MolFromSmiles(molecule)
        mols_no_H.append(mol1)
        mol2 = Chem.AddHs(mol1, addCoords=True)
        mols_addHs.append(mol2)
    print("分子数据集:")
    # ac = Draw.MolsToGridImage(mols_no_H,returnPNG=False,molsPerRow=3, subImgSize=(300, 250))
    # ac = Draw.MolsToGridImage(mols_no_H, maxMols = 200, subImgSize=(300, 250))
    # display(ac)
    
    matches = []
    matches2 = []
    matched_ats = []
    for one in substruct_list:
        for x in mols_addHs:
            match = x.GetSubstructMatch(one)
            if match:
                x_no_H = Chem.RemoveHs(x)
                smiles = Chem.MolToSmiles(x_no_H)
                new_match = [atom.GetIdx() for atom in x_no_H.GetAtoms() if atom.GetIdx() in match]
                matches.append(smiles)
                matches2.append(x_no_H)
                matched_ats.append(new_match)
    
    ic = []
    for match in matches:
        matching_rows = df1[df1[smiles_columns] == match]
        ic.extend(matching_rows[ic_columns].tolist())

    ic = list(map(lambda x: f"{x:.3g} nM",ic))
    
    if matches:
        # print("分子数据集:")maxMols = 50
        ae = Draw.MolsToGridImage(matches2,molsPerRow=3,returnPNG=False, legends = ic, 
                                  highlightAtomLists=matched_ats, subImgSize=(300, 250))
        # ac = Draw.MolsToGridImage(qry_bundle, subImgSize=(200, 200), molsPerRow=len(qry_bundle))
        # display(ac)
        print("搜索结果：")
        display(ae)
    else:
        ae = []
        print("no matches")
    return matches, ae
    

# 读取CSV文件

class Subsearch_pre(BaseTool):
    name = "subsearch_pre"
    description = (
        "the input should be the original strings of the task, returns strings."
    )

    def __init__(self):
        super(Subsearch_pre, self).__init__()

    def _run(self, input1: str) -> str:
        try:
            
            # df1 = read_file('output_39.csv')
            # print("df1 len: ",len(df1))
            # smiles_columns = [col for col in df1.columns if re.search(r'(?i)smiles', col)][0]
            # ic_columns = [col for col in df1.columns if re.search(r'(?i)IC50', col)][0]
            # df2 = df1.groupby(smiles_columns, as_index=False)[ic_columns].mean()
            # print("df2 len: ",len(df2))
            
            
            pattern = r'((\S+\.csv)\s+and\s+(\S+\.mol))|((\S+\.mol)\s+and\s+(\S+\.csv))'
  
            match = re.search(pattern, input1)
    
            if match:
                
                if match.group(2) and match.group(3):
                    csv_file = match.group(2)
                    mol_file = match.group(3)
                else:
                    csv_file = match.group(6)
                    mol_file = match.group(5)
  
            else:
              print("the input is not proper")
            
            df1 = read_file(csv_file)
            print("df1 len: ",len(df1))
            smiles_columns = [col for col in df1.columns if re.search(r'(?i)smiles', col)][0]
            ic_columns = [col for col in df1.columns if re.search(r'(?i)IC50', col)][0]
            df2 = df1.groupby(smiles_columns, as_index=False)[ic_columns].mean()
            print("df2 len: ",len(df2))
            
            current__path = os.getcwd()
            base = 'output'
            extension = '.csv'
            counter=0
            
            while True:
                filename = f"{base}_{counter}{extension}"
                
                if not os.path.exists(filename):
                    break
                
                counter +=1
            # print()
            new_file = current__path + "/" + filename
            df2.to_csv(new_file,index=False)
            # print(new_file)
            # import csv
            # with open(new_file, mode='w', newline='', encoding='utf-8') as file:
            #     csv_writer = csv.writer(file)
    
            #     csv_writer.writerow(column_names)
    
            #     csv_writer.writerows(rows)
            
            # df2.to_csv(file2,index=False)
            # sub, image1 = Substruc_Prepare(mol_file)
            # display(image1)
            
            # print(df2.head())
            # print(sub)
            
            return f"{mol_file} and {new_file}"
        
            # if df1.empty:
            #     return "Error, no data left."
            # elif nn == 0 :
            #     df1.to_csv(file2,index = False)
            #     print("The csv has been preprocessed, none of molecular smiles is invalid")
            #     return f"None of molecular smiles is invalid, the updated file is {file2}."
            # else:
            #     df1.to_csv(file2,index = False)
            #     print(f"The csv has been preprocessed, {nn} of molecular smiles is invalid and the file has been updated")
            #     return f"{nn} of molecular smiles is invalid, the updated file is {file2}."
        except:
            return f"the input is {input1}"
            
    async def _arun(self, path_file: str) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError()
        
class Sub_search(BaseTool):
    name = "sub_search"
    description = (
        "the input should be the original strings of the task or Observation, returns a list and dataframe variables."
    )

    def __init__(self):
        super(Sub_search, self).__init__()
        
    def _run(self, input1) -> str:
        try:
            
            pattern = r'((\S+\.csv)\s+and\s+(\S+\.mol))|((\S+\.mol)\s+and\s+(\S+\.csv))'
  
            match = re.search(pattern, input1)
    
            if match:
                
                if match.group(2) and match.group(3):
                    csv_file = match.group(2)
                    mol_file = match.group(3)
                else:
                    csv_file = match.group(6)
                    mol_file = match.group(5)
  
            else:
                print("the input is not proper")
            
            sub_list = Substruc_Prepare(mol_file)
            
            matches, ae = SubSearch(sub_list, csv_file)
            print(f"result: {len(matches)}")
            # display(ac)
            # display(ae)
            # df1 = read_file('output_39.csv')
            # print("df1 len: ",len(df1))
            # smiles_columns = [col for col in df1.columns if re.search(r'(?i)smiles', col)][0]
            # ic_columns = [col for col in df1.columns if re.search(r'(?i)IC50', col)][0]
            # df2 = df1.groupby(smiles_columns, as_index=False)[ic_columns].mean()
            # print("df2 len: ",len(df2))
            
            # m1, a2, a3 = SubSearch(list1)
            # sub, image1 = Substruc_Prepare(mol_file)
            # display(a2)
            # display(a3)
            
            # print(df2.head())
            # print(sub)
            
            return 'substrcutre search is completed'
        
            # if df1.empty:
            #     return "Error, no data left."
            # elif nn == 0 :
            #     df1.to_csv(file2,index = False)
            #     print("The csv has been preprocessed, none of molecular smiles is invalid")
            #     return f"None of molecular smiles is invalid, the updated file is {file2}."
            # else:
            #     df1.to_csv(file2,index = False)
            #     print(f"The csv has been preprocessed, {nn} of molecular smiles is invalid and the file has been updated")
            #     return f"{nn} of molecular smiles is invalid, the updated file is {file2}."
        except:
            return f"the input is {input1}"
            
    async def _arun(self, path_file: str) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError()
