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
# from DECIMER import predict_SMILES
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

class CX_Preprocess(BaseTool):
    name = "cx_preprocess"
    description = (
        "Input csv or xlsx file within molecular SMILES data, returns preprocessed csv file."
    )

    def __init__(self):
        super(CX_Preprocess, self).__init__()

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


class CX_prediction(BaseTool):
    name = "cx_prediction"
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

class CX_distribution(BaseTool):
    name = "cx_distribution"
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
        
# class CSV2Structre_image(BaseTool):
#     name = "CSV2Structre_image"
#     description = "Input is the csv file, draw the molecular image based smiles, then returns the xlsx file including the image."

#     def __init__(
#         self,
#     ):
#         super(CSV2Structre_image, self).__init__()

#     def _run(self, path_file: str) -> str:
        
#         try:
            
#             df = pd.read_csv(path_file)

#             wb = Workbook()
#             ws = wb.active
#             nn =  len(df.columns) + 1
#             pattern = re.compile(r'.*smiles.*', re.IGNORECASE)
#             smiles_columns = [col for col in df.columns if pattern.match(col)]

#             smiles_column = smiles_columns[0]
#             if len(smiles_columns) != 1:
#                 raise ValueError("no column name includes 'smiles' ")

#             for row in range(1, len(df) + 2):
#                 ws.row_dimensions[row].height = 250

#             center_alignment = Alignment(horizontal='center', vertical='center')

#             for c_idx, col_name in enumerate(df.columns):
#                 cell = ws.cell(row=1, column=c_idx + 1, value=col_name)
#                 # cell.alignment = center_alignment
                

#                 col_letter = chr(65 + c_idx)  
#                 if df.columns[c_idx] == smiles_column:
#                     ws.column_dimensions[col_letter].width = 50
#                 else:
#                     ws.column_dimensions[col_letter].width = 20

#             for r_idx, row in df.iterrows():
#                 for c_idx, value in enumerate(row):
#                     ws.cell(row=r_idx + 2, column=c_idx + 1, value=value)
#                     # cell.alignment = center_alignment

#             img_column_index = len(df.columns) + 1
#             ws.cell(row=1, column=img_column_index, value='structures')

#             for idx, smiles in enumerate(df[smiles_column], start=2):
#                 molecule = Chem.MolFromSmiles(smiles)
#                 img = Draw.MolToImage(molecule, size=(400, 400))
#                 buffered = BytesIO()
#                 img.save(buffered, format="PNG")
#                 buffered.seek(0)
                
#                 image = Image(buffered)
#                 image.width = 300  
#                 image.height = 300 
#                 image.anchor = ws.cell(row=idx, column=nn).coordinate  # 将图片插入到B列
#                 ws.add_image(image)

#             ws.column_dimensions[chr(65 + len(df.columns))].width = 50

#             ws.row_dimensions[1].height = 20

#             for row in ws.iter_rows(min_row=1, max_row=len(df) + 1, max_col=len(df.columns) + 2):
#                 for cell in row:
#                     cell.alignment = center_alignment
                    
#             filename = path_file[:-3] + 'xlsx'
#             wb.save(filename)
#             return filename
            
#         except:
#             return "The task execution is completed, please check the error in this step."
        
       
#     async def _arun(self, path_file: str) -> str:
#         """Use the tool asynchronously."""
#         raise NotImplementedError()
        
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
        
# import os
# from langchain.tools import BaseTool
# from rdkit import Chem, DataStructs
# from rdkit.Chem import AllChem, rdMolDescriptors
# from chemcrow.utils import *
# from DECIMER import predict_SMILES
# from rdkit.Chem.SaltRemover import SaltRemover
# from pydantic import ValidationError
# import pandas as pd
# from rdkit.Chem import Descriptors, QED
# from rdkit.Chem.rdMolDescriptors import CalcTPSA
# from rdkit.Chem.Crippen import MolLogP
# import matplotlib.pyplot as plt
# from rdkit.Chem import Draw
# from openpyxl import Workbook
# from openpyxl.drawing.image import Image
# from openpyxl.utils.dataframe import dataframe_to_rows
# from openpyxl.styles import Alignment
# from io import BytesIO
# import re   

'''    
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
    
    
    if matches:
        # print("分子数据集:")maxMols = 50
        ae = Draw.MolsToGridImage(matches2,molsPerRow=3,returnPNG=False,
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
'''