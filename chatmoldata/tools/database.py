#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 10:34:25 2024

@author: user
"""

from langchain_community.utilities import SQLDatabase
#from sqlalchemy import create_engine
from langchain.tools import BaseTool
import mysql.connector
from langchain.chains import create_sql_query_chain
from langchain.chat_models import ChatOpenAI
from langchain_community.vectorstores import FAISS
from langchain_core.example_selectors import SemanticSimilarityExampleSelector
from langchain_community.embeddings import OpenAIEmbeddings
from langchain_core.prompts import FewShotPromptTemplate, PromptTemplate
import re
import os

def extract_quoted_strings_from_tuple(row):
    row_str = str(row)
    
    quoted_strings = re.findall(r"'(.*?)'", row_str)
    combined_string = ' '.join(quoted_strings)
    return combined_string

connection = mysql.connector.connect(
    host="localhost",
    # port='7908',
    user="root",
    password="pass5959",
    # auth_plugin='mysql_native_password',
    database="chembl24"
)
cursor = connection.cursor()

db = SQLDatabase.from_uri('mysql+mysqlconnector://root:pass5959@localhost/chembl24',
                          include_tables=["compound_structures","activities","assays",
                                          "target_dictionary","target_components",
                                          "component_synonyms"],
                          sample_rows_in_table_info=1)
# print(db.dialect)
# print(db.get_usable_table_names())

openai_api_key = api_keys.get("OPENAI_API_KEY") or os.getenv("OPENAI_API_KEY")

llm = ChatOpenAI(model="gpt-4", temperature=0.0, openai_api_key = openai_api_key,streaming=False)

examples = [
    {'input': 'Retrieve compound smiles and IC50 for the target CHEMBL413, the IC50 < 200nM.',
     'query': '''SELECT cs.canonical_smiles,
     act.standard_value AS IC50 
     FROM target_dictionary td 
     JOIN assays a ON td.tid = a.tid 
     JOIN activities act ON a.assay_id = act.assay_id 
     JOIN compound_structures cs ON act.molregno = cs.molregno 
     WHERE act.standard_relation = '=' 
     AND act.standard_type IN ('IC50') 
     AND act.standard_units = 'nM' 
     AND act.standard_value < 200 
     AND td.chembl_id = 'CHEMBL413';'''},
     {'input':'''Retrieve compound information which are selective to the target CHEMBL5061 over 
      the target CHEMBL413, the IC50 for the former target is < 50nM, and IC50 for the latter target
      is > 200nM, the outputs includes smiles and two IC50s.'''
      ,'query':'''SELECT cs.canonical_smiles,
        act1.standard_value AS IC50_1,
        act2.standard_value AS IC50_2
        FROM compound_structures cs
        JOIN activities act1 ON cs.molregno = act1.molregno
        JOIN assays a1 ON act1.assay_id = a1.assay_id
        JOIN target_dictionary td1 ON a1.tid = td1.tid
        LEFT JOIN activities act2 ON cs.molregno = act2.molregno
        LEFT JOIN assays a2 ON act2.assay_id = a2.assay_id
        LEFT JOIN target_dictionary td2 ON a2.tid = td2.tid
        WHERE td1.chembl_id = 'CHEMBL5061'
        AND act1.standard_relation = '='
        AND act1.standard_type = 'IC50'
        AND act1.standard_units = 'nM'
        AND act1.standard_value < 50
        AND td2.chembl_id = 'CHEMBL413'
        AND act2.standard_relation = '='
        AND act2.standard_type = 'IC50'
        AND act2.standard_units = 'nM'
        AND act2.standard_value > 200;
      '''},
      {'input':'''Retrieve the chembl_id for the target 'MAT2A'.''',
       'query':'''SELECT target_dictionary.chembl_id
        FROM component_synonyms
        JOIN target_components ON component_synonyms.component_id = target_components.component_id
        JOIN target_dictionary ON target_components.tid = target_dictionary.tid
        WHERE component_synonyms.component_synonym = 'MAT2A'
        AND component_synonyms.syn_type = "GENE_SYMBOL"
        AND target_dictionary.target_type = "SINGLE PROTEIN"
        AND target_dictionary.organism = "Homo sapiens"
        LIMIT 1;
       '''
       },
      {'input':'''Search the compound "CC1=NC(C)=C(NC(C)=O)C=N1", show the activities in different targets (IC50 lower than 5000nM ).''',
       'query':'''SELECT act.standard_value,td.chembl_id,
       csy.component_synonym 
       FROM compound_structures cs 
       JOIN activities act ON act.molregno = cs.molregno 
       JOIN assays ays ON act.assay_id = ays.assay_id 
       JOIN target_dictionary td ON td.tid = ays.tid 
       JOIN target_components tc ON tc.tid = ays.tid 
       JOIN component_synonyms csy ON csy.component_id = tc.component_id 
       WHERE cs.canonical_smiles = 'COc1ccc(NC(=S)NC(=O)c2cn(-c3ccccc3)nc2-c2ccc(OC)cc2)cc1' 
       AND csy.syn_type = "GENE_SYMBOL" 
       AND act.standard_relation = '=' 
       AND act.standard_type IN ('IC50') 
       AND act.standard_units = 'nM' 
       AND act.standard_value < 5000;
       '''
       }
    ]

example_selector = SemanticSimilarityExampleSelector.from_examples(
    examples,
    OpenAIEmbeddings(api_key="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"),
    FAISS,
    k=5,
    input_keys=["input"],
)

example_prompt = PromptTemplate.from_template("User input: {input}\nSQL query: {query}")
prompt = FewShotPromptTemplate(
    example_selector=example_selector,
    example_prompt=example_prompt,
    prefix="You are a MySQL expert. Given an input question, create a syntactically correct SQLite query to run. Unless otherwise specificed, do not return more than {top_k} rows.\n\nHere is the relevant table info: {table_info}\n\nBelow are a number of examples of questions and their corresponding SQL queries.",
    suffix="User input: {input}\nSQL query: ",
    input_variables=["input", "top_k", "table_info"],
)

chain = create_sql_query_chain(llm, db, prompt)

class Target2ID(BaseTool):
    name = "Target2ID"
    description = (
        "Input should be the original texts of the task, returns text with the content of target is replaced by the chembl ID."
    )

    def __init__(self):
        super(Target2ID, self).__init__()

    def _run(self, texts: str) -> str:
        # print(path_file)
        
        # print('TD',texts)
        
        pattern = r'"(.*?)"'
        matches = re.findall(pattern, texts)
        
        # print("TD2",matches)
        
        list1 = []
        for ele in matches:
            task = "Retrieve the chembl_id for the target " + ele
            # input1 = '"question"' +": " + task
            # print("start invoking", input1)
            dict1 = {"question": task}
            # print(dict1)
            response = chain.invoke(dict1)
            cursor.execute(response)
            myresult = cursor.fetchone()
            out = extract_quoted_strings_from_tuple(myresult)
            if out is None:
                print(f"Error, unsuccessful replacement for {ele}.")
                
            else:
                
                result = f"replacing {ele} with {out}." 
                print(result)
            list1.append(out)
            
        # list_iterator = iter(list1)
        
        # for id in list1:
        #     replaced_text = re.sub(pattern, i, texts)
        
        # def replace_match(match):

        #     return f'"{next(list_iterator)}"'
            
        try:
            for item in list1:
                replaced_text = re.sub(pattern, item, texts, count = 1)
                texts = replaced_text
            # replaced_text = re.sub(pattern, next(list_iterator), texts)
            # print("result3: ", replaced_text)
            return replaced_text
        
        except:
            replaced_text = 'The task execution is completed, please check the error in this step.'
            print(replaced_text)
        return replaced_text
        
    async def _arun(self, path_file: str) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError()
        
class SQLtext2CSV(BaseTool):
    name = "SQLtext2CSV"
    description = (
        "Input should be the strings of Observation, returns the CSV file after excecuting SQL task."
    )

    def __init__(self):
        super(SQLtext2CSV, self).__init__()

    def _run(self, texts: str) -> str:
        # print(path_file)
        # print("kanzheli",texts)
        # input1 = '"question"' +": " + texts
        # print("start invoking", input1)
        
        dict1 = {"question": texts}
        # print(dict1)
        response = chain.invoke(dict1)
        
        # print(response)
        
        try:
            cursor.execute(response)
    
            column_names = [i[0] for i in cursor.description]
            rows = cursor.fetchall()

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
            # print(new_file)
            import csv
            with open(new_file, mode='w', newline='', encoding='utf-8') as file:
                csv_writer = csv.writer(file)
    
                csv_writer.writerow(column_names)
    
                csv_writer.writerows(rows)
                
            cursor.close()
            connection.close()
            result = f"Query results exported to {new_file} successfully, structural images need to be created."
            
    
        except:
            result = 'The task execution is completed, please check the error in this step.'

        return result
        
    async def _arun(self, path_file: str) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError()
