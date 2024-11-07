import os

from langchain import agents
from langchain.base_language import BaseLanguageModel

from chatmoldata.tools.molecule import *
from chatmoldata.tools.rdkit import SMILES2Weight,MolSimilarity
from chatmoldata.tools.database import Target2ID, SQLtext2CSV
from chatmoldata.tools.dataframe import CSV2Structre_image
from chatmoldata.tools.search import *

def make_tools(llm: BaseLanguageModel, api_keys: dict = {}, verbose=True):
    
    openai_api_key = api_keys.get("OPENAI_API_KEY") or os.getenv("OPENAI_API_KEY")

    all_tools = []
  

    all_tools += [
        MolSimilarity(),
        SMILES2Weight(),
        Target2ID(),
        Sub_search(),
        Subsearch_pre(),
        Text2SQL2CSV(),
        CSV2Structre_image(),
        CX_Preprocess(),
        CX_prediction(),
        CX_distribution(),
    ]
   

    return all_tools
