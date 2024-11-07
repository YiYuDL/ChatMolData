from langchain.tools import BaseTool
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdMolDescriptors
from chatmoldata.utils.utils import *
from pydantic import ValidationError
from langchain_core.prompts import FewShotPromptTemplate, PromptTemplate


class MolSimilarity(BaseTool):
    name = "MolSimilarity"
    description = (
        "Input two molecule SMILES (separated by '.'), returns Tanimoto similarity."
    )

    def __init__(self):
        super(MolSimilarity, self).__init__()

    def _run(self, smiles_pair: str) -> str:
        smi_list = smiles_pair.split(".")
        if len(smi_list) != 2:
            return "Input error, please input two smiles strings separated by '.'"
        else:
            smiles1, smiles2 = smi_list

        similarity = tanimoto(smiles1, smiles2)

        if isinstance(similarity, str):
            return similarity

        sim_score = {
            0.9: "very similar",
            0.8: "similar",
            0.7: "somewhat similar",
            0.6: "not very similar",
            0: "not similar",
        }
        if similarity == 1:
            return "Error: Input Molecules Are Identical"
        else:
            val = sim_score[
                max(key for key in sim_score.keys() if key <= round(similarity, 1))
            ]
            message = f"The Tanimoto similarity between {smiles1} and {smiles2} is {round(similarity, 4)},\
            indicating that the two molecules are {val}."
        return message

    async def _arun(self, smiles_pair: str) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError()


class SMILES2Weight(BaseTool):
    name = "SMILES2Weight"
    description = "Input SMILES, returns molecular weight."

    def __init__(
        self,
    ):
        super(SMILES2Weight, self).__init__()

    def _run(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES string"
        mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
        return mol_weight

    async def _arun(self, smiles: str) -> str:
        """Use the tool asynchronously."""
        raise NotImplementedError()
