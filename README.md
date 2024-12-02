# ChatMolData

Implementation of the Paper "[ChatMolData: a Multimodal Agent for Automatic Molecular Data Processing](https://chemrxiv.org/engage/chemrxiv/article-details/67320b507be152b1d0bcf7f3)" (upload chemrxiv) by Yi Yu et al.. We assumed that the ChatMolData will bridge the gap between chemical experimenters and algorithm developers. 

<img src="example/TOC.png" width="100%" height="100%">

## Install via Anaconda
Create a new envioronment:
```bash
cd chatmoldata
conda env create -f environment.yml
conda activate chatmoldata
```
Install the chatmoldata package
```bash
python setup.py install
```

## Getting start
First set up your API keys in your environment.
```
export OPENAI_API_KEY=your-openai-api-key
```

In a Python:
```python
from chatmoldata.agents import ChatMolData

CMD = ChatMolData(model="gpt-4", temp=0.1)
CMD.run('''For ./data/mol_smiles.csv, calculate the multiple proper-tiesofthe molecules and the plot histogram showing distribution ofdifferent properties, properties inelude MW, ALogP, tPSA and QED.''') 
```
The output is:
```bash
Final Answer: The multiple properties of the molecules in the provided CSV file have been calculated and histograms showing the distributions of these properties have been plotted. The final CSV file with the calculated properties is 'mol smiles cleanedpred.csv`
```
<img src="example/output3.png" width="40%" height="40%">
