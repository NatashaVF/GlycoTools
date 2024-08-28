# GlycoTools
Useful tools for incorporating ML in carbohydrate chemistry. Currently consisting of the GlycoPredicter which can be used for predicting anomeric selectivitites. 

## Installation

To install the project, clone the repository and install the dependencies:

```bash
git clone https://github.com/NatashaVF/GlycoTools.git
cd GlycoTools
conda env create -f environment.yml
conda activate glycotools
pip install -e .

## Instruction
cd GlycoPredict
 GlycoPredict --input_path 'data/my_glycosylation.csv' --save_name 'my_glycosylation'
