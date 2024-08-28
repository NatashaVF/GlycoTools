# GlycoTools
Useful tools for incorporating ML in carbohydrate chemistry. Currently consisting of the GlycoPredicter which can be used for predicting anomeric selectivitites. 

## Installation

To install the project using conda, clone the repository, create a conda environment, and install the packages:

1. `git clone https://github.com/NatashaVF/GlycoTools.git`
2. `cd GlycoTools`
3. `conda env create -f environment.yml`
4. `conda activate glycotools`
4. `pip install -e .`

## Instruction
The GlycoPredicter can either be used in the terminal in the following way:


´cd GlycoPredict
GlycoPredict --input_path 'data/my_glycosylation.csv' --save_name 'my_glycosylation'´

or by using the `my_predictions.ipynb` jypyter notebook provided in the `/GlycoPredicter/` directory.

The csv file `my_glycosylation.csv` in the `GlycoPredicter/data/` directory serves as an illustration of the input dataformat.

TIP: ReactionSMILES can be created by drawing the reaction in ChemDraw and Copiyng As SMILES.

The predictions can be found in the `GlycoPredicter/results/` directory under the specified `save_name`.
