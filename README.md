## Flavor Analysis and Recognizion Transformer (FART)

## Description
FART is a tool to predict the taste of a molecule, based on its SMILES representation. The model utilizes a fine-tuned LLM, ChemBERTa, trained on SMILES....
The FART dataset combines several sources into a FAIR dataset.

## Installation
All .ipynb notebooks can be run on GoogleColab without any further modifications. 

## Overview of files

* the entire raw, curated and enriched dataset is found in /dataset
* the dataset split along sources can be found in /dataset/Individual-Datasets
* Files for the extraction, curation and enrichment of the dataset are found in /dataset
* Files for the training of the random forest and transformer models are found in /models 

## Dataset

`FART_Data_Extraction.ipynb` extracts data from five different online sources and produces the dataset `fart_uncurated.csv`. 

`FART_Data_curation.ipynb` curates the extracted data by for example removing duplicates through standardized SMILES. This scripts produces the dataset `fart_curated.csv` which was used in the training of the machine learning models. 

`FART_dataset_enrichment.ipynb` can be optionally used to retrieve more features for molecules which are also listed on PubChem. This script produces the `fart_enriched.csv` dataset which additionally includes the columns `PubChemID`, `IUPAC Name`, `Molecular Formula`, `Molecular Weight`, `InChI` and `InChiKey`. 

## Random Forest Models 

All three random forest models were trained in `model/FART-Random-Forest-Models.ipynb` using `fart_uncurated.csv` which needs to be downloaded for the script to function. The path in Cell 4 needs to be adjusted for the pandas import to succeed. 

## Transformer Models

The transformer models were trained in `model/FART-Transformer-Models.ipynb`. The data is loaded using the hugging face api. For different pretrained models one needs to adjust the `model_checkpoint` parameter. To use a weighted loss fontion, one needs to use `trainer = CustomTrainer` instead of `trainer = Trainer`. To use augmentation on needs to set `augmentation = True`. The attention weights were extracted using `model/attention_extractor.ipynb`. Trained transformer models are available at https://huggingface.co/FartLabs.

