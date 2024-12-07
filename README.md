## Flavor Analysis and Recognizion Transformer (FART)

## Description
The Flavor Analysis and Recognition Transformer (FART) is a state-of-the-art machine learning model designed to predict molecular taste from chemical structures encoded as SMILES. Developed using the pre-trained foundation model ChemBERTa, FART leverages a transformer architecture to classify molecules across four key taste categories—sweet, bitter, sour, and umami—while also accommodating tasteless or undefined compounds. Unlike previous approaches that relied on binary classification, FART performs parallel multi-class predictions with an accuracy exceeding 91%, offering interpretability through gradient-based visualizations of molecular features. This novel approach facilitates the identification of key structural elements influencing taste, enabling applications in both flavor compound discovery and rational food design.

## Installation
All .ipynb notebooks can be run on GoogleColab without any further modifications. 

## Overview of files

* The entire raw, curated and enriched dataset is found in /dataset
* The source databases are found in dataset/individual-datasets
* Files for the extraction, curation and enrichment of the dataset are found in /dataset/scripts
* Files for the training of the tree-based and transformer models are found in /models
* Script for the generation of the interpretability heatmaps are found in /plots

## Dataset

The FART dataset is the largest publicly available collection of molecular tastants to date, comprising 15,025 curated entries derived from six independent sources. Each molecule is annotated with one or more taste labels (sweet, bitter, sour, umami, or undefined) following rigorous curation protocols. Duplicates were removed based on canonicalized SMILES, reducing overlap among entries. Additional enrichment with metadata such as PubChem ID, IUPAC name, molecular formula, and molecular weight was performed using the PubChem API. The dataset adheres to the FAIR principles, ensuring accessibility and reusability, and is hosted publicly to support further research. Its chemical diversity spans a molecular weight range centered at 374 Da ± 228, making it suitable for small molecule taste prediction.

`FART_Data_Extraction.ipynb` extracts data from five different online sources and produces the dataset `fart_uncurated.csv`. 

`FART_Data_curation.ipynb` curates the extracted data by for example removing duplicates through standardized SMILES. This scripts produces the dataset `fart_curated.csv` which was used in the training of the machine learning models. 

`FART_dataset_enrichment.ipynb` can be optionally used to retrieve more features for molecules which are also listed on PubChem. This script produces the `fart_enriched.csv` dataset which additionally includes the columns `PubChemID`, `IUPAC Name`, `Molecular Formula`, `Molecular Weight`, `InChI` and `InChiKey`. 

## Random Forest Models 

All three tree-based classifiers were trained in `model/Tree-Baseline-Models.ipynb`.

## Transformer Models

The transformer models were trained in `model/FART_Models.ipynb`. The data is loaded using the hugging face api. For different pretrained models one needs to adjust the `model_checkpoint` parameter. To use a weighted loss fontion, one needs to use `trainer = CustomTrainer` instead of `trainer = Trainer`. To use augmentation on needs to set `augmentation = True`.

