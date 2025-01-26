# digID: disease gene IDentification computational framework

This is a pipeline for identifying candidate drug targets for diseases. This project in particular has been implemented to identify candidate drug targets for Alzheimer's disease (AD). Of the top ten ranked candidates, seven are well-studied AD genes, including APP (amyloid precursor protein). The remaining three, GNB1, GNAI1, and KNG1, were barely studied in the context of AD, and therefore could be novel drug targets. This hypothesis has support from experimental evidence produced from experiments conducted at the Brody School of Medicine.

This project is currently in the process of submitting for publication, in collaboration with the Brody School of Medicine.

## Dependencies
Both Python and R are required to run this pipeline. The packages required for both Python and R are specified in `requirements_python.txt` and `requirements_r.txt`, respectively.

## Data
The data for this project is currently stored in [Google Drive](https://drive.google.com/drive/u/1/folders/1AeBow8IXNZIheAO4w3kvgUhp04bW3kYA). This is a folder that contains:
- `GO_PATHWAY_EXP_PPI_combined_data.csv`: gene ontology, pathway annotation, gene expression, and protein-protein interaction data for training the semi-supervised machine learning model to predict AD associated genes
- `AD_labels.csv`: labels for machine learning
- `9606.protein.actions.v10.5.txt`: protein-protein interaction data, used for constructing the interaction network for network analysis to rank the genes predicted by the semi-supervised machine learning model
- `gene_ID_name.csv`: map between different gene naming conventions
- `annotation_string.txt`: STRING protein annotation data (mostly used for switching between naming conventions)

## Predicting AD-associated genes using ladder networks
`predict_AD_associated_genes.py` uses the ladder networks model defined in `ladder_networks.py` in order to predict AD-associated genes after training on `GO_PATHWAY_EXP_PPI_combined_data.csv` and `AD_labels.csv`, and returns the predicted genes into a `results` folder.

The machine learning model is a [Keras implementation by divamgupta](https://github.com/divamgupta/ladder_network_keras) of the ladder networks model proposed in [Rasmus et al., 2015](https://arxiv.org/abs/1507.02672). 

## Prioritize AD-associated genes using network analysis
The machine learning model cannot infer a causal relationship between the features and the labels. Therefore, network analysis is used to prioritize the predicted AD-associated genes that are most central to the disease gene network. These central genes are most likely to be promising drug targets. The `get_centrality_scores.R` script finds these most central genes and returns them into a `results` folder. 
