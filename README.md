# # Genetic Epilepsy Interaction Networks
Code for constructing a human brain-expression-filtered molecular interaction network from WikiPathways and analysing  their subnetworks from propagation algorithm.


![Graphical abstract](figures/Abstract.png)


## Repository contents

| File | Purpose |
| --- | --- |
| [WP.py](WP.py) | Retrieve human WikiPathways RDF interactions, map gene products to Ensembl gene IDs, apply brain-expression and metabolite filters, and export interaction tables. |
| [metabolites to delete.xlsx](metabolites%20to%20delete.xlsx) | List of chemicals designated for exclusion from the interaction network to avoid overrepresentation. |
| [datavisual.ipynb](datavisual.ipynb) | Jupyter notebook for data analysis and visualization. |

## Requirements

The Python scripts require Python 3 and the following packages:
```bash
python -m pip install pandas requests openpyxl
```

### Run

```bash
python WP.py
```

By default, the script processes all human pathways returned by the WikiPathways endpoint, writes results to `wp_brain_output/`, and stores reusable data in `wp_brain_cache/`.

To use a saved HPA regional-expression table and specify an output directory:

```bash
python WP.py --brain-expression "data/rna_brain_region_hpa.tsv.zip" --brain-threshold 1.0 --output-folder results
```

### Processing and filtering

1. **Retrieve human interactions.** Query WikiPathways RDF for human pathways, interaction participants, source and target relationships, identifiers, and available reference links.
2.**Map gene products.** Use Ensembl cross-references from WikiPathways, with MyGene.info mapping for supported identifiers that lack an Ensembl cross-reference. Multiple gene mappings are expanded and evaluated individually.
3. **Exclude selected metabolites.** Match participants against the supplied chemical identifiers and labels. An interaction is removed if any participant matches the exclusion list.
4. **Apply the expression filter.** Retain gene products with expression of at least **1 nTPM in one or more HPA regions**, unless a different threshold is supplied. Genes below the threshold or without expression data are excluded. Metabolites are exempt from this expression filter.
5. **Export the network.** Remove self-loops and unresolved or unclassified endpoints, preserve available interaction metadata, and generate a separate table of unique undirected node pairs.



## Downstream analysis

The [datavisual.ipynb](datavisual.ipynb) notebook analyzes input genes and interactions through the following steps:

1. Check gene sources and cross-reference the input gene list.
2. Perform GO enrichment and pathway overrepresentation analyses.
3. Examine the distribution of interaction sources.
4. Visualize HotNet2 results.
