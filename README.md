# biochemical-data-connectors

`biochemical-data-connectors` is a Python package for extracting chemical, biochemical, and bioactivity data from public databases like ORD, ChEMBL and PubChem.

## Overview
`biochemical-data-connectors` provides a simple and consistent interface to query major cheminformatics bioinformatics databases for compounds. It is designed to be a modular and reusable tool for researchers and developers in cheminformatics and drug discovery.

### Key Features
1. **Bioactive Compounds**
   * **Unified Interface**: A single, easy-to-use abstract base class for fetching bioactives for a given target.
   * **Multiple Data Sources**: Includes concrete connectors for major public databases:
     1. ChEMBL (`ChEMBLBioactivesExtractor`)
     2. PubChem (`PubChemBioactivesExtractor`)
   * **Powerful Filtering**: Filter compounds by bioactivity type (e.g., Kd, IC50) and potency value.
   * **Efficient Fetching**: Uses concurrency to fetch data from APIs efficiently.
2. **Chemical Reactions**
   * **Local ORD Processing**: Includes a connector (`OpenReactionDatabaseConnector`) to efficiently process a local copy of the Open Reaction Database.
   * **Reaction Role Correction**: Uses RDKit to automatically correct and reassign reactant/product roles from the source data, improving data quality.
   * **Robust SMILES Extraction**: Canonicalizes and validates SMILES strings for both reactants and products to ensure high-quality, standardized output.
   * **Memory-Efficient Processing**: Employs a generator-based extraction method, allowing for iteration over massive reaction datasets with a low memory footprint.

## Installation
You can install this package locally via:
```
pip install biochemical-data-connectors
```

## Quick Start
Here is a simple example of how to retrieve all compounds from ChEMBL with a measured Kd of less than or equal to 1000 nM for the EGFR protein (UniProt ID: `P00533`).
```
import logging
from biochemical_data_connectors import ChEMBLConnector

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# 1. Instantiate the connector for the desired database
chembl_connector = ChEMBLConnector(
    bioactivity_measure='Kd',
    bioactivity_threshold=1000.0, # in nM
    logger=logger
)

# 2. Specify the target's UniProt ID
target_uniprot_id = "P00533" # EGFR

# 3. Get the bioactive compounds
print(f"Fetching bioactive compounds for {target_uniprot_id} from ChEMBL...")
smiles_list = chembl_connector.get_bioactive_compounds(target_uniprot_id)

# 4. Print the results
if smiles_list:
    print(f"\nFound {len(smiles_list)} compounds.")
    print("First 5 compounds:")
    for smiles in smiles_list[:5]:
        print(smiles)
else:
    print("No compounds found matching the criteria.")
```

## Package Structure
```
biochemical-data-connectors/
├── pyproject.toml
├── requirements-dev.txt
├── src/
│   └── biochemical_data_connectors/
│       ├── __init__.py
│       ├── constants.py
│       ├── models.py
│       ├── connectors/
│       │   ├── __init__.py
│       │   ├── ord_connectors.py
│       │   └── bioactive_compounds
│       │       ├── __init__.py
│       │       ├── base_bioactives_connector.py
│       │       ├── chembl_bioactives_connector.py
│       │       └── pubchem_bioactives_connector.py
│       └── utils/
│           ├── __init__.py
│           ├── iter_utils.py
│           ├── standardization_utils.py
│           └── api/
│               ├── __init__.py
│               ├── mappings.py
│               └── pubchem_api.py
├── tests/
│   └── ...
└── README.md
```

## License
This project is licensed under the terms of the [MIT License](https://opensource.org/license/mit).