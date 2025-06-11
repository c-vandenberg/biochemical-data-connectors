"""
BiochemicalDataConnectors: A Python package to extract chemical
and biochemical from public databases.
"""
from .connectors.bioactive_compounds_connectors import (
    BaseBioactivesConnector,
    ChEMBLBioactivesConnector,
    PubChemBioactivesConnector
)

from .connectors.ord_connectors import OpenReactionDatabaseConnector

from .utils.api.mappings import uniprot_to_gene_id_mapping, pdb_to_uniprot_id_mapping

__all__ = [
    # --- Base Classes ---
    "BaseBioactivesConnector",

    # --- Concrete Connectors / Extractors ---
    "ChEMBLBioactivesConnector",
    "PubChemBioactivesConnector",
    "OpenReactionDatabaseConnector",

    # --- Public Utility Functions ---
    "uniprot_to_gene_id_mapping",
    "pdb_to_uniprot_id_mapping",
]
