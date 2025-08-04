from enum import Enum

CONVERSION_FACTORS_TO_NM = {
    'NM': 1.0,
    'NANOMOLAR': 1.0,
    'UM': 1000.0,
    'MICROMOLAR': 1000.0,
    'MM': 1_000_000.0,
    'MILLIMOLAR': 1_000_000.0,
    'PM': 0.001,
    'PICOMOLAR': 0.001
}


class RestApiEndpoints(Enum):
    PDB_UNIPROT_ID_MAPPING = 'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}'

    UNIPROT_PDB_ID_MAPPING = 'https://www.ebi.ac.uk/pdbe/api/mappings/pdb/{uniprot_id}'

    UNIPROT_MAPPING = 'https://rest.uniprot.org/idmapping/run'

    UNIPROT_MAPPING_STATUS = 'https://rest.uniprot.org/idmapping/status/{job_id}'

    CHEMBL_ACTIVITY = 'https://www.ebi.ac.uk/chembl/api/data/activity.json'

    PUBCHEM_ASSAYS_IDS_FROM_GENE_ID = (
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/'
        'target/geneid/{target_gene_id}/aids/JSON'
    )

    PUBCHEM_COMPOUND_ID_FROM_ASSAY_ID = (
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/'
        'aid/{aid}/cids/JSON'
    )

    PUBCHEM_ASSAY_SUMMARY_FROM_CID = (
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'
        'cid/{cid}/assaysummary/JSON'
    )

    BINDINGDB_LIGANDS_FROM_UNIPROT_ID = (
        'https://bindingdb.org/rest/getLigandsByUniprot?uniprot={uniprot_id}{cutoff_str}&response=application/json'
    )

    IUPHAR_TARGET_ID_FROM_UNIPROT = (
        'https://www.guidetopharmacology.org/services/targets?accession={uniprot_id}'
    )

    IUPHAR_INTERACTIONS_FROM_TARGET_ID_FILTERED = (
        'https://www.guidetopharmacology.org/services/targets/{target_id}/interactions?affinityType={p_measure}'
    )

    IUPHAR_LIGAND_STRUCTURE_FROM_LIGAND_ID = (
        'https://www.guidetopharmacology.org/services/ligands/{ligand_id}/structure'
    )

    IUPHAR_LIGAND_MOLECULAR_PROPERTIES_FROM_LIGAND_ID = (
        'https://www.guidetopharmacology.org/services/ligands/{ligand_id}/molecularProperties'
    )

    def url(self, **kwargs) -> str:
        """
        Return the fully‐qualified URL, substituting any placeholders
        in the template with the keyword arguments provided.
        """
        return self.value.format(**kwargs)