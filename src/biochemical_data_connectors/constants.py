from enum import Enum

CONVERSION_FACTORS_TO_NM = {
    "NM": 1.0,
    "NANOMOLAR": 1.0,
    "UM": 1000.0,
    "MICROMOLAR": 1000.0,
    "MM": 1_000_000.0,
    "MILLIMOLAR": 1_000_000.0,
    "PM": 0.001,
    "PICOMOLAR": 0.001
}


class RestApiEndpoints(Enum):
    PDB_ID_UNIPROT_MAPPING = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"

    UNIPROT_MAPPING = "https://rest.uniprot.org/idmapping/run"

    UNIPROT_MAPPING_STATUS = "https://rest.uniprot.org/idmapping/status/{job_id}"

    CHEMBL_ACTIVITY = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

    PUBCHEM_ASSAYS_IDS_FROM_GENE_ID = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/"
        "target/geneid/{target_gene_id}/aids/JSON"
    )

    PUBCHEM_COMPOUND_ID_FROM_ASSAY_ID = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/"
        "aid/{aid}/cids/JSON"
    )

    PUBCHEM_ASSAY_SUMMARY_FROM_CID = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        "cid/{cid}/assaysummary/JSON"
    )

    def url(self, **kwargs) -> str:
        """
        Return the fully‐qualified URL, substituting any placeholders
        in the template with the keyword arguments provided.
        """
        return self.value.format(**kwargs)