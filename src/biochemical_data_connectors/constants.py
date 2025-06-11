from enum import Enum


class RestApiEndpoints(Enum):
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