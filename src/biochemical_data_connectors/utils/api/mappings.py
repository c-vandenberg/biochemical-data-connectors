import time
import requests
from typing import Optional


def pdb_to_uniprot_id_mapping(pdb_id: str) -> Optional[str]:
    """
    Maps a PDB ID to a UniProt accession using the PDBe API.

    Parameters
    ----------
    pdb_id : str
        The PDB ID (e.g., "1A2B").

    Returns
    -------
    Optional[str]
        The first UniProt accession found corresponding to the PDB ID, or None if not found.

    Examples
    --------
    >>> pdb_to_uniprot_id_mapping("1A2B")
    'P12345'
    """
    pdb_id = pdb_id.lower()
    pdb_uniprot_mapping_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
    try:
        mapping_response = requests.get(pdb_uniprot_mapping_url, timeout=10)
        mapping_response.raise_for_status()
        mapping_response_json = mapping_response.json()
        if pdb_id not in mapping_response_json:
            return None

        uniprot_mappings = mapping_response_json[pdb_id].get("UniProt", {})
        if not uniprot_mappings:
            return None

        return next(iter(uniprot_mappings.keys()))
    except Exception as e:
        print(f"Error mapping PDB ID {pdb_id} to UniProt ID: {e}")
        return None


def uniprot_to_gene_id_mapping(uniprot_id: str) -> Optional[str]:
    """
    Map a UniProt accession to an NCBI GeneID using the UniProt ID mapping API.

    Parameters
    ----------
    uniprot_id : str
        The UniProt accession (e.g., "P00533").

    Returns
    -------
    Optional[str]
        The corresponding NCBI GeneID as a string if found, otherwise None.

    Notes
    -----
    This function uses the asynchronous UniProt mapping service.
    """
    uniprot_mapping_url = "https://rest.uniprot.org/idmapping/run"
    uniprot_mapping_params = {
        "from": "UniProtKB_AC-ID",
        "to": "GeneID",
        "ids": uniprot_id
    }
    uniprot_mapping_response = requests.post(uniprot_mapping_url, data=uniprot_mapping_params, timeout=10)
    if uniprot_mapping_response.status_code != 200:
        print(f"Error starting mapping job for {uniprot_id}: {uniprot_mapping_response.text}")
        return None

    job_id = uniprot_mapping_response.json().get("jobId")
    if not job_id:
        print(f"No job ID returned for {uniprot_id}")
        return None

    uniprot_mapping_status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    for _ in range(30):
        status_response = requests.get(uniprot_mapping_status_url, timeout=10)
        status_data = status_response.json()
        if "results" in status_data:
            break
        time.sleep(1)
    else:
        print(f"Mapping job for {uniprot_id} timed out.")
        return None

    result_data = status_data.get('results', [])
    if result_data:
        return result_data[0].get("to", None)

    return None
