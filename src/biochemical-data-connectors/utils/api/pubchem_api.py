import requests
import logging
from typing import List, Iterator, Optional

import pubchempy as pcp

from src.constants import RestApiEndpoints


def get_active_aids(target_gene_id: str) -> List[str]:
    """Query PubChem to get assay IDs for a given target GeneID."""
    assay_id_url = RestApiEndpoints.PUBCHEM_ASSAYS_IDS_FROM_GENE_ID.url(
        target_gene_id=target_gene_id
    )
    response = requests.get(assay_id_url, timeout=10)
    response.raise_for_status()
    data = response.json()

    return data.get("IdentifierList", {}).get("AID", [])


def get_active_aids_wrapper(target_gene_id: str, logger: logging.Logger) -> List[str]:
    try:
        return get_active_aids(target_gene_id)
    except Exception as e:
        if logger:
            logger.error(f"Error retrieving assay IDs for GeneID {target_gene_id}: {e}")
        else:
            print(f"Error retrieving assay IDs for GeneID {target_gene_id}: {e}")

        return []


def get_active_cids(aid: str) -> List[int]:
    """Query PubChem assay details to get active compound IDs for a given assay ID."""
    compound_id_url = RestApiEndpoints.PUBCHEM_COMPOUND_ID_FROM_ASSAY_ID.url(aid=aid)
    response = requests.get(compound_id_url, timeout=10)
    response.raise_for_status()
    data = response.json()

    return data.get("InformationList", {}).get("Information", [])[0].get("CID", [])


def get_active_cids_wrapper(aid: str, logger: logging.Logger = None) -> List[int]:
    try:
        return get_active_cids(aid)
    except Exception as e:
        if logger:
            logger.error(f"Error processing assay {aid}: {e}")
        else:
            print(f"Error processing assay {aid}: {e}")

        return []


def batch_iterable(iterable: List, n: int = 1000) -> Iterator[List]:
    """
    Yield successive n-sized batches from the iterable.

    Parameters
    ----------
    iterable : List
        The list to split into batches.
    n : int, optional
        The batch size, by default 100.

    Yields
    ------
    Iterator[List]
        A list containing a batch of elements from the original iterable.
    """
    iter_len = len(iterable)
    for idx in range(0, iter_len, n):
        yield iterable[idx:idx + n]


def get_compounds_in_batches(
    cids: List[int],
    batch_size: int = 1000,
    logger: logging.Logger = None
) -> List[pcp.Compound]:
    """
    Retrieve compound details from PubChem for a list of compound IDs (CIDs) in batches.

    Parameters
    ----------
    cids : List[int]
        List of PubChem compound IDs.
    batch_size : int, optional
        Number of CIDs per batch to query. Default is 100.
    logger : logging.Logger, optional
        Logger to log error or informational messages.

    Returns
    -------
    List[pcp.Compound]
        A list of compound objects retrieved from PubChem.
    """
    compounds = []
    for cid_batch in batch_iterable(cids, batch_size):
        try:
            batch_compounds = pcp.get_compounds(cid_batch, 'cid')
            compounds.extend(batch_compounds)
        except Exception as e:
            if logger:
                logger.error(f"Error retrieving compounds for batch {cid_batch}: {e}")
            else:
                print(f"Error retrieving compounds for batch {cid_batch}: {e}")

    return compounds


def get_compound_potency(
    compound: pcp.Compound,
    target_gene_id: str,
    bioactivity_measure: str,
    logger: logging.Logger = None
) -> Optional[float]:
    """
    Retrieve a potency value (e.g., Kd in nM) for a compound by querying the
    PubChem bioassay endpoint.

    Parameters
    ----------
    compound : pcp.Compound
        A compound object from PubChem.

    Returns
    -------
    Optional[float]
        The potency value if available, otherwise None.
    """
    cid = compound.cid
    assay_summary_url =  RestApiEndpoints.PUBCHEM_ASSAY_SUMMARY_FROM_CID.url(cid=cid)
    try:
        response = requests.get(assay_summary_url, timeout=10)
        response.raise_for_status()
        response_json = response.json()

        response_table = response_json.get('Table')
        if not response_table:
            return

        response_columns = response_table.get('Columns')
        response_rows = response_table.get('Row')
        if not response_columns or not response_rows:
            return None

        try:
            columns_list = response_columns.get('Column', [])
            target_gene_idx = columns_list.index('Target GeneID')
            activity_name_idx = columns_list.index('Activity Name')
            activity_value_idx = columns_list.index('Activity Value [uM]')
        except ValueError as e:
            logger.error(f'Column not found in bioassay data: {e}')
            return None

        ic50_values = []
        for row in response_rows:
            row_cell = row.get('Cell', [])
            if not row_cell:
                continue

            row_target_gene = row_cell[target_gene_idx]
            row_activity_name = row_cell[activity_name_idx]
            if str(row_target_gene).strip() != str(target_gene_id).strip():
                continue
            if row_activity_name.strip().upper() != bioactivity_measure:
                continue

            # Extract the activity value (in µM) and convert it to nM
            try:
                value_um = float(row_cell[activity_value_idx])
                value_nm = value_um * 1000.0
                ic50_values.append(value_nm)
            except (ValueError, TypeError):
                continue

        if ic50_values:
            return min(ic50_values)

    except Exception as e:
        logger.error(f'Error retrieving potency for CID {cid}: {e}')
        return None

    return None
