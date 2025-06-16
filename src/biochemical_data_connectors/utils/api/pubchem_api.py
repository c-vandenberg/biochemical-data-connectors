import statistics
import requests
import time
import logging
from typing import List, Dict, Any, Optional

import pubchempy as pcp

from biochemical_data_connectors.constants import RestApiEndpoints
from biochemical_data_connectors.utils.iter_utils import batch_iterable


def get_active_aids(target_gene_id: str, logger: Optional[logging.Logger] = None) -> List[str]:
    """
    Query PubChem's BioAssay database to get all assay IDs (AIDs) associated
    with a specific target, identified by its NCBI GeneID.

    Parameters
    ----------
    target_gene_id : str
        The NCBI GeneID of the target protein.
    logger : logging.Logger, optional
        A logger instance for logging potential errors. If None, errors will
        be printed to standard output. Default is None.

    Returns
    -------
    List[str]
        A list of assay ID strings, or an empty list if an error occurs or
        no AIDs are found.
    """
    assay_id_url = RestApiEndpoints.PUBCHEM_ASSAYS_IDS_FROM_GENE_ID.url(
        target_gene_id=target_gene_id
    )
    try:
        response = requests.get(assay_id_url, timeout=10)
        response.raise_for_status()
        data = response.json()

        return data.get("IdentifierList", {}).get("AID", [])
    except Exception as e:
        message = f"Error retrieving assay IDs for GeneID {target_gene_id}: {e}"
        logger.error(message) if logger else print(message)

        return []


def get_active_cids(aid: str, logger: Optional[logging.Logger] = None) -> List[int]:
    """
    Query a PubChem assay by its assay ID to get the CIDs of all active compounds.

    Parameters
    ----------
    aid : str
        The PubChem Assay ID (AID) to query.
    logger : logging.Logger, optional
        A logger instance for logging potential errors. If None, errors will
        be printed to standard output. Default is None.

    Returns
    -------
    List[int]
        A list of integer Compound IDs (CIDs) for active compounds, or an
        empty list if an error occurs.
    """
    compound_id_url = RestApiEndpoints.PUBCHEM_COMPOUND_ID_FROM_ASSAY_ID.url(aid=aid)
    try:
        response = requests.get(compound_id_url, timeout=10)
        response.raise_for_status()
        data = response.json()

        return data.get("InformationList", {}).get("Information", [])[0].get("CID", [])
    except Exception as e:
        message = f"Error processing assay {aid}: {e}"
        logger.error(message) if logger else print(message)

        return []


def get_compounds_in_batches(
    cids: List[int],
    batch_size: int = 1000,
    logger: logging.Logger = None
) -> List[pcp.Compound]:
    """
    Retrieve full compound details from PubChem for a list of CIDs.

    This function processes the CIDs in batches to avoid creating overly
    long API requests and to handle errors gracefully for individual batches.

    Parameters
    ----------
    cids : List[int]
        A list of PubChem Compound IDs (CIDs) to retrieve.
    batch_size : int, optional
        The number of CIDs to include in each batch request to PubChemPy.
        Default is 1000.
    logger : logging.Logger, optional
        A logger instance for logging potential errors during batch processing.
        If None, errors are printed to standard output. Default is None.

    Returns
    -------
    List[pcp.Compound]
        A list of `pubchempy.Compound` objects. This list may be smaller than
        the input list if some CIDs were invalid or if errors occurred.

    Notes
    -----
    Errors encountered during the processing of a specific batch are logged
    and that batch is skipped, allowing the function to continue with the
    remaining batches.
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


def get_compound_bioassay_data(
    compound: pcp.Compound,
    target_gene_id: str,
    bioactivity_measures: List[str],
    logger: logging.Logger = None,
    max_retries: int = 3,
    delay: float = 1.0
) -> Optional[Dict[str, Any]]:
    """
    Retrieve a potency value (e.g., Kd in nM) for a compound by querying the
    PubChem bioassay endpoint.

    This function queries the PubChem BioAssay summary for a given compound's
    CID. It parses the results to find activity data that matches the specified
    `target_gene_id` and `bioactivity_measure`. If multiple values are found,
    the lowest (most potent) value is returned.

    Parameters
    ----------
    compound : pcp.Compound
        A `pubchempy.Compound` object for which to retrieve potency.
    target_gene_id : str
        The NCBI GeneID of the target protein, used to filter bioassays.
    bioactivity_measures : str
        The type of activity to search for (e.g., 'Kd', 'IC50'). The match
        is case-insensitive.
    logger : logging.Logger, optional
        A logger instance for logging potential errors during the process.
        If None, errors will not be logged. Default is None.

    Returns
    -------
    Optional[float]
        The lowest potency value found (in nM), or `None` if no matching
        activity data is found or if an error occurs.

    Notes
    -----
    - The function expects activity values in micromolar (µM) from the PubChem
      API and converts them to nanomolar (nM) before returning.
    """
    cid = compound.cid
    assay_summary_url =  RestApiEndpoints.PUBCHEM_ASSAY_SUMMARY_FROM_CID.url(cid=cid)
    for attempt in range(max_retries):
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

            grouped_activities = {measure.upper(): [] for measure in bioactivity_measures}

            for row in response_rows:
                row_cell = row.get('Cell', [])
                if not row_cell:
                    continue

                row_target_gene = row_cell[target_gene_idx]
                row_activity_name_upper = row_cell[activity_name_idx].strip().upper()

                if (str(row_target_gene).strip() == str(target_gene_id) and row_activity_name_upper in grouped_activities.keys()
                        and row_cell[activity_value_idx]):
                    try:
                        # Extract the activity value (in µM) and convert it to nM
                        value_um = float(row_cell[activity_value_idx])
                        value_nm = value_um * 1000.0

                        # Store the measure name, value, and its priority index
                        grouped_activities[row_activity_name_upper].append(value_nm)
                    except (ValueError, TypeError):
                        continue

            # Find the best list of values based on the priority list
            final_measure_type = None
            final_values = []
            for measure in bioactivity_measures:
                if grouped_activities[measure.upper()]:
                    final_measure_type = measure.upper()
                    final_values = grouped_activities[final_measure_type]
                    break  # Stop at the first (highest-priority) measure found

            if not final_values:
                return None

            # Calculate statistics on the chosen list of values
            count = len(final_values)
            return {
                "activity_type": final_measure_type,
                "best_value": min(final_values),
                "n_measurements": count,
                "mean_value": statistics.mean(final_values) if count > 0 else None,
                "median_value": statistics.median(final_values) if count > 0 else None,
                "std_dev_value": statistics.stdev(final_values) if count > 1 else 0.0,
            }

        except Exception as e:
            message = f"Attempt {attempt + 1}/{max_retries} failed for CID {cid} bioassay data: {e}"
            logger.warning(message) if logger else print(message)
            if attempt < max_retries - 1:
                time.sleep(delay)
            else:
                final_message = f"All {max_retries} retries failed for CID {cid} bioassay data."
                logger.error(final_message) if logger else print(final_message)

    return None
