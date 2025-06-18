import time
import statistics
import requests
import logging
from typing import List, Dict, Any, Optional

import pubchempy as pcp

from biochemical_data_connectors.utils.api.base_api import BaseAPIClient
from biochemical_data_connectors.constants import RestApiEndpoints, CONVERSION_FACTORS_TO_NM
from biochemical_data_connectors.utils.iter_utils import batch_iterable


class PubChemAPIClient(BaseAPIClient):
    """
    A client for interacting with the PubChem API.

    This class encapsulates all direct interactions with the PubChem PUG REST
    API. A persistent `requests.Session` instance with a retry strategy handles
    transient network errors and API rate-limiting.
    """
    def __init__(self, logger: Optional[logging.Logger] = None):
        super().__init__()
        self._logger = logger if logger else logging.getLogger(__name__)

    def get_active_aids(self, target_gene_id: str) -> List[str]:
        """
        Query PubChem's BioAssay database to get all assay IDs (AIDs) associated
        with a specific target, identified by its NCBI GeneID.

        Parameters
        ----------
        target_gene_id : str
            The NCBI GeneID of the target protein.

        Returns
        -------
        List[str]
            A list of assay ID strings, or an empty list if an error occurs.
        """
        # 1) Throttle the request to be polite to the API server.
        time.sleep(0.1)

        # 2) Construct the API URL for finding assay IDs from a GeneID.
        assay_id_url = RestApiEndpoints.PUBCHEM_ASSAYS_IDS_FROM_GENE_ID.url(
            target_gene_id=target_gene_id
        )
        try:
            # 3) Make the API call using the persistent session and parse the JSON to extract the list of AIDs.
            response = self._session.get(assay_id_url, timeout=10)
            response.raise_for_status()
            data = response.json()

            return data.get("IdentifierList", {}).get("AID", [])
        except requests.exceptions.RequestException as e:
            self._logger.error(f"API error retrieving AIDs for GeneID {target_gene_id}: {e}")

            return []

    def get_active_cids(self, aid: str) -> List[int]:
        """
        Query a PubChem assay to get the Compound IDs (CIDs) of all active compounds.

        Parameters
        ----------
        aid : str
            The PubChem Assay ID (AID) to query.

        Returns
        -------
        List[int]
            A list of integer CIDs for active compounds, or an empty list if an error
            occurs.
        """
        # 1) Throttle the request and construct the API URL for finding active CIDs from an Assay ID.
        time.sleep(0.1)
        compound_id_url = RestApiEndpoints.PUBCHEM_COMPOUND_ID_FROM_ASSAY_ID.url(aid=aid)
        try:
            # 2) Make the API call using the persistent session and parse the nested JSON to extract the list of CIDs.
            response = self._session.get(compound_id_url, timeout=10)
            response.raise_for_status()
            data = response.json()
            info_list = data.get("InformationList", {}).get("Information", [])
            if info_list:
                return info_list[0].get("CID", [])
        except requests.exceptions.RequestException as e:
            self._logger.error(f"API error processing assay {aid}: {e}")

            return []

    def get_compound_bioassay_data(
        self,
        compound: pcp.Compound,
        target_gene_id: str,
        bioactivity_measures: List[str]
    ) -> Optional[Dict[str, Any]]:
        """
        Retrieve a potency value (e.g., Kd in nM) for a compound by querying the
        PubChem bioassay endpoint.

        This function queries the PubChem BioAssay summary for a given compound's
        CID. It finds all activity values matching a prioritized list of measures,
        selects the highest-priority measure found, calculates statistics on its
        values, and returns the results.

        If multiple values are found for highest-priority measure the lowest
        (most potent) value is returned.


        Parameters
        ----------
        compound : pcp.Compound
            A `pubchempy.Compound` object for which to retrieve potency.
        target_gene_id : str
            The NCBI GeneID of the target protein, used to filter bioassays.
        bioactivity_measures : str
            A prioritized list of activity types to search for (e.g., ['Kd', 'Ki']).
            The search is case-insensitive.

        Returns
        -------
        Optional[Dict[str, Any]]
            A dictionary of potency statistics, or None if no matching activity
            data is found or an error occurs. The dictionary contains the following
            keys:
                - 'activity_type'
                - 'best_value'
                - 'n_measurements'
                - 'best_value'
                - 'mean_value'
                - 'median_value'
                - 'std_dev_value'

        Notes
        -----
        - All returned activity values and statistics are standardized to the
          nanomolar (nM) scale for consistency.
        - The function robustly handles multiple common concentration units from the
          source data, including picomolar (pM), nanomolar (nM), and micromolar (µM).
        - It first attempts to find an explicit 'Activity Unit' column in the API
          response. If not found, it falls back to assuming a column named
          'Activity Value [uM]' contains micromolar values.
        - Activity data with unrecognized or unsupported units (e.g., '%', 'ppm')
          are logged and gracefully skipped.
        """
        # 1) Throttle the request and build the API URL for the compound's assay summary.
        time.sleep(0.1)
        cid = compound.cid
        assay_summary_url = RestApiEndpoints.PUBCHEM_ASSAY_SUMMARY_FROM_CID.url(cid=cid)
        try:
            # 2) Make the API call and perform initial validation on the JSON response.
            response = self._session.get(assay_summary_url, timeout=10)
            response.raise_for_status()
            response_json = response.json()

            response_table = response_json.get('Table')
            if not response_table:
                return None

            response_columns = response_table.get('Columns').get('Column', [])
            response_rows = response_table.get('Row')
            if not response_columns or not response_rows:
                return None

            # 3) Defensively find the column indices for the data that needs to be extracted.
            try:
                target_gene_idx = response_columns.index('Target GeneID')
                activity_name_idx = response_columns.index('Activity Name')
                activity_value_idx = None
                activity_unit_idx = None
                unit_is_explicit = False

                # Prioritise explicit unit columns, but fallback to assuming µM from the header.
                if 'Activity Unit' in response_columns and 'Activity Value' in response_columns:
                    activity_value_idx = response_columns.index('Activity Value')
                    activity_unit_idx = response_columns.index('Activity Unit')
                    unit_is_explicit = True
                elif 'Activity Value [uM]' in response_columns:
                    activity_value_idx = response_columns.index('Activity Value [uM]')
                    activity_unit_idx = -1
                else:
                    self._logger.error(f"Could not find a valid activity value column in CID {cid} bioassay data")
            except ValueError as e:
                self._logger.error(f'Required column not found in CID {cid} bioassay data: {e}')
                return None

            grouped_activities = {measure.upper(): [] for measure in bioactivity_measures}

            # 4) Iterate through all rows in the bioassay table to find and parse relevant data.
            for row in response_rows:
                row_cell = row.get('Cell', [])

                if not activity_value_idx or not activity_unit_idx:
                    continue

                if not row_cell or len(row_cell) <= max(target_gene_idx, activity_name_idx, activity_value_idx):
                    continue

                row_target_gene = row_cell[target_gene_idx]
                row_activity_name_upper = row_cell[activity_name_idx].strip().upper()

                if not (str(row_target_gene).strip() == str(
                        target_gene_id) and row_activity_name_upper in grouped_activities.keys()):
                    continue

                # 5) For each valid row, check the activity unit and convert the value to nM.
                try:
                    value = float(row_cell[activity_value_idx])
                    unit_str = str(row_cell[activity_unit_idx]).upper() if unit_is_explicit else "UM"

                    conversion_factor = CONVERSION_FACTORS_TO_NM.get(unit_str)
                    if conversion_factor:
                        value_nm = value * conversion_factor
                        grouped_activities[row_activity_name_upper].append(value_nm)
                    else:
                        self._logger.debug(f"Skipping unsupported unit '{unit_str}' for CID {cid}")
                except (ValueError, TypeError):
                    continue

            # 6) Find the highest-priority activity type that has data.
            final_measure_type = None
            final_values = []
            for measure in bioactivity_measures:
                if grouped_activities[measure.upper()]:
                    final_measure_type = measure.upper()
                    final_values = grouped_activities[final_measure_type]
                    break  # Stop at the first (highest-priority) measure found

            if not final_values:
                return None

            # 7) Calculate statistics on the final list of values.
            count = len(final_values)
            return {
                "activity_type": final_measure_type,
                "best_value": min(final_values),
                "n_measurements": count,
                "mean_value": statistics.mean(final_values) if count > 0 else None,
                "median_value": statistics.median(final_values) if count > 0 else None,
                "std_dev_value": statistics.stdev(final_values) if count > 1 else 0.0,
            }

        except requests.exceptions.RequestException as e:
            self._logger.error(f"API error retrieving bioassay data for {cid}: {e}")

            return None


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
    # 1) Iterate through the provided CIDs in batches using the utility function.
    compounds = []
    for cid_batch in batch_iterable(cids, batch_size):
        try:
            # 2) For each batch, call the pubchempy library to retrieve compound data.
            #    pubchempy handles its own connection and retry logic.
            batch_compounds = pcp.get_compounds(cid_batch, 'cid')
            compounds.extend(batch_compounds)
        except Exception as e:
            # 3) Handle any errors during the batch fetch to ensure the process can continue.
            message = f"Error retrieving compounds for batch {cid_batch}: {e}"
            if logger:
                logger.error(message)
            else:
                print(message)

    return compounds
