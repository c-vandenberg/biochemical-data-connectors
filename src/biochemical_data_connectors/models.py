from dataclasses import dataclass, field
from typing import Optional, Any


@dataclass
class BioactiveCompound:
    """
    A dataclass to hold standardized information about a bioactive compound.

    Attributes
    ----------
    source_db : str
        The name of the database where the data was sourced, e.g., 'ChEMBL'.
    source_id : str
        The compound's primary identifier from the source database, e.g.,
        a ChEMBL ID or a PubChem CID.
    smiles : str
        The canonical SMILES string representing the molecule.
    target_uniprot : str
        The UniProt accession ID of the biological target.
    activity_type : str
        The type of bioactivity measurement, e.g., 'Kd', 'IC50'.
    activity_value : float
        The numeric value of the bioactivity, typically in nM.
    n_measurements : int
        The total number of measurements found for the primary `activity_type`.
    source_inchikey : str
        The standard 27-character InChIKey for the compound's structure from the source database.
    standardized_inchikey : Optional[str]
        The canonical InChIKey calculated from the standardized molecular
        structure.
    iupac_name : Optional[str]
        The IUPAC name of the compound.
    molecular_formula : Optional[str]
        The molecular formula of the compound, e.g., 'C9H8O4'.
    molecular_weight : Optional[float]
        The molecular weight of the compound's free base.
    raw_data : Optional[Any]
        The original, unprocessed data object from the source library
        (e.g., a `pubchempy.Compound` object) for advanced use cases.
        The string representation of this field is hidden for clarity.
    """
    source_db: str
    source_id: str
    smiles: str
    target_uniprot: str
    activity_type: str
    activity_value: float
    n_measurements: int
    source_inchikey: Optional[str] = None
    standardized_inchikey: Optional[str] = None
    iupac_name: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    mean_activity: Optional[float] = None
    median_activity: Optional[float] = None
    std_dev_activity: Optional[float] = None
    raw_data: Optional[Any] = field(default=None, repr=False)
