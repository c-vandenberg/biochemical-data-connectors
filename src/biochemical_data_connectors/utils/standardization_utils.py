import logging
from typing import Optional

from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from dimorphite_dl import protonate_smiles

salt_remover = SaltRemover()


class CompoundStandardizer:
    def __init__(self, logger: Optional[logging.Logger]):
        self._logger = logger if logger else logging.getLogger(__name__)
        self._salt_remover = salt_remover

    def standardize_smiles(
        self,
        smiles: str,
        ph: float = 7.4
    ) -> Optional[dict]:
        """
        Performs a full standardization workflow on a single SMILES string.

        The workflow includes:
        1. Desalting the molecule to keep the parent structure.
        2. Standardizing the protonation state to a specific pH using `dimorphite-dl`.
        3. Generating a canonical SMILES string and an InChIKey of the standardized molecule.

        Parameters
        ----------
        smiles : str
            The input SMILES string to be standardized.
        ph : float, optional
            The pH at which to standardize the protonation state. Default is 7.4.

        Returns
        -------
        Optional[dict]
            A dictionary containing the 'smiles', and 'inchi_key' of the
            standardized molecule, or None if the input SMILES is invalid.
        """
        try:
            # 1. Create initial RDKit `mol` object
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self._logger.warning(f"Invalid input SMILES string: '{smiles}'")

            # 2. Remove salts to get the parent molecule
            parent_mol = self._salt_remover.StripMol(mol=mol)
            parent_smiles = Chem.MolToSmiles(parent_mol, canonical=True)

            # 3. Standardize protonation state using the imported function
            #    Ask for only one variant at the specified pH
            protonated_smiles_list: list[str] = protonate_smiles(
                parent_smiles, ph_min=ph, ph_max=ph, max_variants=1
            )

            # If dimorphite returns nothing, use the desalted parent
            final_smiles = protonated_smiles_list[0] if protonated_smiles_list else parent_smiles

            # 4. Create the final, fully standardized RDKit Mol object
            final_mol = Chem.MolFromSmiles(final_smiles)
            if final_mol is None:
                self._logger.warning(f"Failed to create compound after standardization for SMILES string: '{smiles}'")

            # 5. Generate final canonical representations
            canonical_smiles = Chem.MolToSmiles(final_mol, canonical=True, isomericSmiles=True)
            canonical_inchi_key = Chem.MolToInchiKey(final_mol)

            return {
                "smiles": canonical_smiles,
                "inchi_key": canonical_inchi_key
            }

        except Exception as e:
            self._logger.warning(f"Failed to standardize SMILES '{smiles}': {e}")

            return None
