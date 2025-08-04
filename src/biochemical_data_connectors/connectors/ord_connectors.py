import os
import glob
import warnings
import logging
from typing import Generator, Tuple, List, Callable, Sequence

from ord_schema.message_helpers import load_message
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Contrib.RxnRoleAssignment import identifyReactants

logging.basicConfig(level=logging.INFO)

warnings.filterwarnings("ignore", message="DEPRECATION WARNING: please use MorganGenerator")


class OpenReactionDatabaseConnector:
    """
    A connector for processing the Open Reaction Database (ORD).

    This class provides methods to stream and process reaction data from a
    local copy of the ORD dataset, handling the extraction and standardization
    of reactant and product SMILES strings.

    Attributes
    ----------
    _ord_data_dir : str
        Path to the local directory containing the ORD data files.
    _logger : logging.Logger
        A logger instance for logging messages.
    """
    def __init__(
        self,
        ord_data_dir: str,
        logger: logging.Logger
    ):
        self._ord_data_dir: str = ord_data_dir
        self._logger = logger

    def extract_all_reactions(self)-> Generator[Tuple[str, str], None, None]:
        """
        Generator that yields processed reactant and product SMILES for all
        reactions contained in the open-reaction-database/ord-data repository dataset.

        This method scans a local directory for ORD data files (`.pb.gz`),
        parses them, cleans the reaction SMARTS, reassigns roles, and yields
        a tuple containing a list of reactant SMILES and a list of product SMILES.

        Yields
        ------
        Generator[Tuple[List[str], List[str]], None, None]
            A generator that yields tuples, where the first element is a list
            of reactant SMILES and the second is a list of product SMILES.
        """
        pb_files = glob.glob(os.path.join(self._ord_data_dir, '**', '*.pb.gz'), recursive=True)

        for pb_file in pb_files:
            dataset = load_message(pb_file, dataset_pb2.Dataset)

            for rxn in dataset.reactions:
                rxn_smarts = None
                for identifier in rxn.identifiers:
                    if identifier.type == reaction_pb2.ReactionIdentifier.REACTION_CXSMILES:
                        rxn_smarts = identifier.value
                        break
                if rxn_smarts is None:
                    continue

                try:
                    cleaned_rxn_smiles = identifyReactants.reassignReactionRoles(rxn_smarts)
                except ValueError:
                    continue

                if not cleaned_rxn_smiles:
                    continue

                try:
                    cleaned_rxn = AllChem.ReactionFromSmarts(cleaned_rxn_smiles, useSmiles=True)
                except ValueError:
                    continue

                _, unmodified_reactants, unmodified_products = identifyReactants.identifyReactants(
                    cleaned_rxn
                )

                reactant_smiles = self._get_reactant_smiles_from_cleaned_rxn(cleaned_rxn, unmodified_reactants)
                product_smiles = self._get_product_smiles_from_cleaned_rxn(cleaned_rxn, unmodified_products)

                if len(reactant_smiles) == 0 or len(product_smiles) == 0:
                    continue

                yield reactant_smiles, product_smiles

    def _get_reactant_smiles_from_cleaned_rxn(
        self,
        cleaned_rxn: ChemicalReaction,
        unmodified_reactants=None
    ) -> List[str]:
        """
        Helper to extract and validate reactant SMILES from a cleaned RDKit reaction.

        Parameters
        ----------
        cleaned_rxn : ChemicalReaction
            The RDKit reaction object after role assignment.
        unmodified_reactants : Optional[List[int]], optional
            A list of indices for molecules that act as reagents or catalysts,
            which should be excluded from the final reactant list.

        Returns
        -------
        List[str]
            A list of canonical SMILES strings for the main reactants.
        """
        return self._get_smiles_from_templates(
            get_template_count=cleaned_rxn.GetNumReactantTemplates,
            get_template=cleaned_rxn.GetReactantTemplate,
            unmodified_indices=unmodified_reactants
        )

    def _get_product_smiles_from_cleaned_rxn(
        self,
        cleaned_rxn: ChemicalReaction,
        unmodified_products=None
    ) -> List[str]:
        """
        Helper to extract and validate product SMILES from a cleaned RDKit reaction.

        Parameters
        ----------
        cleaned_rxn : ChemicalReaction
            The RDKit reaction object after role assignment.
        unmodified_products : Optional[List[int]], optional
            A list of indices for product-side molecules that should be excluded.

        Returns
        -------
        List[str]
            A list of canonical SMILES strings for the main products.
        """
        return self._get_smiles_from_templates(
            get_template_count=cleaned_rxn.GetNumProductTemplates,
            get_template=cleaned_rxn.GetProductTemplate,
            unmodified_indices=unmodified_products
        )

    def _get_smiles_from_templates(
        self,
        get_template_count: Callable[[], int],
        get_template: Callable[[int], Chem.Mol],
        unmodified_indices=None
    ) -> List[str]:
        """
        Generic helper to extract and validate SMILES from reaction templates.

        This method removes atom mapping, canonicalizes the SMILES, and validates
        the resulting structure to ensure data quality.

        Parameters
        ----------
        get_template_count : Callable[[], int]
            A function that returns the number of templates (e.g., `rxn.GetNumReactants`).
        get_template : Callable[[int], Chem.Mol]
            A function that returns a molecule template at a given index.
        unmodified_indices : Optional[List[int]], optional
            A list of indices to exclude from the final output.

        Returns
        -------
        List[str]
            A list of validated, canonical SMILES strings.
        """
        num_templates = get_template_count()
        mols = [get_template(i) for i in range(get_template_count())]

        if unmodified_indices is not None:
            main_indices = [i for i in range(num_templates) if i not in unmodified_indices]
            mols = [mols[i] for i in main_indices]

        valid_smiles_list: List = []
        for mol in mols:
            # 1. Remove atom mapping
            self._remove_atom_mapping_from_mol(mol)

            # 2. Convert SMARTS or partial SMILES to SMILES
            raw_smiles = Chem.MolToSmiles(mol, canonical=True)

            # 3. Reparse SMILES to Mol to strictly validate
            parsed = Chem.MolFromSmiles(raw_smiles)
            if parsed is None:
                self._logger.error(f'Invalid SMILES after SMARTS -> SMILES conversion: {raw_smiles}')
                continue

            # 4. If valid, append final canonical SMILES to valid SMILES list
            final_smiles = Chem.MolToSmiles(parsed, canonical=True, isomericSmiles=True)
            valid_smiles_list.append(final_smiles)

        return valid_smiles_list

    def _extract_ord_reaction_smiles(
        self,
        rxn: reaction_pb2.Reaction,
        role_identifier: int
    ) -> List[str]:
        """
        Extracts SMILES from an ORD protobuf message for a specific role.

        This is a legacy helper for parsing older ORD formats. The primary
        extraction method now relies on `REACTION_CXSMILES`.

        Parameters
        ----------
        rxn : reaction_pb2.Reaction
            The reaction protobuf message.
        role_identifier : int
            The role to extract (e.g., `ReactionRole.REACTANT`).

        Returns
        -------
        List[str]
            A list of SMILES strings for the specified role.
        """
        compound_smiles = []

        if role_identifier == reaction_pb2.ReactionRole.REACTANT:
            for rxn_input in rxn.inputs.values():
                for component in rxn_input.components:
                    if component.reaction_role == role_identifier:
                        self._extract_smiles_from_ord_identifiers(component.identifiers, compound_smiles)

            return compound_smiles

        elif role_identifier == reaction_pb2.ReactionRole.PRODUCT:
            for outcome in rxn.outcomes:
                for product in outcome.products:
                    if product.reaction_role == role_identifier:
                        self._extract_smiles_from_ord_identifiers(product.identifiers, compound_smiles)

            return compound_smiles

    @staticmethod
    def _extract_smiles_from_ord_identifiers(
        identifiers: Sequence[reaction_pb2.CompoundIdentifier],
        smiles_list: List
    ):
        """
        Iterates through compound identifiers to find and append valid SMILES.

        Parameters
        ----------
        identifiers : Sequence[reaction_pb2.CompoundIdentifier]
            A sequence of identifier messages from an ORD protobuf.
        smiles_list : List[str]
            The list to which valid canonical SMILES will be appended.

        Returns
        -------
        List[str]
            The updated list of SMILES strings.
        """
        for identifier in identifiers:
            if identifier.type == reaction_pb2.CompoundIdentifier.SMILES:
                mol = Chem.MolFromSmiles(identifier.value)
                if mol:
                    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
                    smiles_list.append(canonical_smiles)

        return smiles_list

    @staticmethod
    def _remove_atom_mapping_from_mol(mol: Chem.Mol):
        """
        Removes atom map numbers from an RDKit molecule object in-place.

        Parameters
        ----------
        mol : Chem.Mol
            The RDKit molecule object to modify.
        """
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
