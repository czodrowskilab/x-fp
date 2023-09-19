from collections import defaultdict
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from rdkit.Chem import Draw, DataStructs, BondType, AllChem as Chem
from tqdm import tqdm

BOND_DICT = {
    BondType.SINGLE: "-",
    BondType.DOUBLE: "=",
    BondType.TRIPLE: "#",
    BondType.AROMATIC: ":",
}


class FingerprintManager:
    """
    Fingerprint class
    """

    def __init__(self):
        self.subs_per_on_bits = None
        self.spo_filtered = None
        self.frags = None
        self.fp = None
        self.bit_info_name = None
        self.fp_name = None
        self.n_bits = None
        self.radius = None
        self.use_features = None
        self.input_df = None
        self.mol_col = None

    @staticmethod
    def fetch_fp_from_mol(
        in_mol: Chem.Mol, radius: int = 2, n_bits: int = 2048, use_features: bool = False, use_chirality: bool = False
    ) -> tuple[np.ndarray, dict[str, list[int]]]:
        """
        Generate Morgan fingerprint from a RDKit molecule object.

        Parameters
        ----------
        in_mol : Chem.Mol
            RDKit molecule object
        radius : int, optional
            Radius of the circular fingerprint (default is 2).
        n_bits : int, optional
            Number of bits in the fingerprint (default is 2048).
        use_features : bool, optional
            Use Feature Morgan fingerprint or not (default is False).
        use_chirality : bool, optional
            Use chirality or not (default is False).

        Returns
        -------
        tuple[np.ndarray, dict[str, list[int]]]
            Returning a two-element tuple containing the fingerprint as a numpy array as the first element and
            bit info as a dictionary as the second element.

        """
        bit_info = {}
        morgan_fp = Chem.GetMorganFingerprintAsBitVect(
            in_mol,
            radius,
            nBits=n_bits,
            useFeatures=use_features,
            bitInfo=bit_info,
            useChirality=use_chirality,
        )
        morgan_fp_as_array = np.empty(n_bits, dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(morgan_fp, morgan_fp_as_array)
        return morgan_fp_as_array, bit_info

    @staticmethod
    def fetch_fp_from_smiles(
        in_smiles: str,
        radius: int = 2,
        n_bits: int = 2048,
        use_features: bool = False,
        use_chirality: bool = False,
    ) -> tuple[np.ndarray, dict[str, list[int]]]:
        """
        Generate Morgan fingerprint from a SMILES string.

        Parameters
        ----------
        in_smiles : str
            SMILES string
        radius : int, optional
            Radius of the circular fingerprint (default is 2).
        n_bits : int, optional
            Number of bits in the fingerprint (default is 2048).
        use_features : bool, optional
            Use Feature Morgan fingerprint or not (default is False).
        use_chirality : bool, optional
            Use chirality or not (default is False).

        Returns
        -------
        morgan_fp_as_array : tuple[np.ndarray, dict[str, list[int]]]
            Returning a two-element tuple containing the fingerprint as a numpy array as the first element and
            bit info as a dictionary as the second element.

        Raises
        ------
        ValueError
            If the SMILES string is invalid.

        """
        # Catch errors based on invalid SMILES
        try:
            mol = Chem.MolFromSmiles(Chem.CanonSmiles(in_smiles))
            return FingerprintManager.fetch_fp_from_mol(mol, radius, n_bits, use_features, use_chirality)
        except:
            raise ValueError(f"Invalid SMILES present: {in_smiles}")

    def fetch_fp_from_df(
        self,
        input_df: pd.DataFrame,
        smiles_col: str = "SMILES",
        radius: int = 2,
        n_bits: int = 2048,
        use_features: bool = False,
        use_chirality: bool = False,
        mol_col: bool = None,
        canonicalize: bool = True,
    ) -> None:
        """
        Generate either Morgan Fingerprints or Feature Morgan Fingerprints. Bit info is always included.

        Parameters
        ----------
        input_df : pd.DataFrame
            Pandas DataFrame containing SMILES strings.
        smiles_col : str, optional
            Name of the column containing SMILES strings (default is "SMILES").
        radius : int, optional
            Radius of the circular fingerprint (default is 2).
        n_bits : int, optional
            Number of bits in the fingerprint (default is 2048).
        use_features : bool, optional
            Use Feature Morgan fingerprint or not (default is False).
        use_chirality : bool, optional
            Use chirality or not (default is False).
        mol_col : str, optional
            Name of the column containing RDKit molecules (default is None).
        canonicalize : bool, optional
            If True, molecules are canonicalized (default is True).

        Raises
        ------
        ValueError
            If either of the SMILES or Molecules column is not present in the input DataFrame.
        """
        self.input_df = input_df
        self.use_features = use_features
        self.radius = radius
        self.n_bits = n_bits
        self.mol_col = mol_col

        # Check whether the SMILES or mol column is present in the input DataFrame
        if smiles_col not in self.input_df.columns and (
            self.mol_col is None or self.mol_col not in self.input_df.columns
        ):
            raise ValueError("Target column not found in the input DataFrame")

        # Check whether the Molecules column is present in the input DataFrame, if no, then make it
        if self.mol_col is None:
            self.mol_col = "Molecules"

            # generate canonical RDKit molecules in the input DataFrame. These molecules will be used many times.
            self.input_df[self.mol_col] = self.input_df[smiles_col].apply(
                lambda x: Chem.MolFromSmiles(Chem.CanonSmiles(x))
            )
        elif canonicalize:  # molecules are present, will canonicalize them
            self.input_df[self.mol_col] = self.input_df[self.mol_col].apply(
                lambda x: Chem.MolFromSmiles(Chem.MolToSmiles(x))
            )

        # Generate Morgan fingerprints and bit info
        self.fp_name = f"Morgan_FP_{self.radius}_{self.n_bits}"
        self.bit_info_name = f"Bit_Info_{self.radius}_{self.n_bits}"
        self.input_df[self.fp_name], self.input_df[self.bit_info_name] = zip(
            *self.input_df[self.mol_col].apply(
                lambda x: self.fetch_fp_from_mol(x, self.radius, self.n_bits, self.use_features, use_chirality)
            )
        )

        # Saving the fingerprints as an array
        self.fp = np.stack(self.input_df[self.fp_name])

    @staticmethod
    def _gen_frag_for_single_bit(
        subs_per_on_bits: list[tuple[Chem.Mol, int, int]], include_connectivity_invariants: bool = True
    ) -> tuple[list[str], list[tuple[Chem.Mol, int, int]]]:
        """
        Function generates all fragments for a single bit including canonical SMARTS patterns of each fragment.

        Parameters
        ----------
        subs_per_on_bits : list[tuple[Chem.Mol, int, int]]
            A list of tuples containing the molecule, atom index and radius each. The list contains all molecules for
            which one particular bit is set.
        include_connectivity_invariants : bool, optional
            If True, connectivity invariants are included in the fragments (default is True).

        Returns
        -------
        tuple[list[str], list[tuple[Chem.Mol, int, int]]]
            Returning a two-element tuple containing a list of unique canonical SMARTS patterns as first element and a
            list of tuples containing a molecule, atom index and radius each as second element.
        """
        # Check if global "ici" variable is set, if so, use it for the invariant parameter
        if "_ici" in globals():
            include_connectivity_invariants = _ici

        seen = []
        spo_filtered = []
        for mol, atom_index, radius in subs_per_on_bits:
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_index, useHs=True)
            atoms = {atom_index}
            for bidx in env:
                atoms.add(mol.GetBondWithIdx(bidx).GetBeginAtomIdx())
                atoms.add(mol.GetBondWithIdx(bidx).GetEndAtomIdx())

            outer_atoms = set()
            if include_connectivity_invariants:
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius + 1, atom_index, useHs=True)
                atoms_rpp = {atom_index}
                for bidx in env:
                    atoms_rpp.add(mol.GetBondWithIdx(bidx).GetBeginAtomIdx())
                    atoms_rpp.add(mol.GetBondWithIdx(bidx).GetEndAtomIdx())
                outer_atoms = atoms_rpp - atoms
                atoms = atoms_rpp

            atom_symbols = []
            for atom in mol.GetAtoms():
                if atom.GetIdx() in outer_atoms:
                    atom_symbols.append("*")
                    continue
                atom_symbol = Chem.MolFragmentToSmiles(mol, atomsToUse=[atom.GetIdx()], allHsExplicit=True)
                atom_symbol = f"{atom_symbol[:-1]}{';R' if atom.IsInRing() else ''}]"
                atom_symbols.append(atom_symbol)

            bond_symbols = [
                "~"
                if bond.GetBeginAtomIdx() in outer_atoms or bond.GetEndAtomIdx() in outer_atoms
                else BOND_DICT[bond.GetBondType()]
                for bond in mol.GetBonds()
            ]

            smi = Chem.MolFragmentToSmiles(
                mol,
                atomsToUse=list(atoms),
                bondsToUse=env,
                atomSymbols=atom_symbols,
                bondSymbols=bond_symbols,
                allBondsExplicit=True,
                rootedAtAtom=atom_index,
                canonical=True,
            )
            if smi not in seen:
                # Check if substructure encoded by "smi" is already present in "seen" but with different atom order
                matches = set()
                for known_smi in seen:
                    matches.add(tuple(sorted(mol.GetSubstructMatch(Chem.MolFromSmarts(known_smi)))))
                if not tuple(sorted(mol.GetSubstructMatch(Chem.MolFromSmarts(smi)))) in matches:
                    spo_filtered.append((mol, atom_index, radius))
                    seen.append(smi)
        return seen, spo_filtered

    @staticmethod
    def _init_pool(ici: bool) -> None:
        """
        Helper function to be able to use the Pool method imap() instead of starmap() because of better performance. It
        takes a boolean specifying if connectivity invariants should be included or not and sets it to the global
        variable `_ici`.

        Parameters
        ----------
        ici : bool
            Boolean specifying if connectivity invariants should be included or not.
        """
        global _ici
        _ici = ici

    def generate_frags(
        self, verbose: bool = True, n_jobs: int = 1, include_connectivity_invariants: bool = True
    ) -> None:
        """
        Generates fragments for all molecules and all nonzero bits in the dataset saved to this instance. The dataset
        must be saved to this instance using the `fetch_fp_from_df()` method.

        Parameters
        ----------
        verbose : bool, optional
            If True, status messages are printed (default is True).
        n_jobs : int, optional
            Number of cores to use for parallelization (default is 1).
        include_connectivity_invariants : bool, optional
            If True, connectivity invariants are included in the fragments (default is True).

        Raises
        ------
        ValueError
            If no input data is present, i.e. `fetch_fp_from_df()` was not called already.
        """
        # Check if there is input data
        if self.input_df is None:
            raise ValueError("Input data is not present. Please use fetch_fp_from_df() method to add input data.")

        # Check n_jobs parameter
        if n_jobs < 1:
            n_jobs = 1
        elif n_jobs > cpu_count():
            n_jobs = cpu_count()

        if verbose:
            print("Generating fragments. This step may take a while.")

        # Gather information for all the bit sets over the entire dataset
        self.subs_per_on_bits = defaultdict(list)

        for i, (_, row) in enumerate(self.input_df.iterrows()):
            single_mol = row[self.mol_col]
            single_mol_bit_info = row[f"Bit_Info_{self.radius}_{self.n_bits}"]
            for bit in single_mol_bit_info:
                for atom_index, radius in single_mol_bit_info[bit]:
                    self.subs_per_on_bits[bit].append((single_mol, atom_index, radius))

        # Print the bits which are set
        if verbose:
            print(f"{len(self.subs_per_on_bits)} / {self.n_bits} bits are set in the dataset.")

        # Filtering duplicates per bit and extract fragment SMARTS
        if n_jobs == 1:
            self.spo_filtered = {}
            self.frags = {}
            for bit in tqdm(self.subs_per_on_bits):
                self.frags[bit], self.spo_filtered[bit] = self._gen_frag_for_single_bit(
                    self.subs_per_on_bits[bit], include_connectivity_invariants
                )
        else:
            if verbose:
                print(f"Parallelize using {n_jobs} cores.")

            with Pool(n_jobs, initializer=self._init_pool, initargs=(include_connectivity_invariants,)) as p:
                frags, spo_filtered = zip(*p.imap(self._gen_frag_for_single_bit, self.subs_per_on_bits.values()))

            self.frags = dict(zip(self.subs_per_on_bits.keys(), frags))
            self.spo_filtered = dict(zip(self.subs_per_on_bits.keys(), spo_filtered))

    def make_substruct_per_bit_plot(self, dpi: int = 100, save_fig: bool = False, fig_name: str = None) -> None:
        """
        Generate substructures per bit plot. This plot is useful to understand the distribution of substructures encoded
        in each bit.

        Parameters
        ----------
        dpi : int, optional
            Dots per inch (default is 100).
        save_fig : bool, optional
            If True, the figure is saved (default is False).
        fig_name : str, optional
            Name of the plot (default is None).
        """
        # Check for frags and spo_filtered
        if self.frags is None or self.spo_filtered is None:
            print("WARNING: Fragments were not generated, yet. Try to generate now with default settings...")
            self.generate_frags()

        plt.figure(figsize=(2, 4), dpi=dpi)
        fig = sns.boxplot([len(self.frags[bit]) for bit in self.frags])
        fig.set_ylabel("Substructures per bit")
        fig.set_xlabel(f'{"F" if self.use_features else "E"}CFP{self.radius * 2} {self.n_bits} bits')

        if save_fig:
            if fig_name is None:
                fig_name = f"Substructures_per_bit_plot.png"
            plt.savefig(fig_name, dpi=dpi)
        else:
            plt.show()

    def substruct_counter_per_bit(self, query_bit: int, verbose: bool = False) -> pd.DataFrame:
        """
        Counts the frequency of the substructures present in a query bit. Total occurrences will calculate all the
        occurrences of a substructure in the dataset, including multiple occurrences in a molecule. Unique occurrences
        will calculate the number of molecules in which a substructure is present.

        Parameters
        ----------
        query_bit : int
            Bit for which the substructure frequency is to be calculated.
        verbose : bool, optional
            If True, status messages are printed (default is False).

        Returns
        -------
        pd.DataFrame
            Pandas DataFrame where substructures (as SMARTS) are present with their total and unique occurrences.
        """
        # Check for frags and spo_filtered
        if self.frags is None or self.spo_filtered is None:
            print("WARNING: Fragments were not generated, yet. Try to generate now with default settings...")
            self.generate_frags()

        # Initialize an empty list to store unique occurrences.
        substruct_unique_occurrences = []

        # Initialize an empty list to store total occurrences.
        substruct_total_occurrences = []

        # Print wait-time warning
        if verbose:
            print(f"Counting substructure frequencies of the query bit: Bit {query_bit}. This step may take a while.")

        # Run a loop of all the substructures over the data to prepare their count
        for substructure in self.frags[query_bit]:
            frequency_substructure_per_mol = 0
            frequency_substructure_per_mol_unique = 0
            substructure_as_mol = Chem.MolFromSmarts(substructure)

            # Loop over molecules in the dataset
            for single_mol in self.input_df[self.mol_col]:
                matches = single_mol.GetSubstructMatches(substructure_as_mol)
                frequency_substructure_per_mol += len(matches)

                # To calculate unique occurrences
                if len(matches) > 0:
                    frequency_substructure_per_mol_unique += 1

            substruct_total_occurrences.append(frequency_substructure_per_mol)
            substruct_unique_occurrences.append(frequency_substructure_per_mol_unique)

        # Creating the substructure dataframe
        substruct_df = pd.DataFrame(
            {
                "Serial": range(1, len(self.frags[query_bit]) + 1),
                "Substructures (as SMARTS)": self.frags[query_bit],
                "Total occurrences": substruct_total_occurrences,
                "Unique occurrences": substruct_unique_occurrences,
            }
        )
        return substruct_df

    def get_substruct_images_per_bit(self, query_bit: int, smarts_as_legends: bool = False) -> str:
        """
        Generates an SVG image of substructures present in a query bit.

        Parameters
        ----------
        query_bit : int
            Bit for which the substructure images are to be generated.
        smarts_as_legends : bool, optional
            If True, substructure SMARTS are used as legends (default is False).

        Returns
        -------
        str
            SVG image as a string of substructures present in a query bit. This image should be saved separately.

        """
        # Check for frags and spo_filtered
        if self.frags is None or self.spo_filtered is None:
            print("WARNING: Fragments were not generated, yet. Try to generate now with default settings...")
            self.generate_frags()

        # Fix Draw Option issue:
        draw_options = Draw.rdMolDraw2D.MolDrawOptions()
        draw_options.prepareMolsBeforeDrawing = False

        # legends
        if smarts_as_legends:
            legends = self.frags[query_bit]
        else:
            legends = [str(i + 1) for i in range(len(self.frags[query_bit]))]

        # Get substructure images
        substruct_images = Draw.DrawMorganEnvs(
            self.spo_filtered[query_bit],
            legends=legends,
            subImgSize=(300, 300),
            molsPerRow=4,
            drawOptions=draw_options,
        )
        return substruct_images
