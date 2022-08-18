import numpy as np
import pandas as pd
import multiprocessing
from multiprocessing import Pool
from itertools import chain
from pathlib import Path
from .modules.olga.sequence_generation import (
    SequenceGenerationVJ,
    SequenceGenerationVDJ,
)
from typing import List


class PGen:
    """
    Class used to calcuate generation probabilities, wrapper for OLGA

    Parameters
    ----------
    chain : str {TRA, TRB}
        alpha or beta t-cell receptor chain
    parallel : Bool, optional
        use parallelization when computing multiple Pgens, default True
    n_cpus : int, optional
        number of processes to use in parallel, default = 8
    """

    def __init__(
        self,
        chain: str = "beta",
        parallel: bool = True,
        n_cpus: int = 8,
    ):
        chain_dict = {
            "alpha": "human_T_alpha",
            "TRA": "human_T_alpha",
            "human_T_alpha": "human_T_alpha",
            "beta": "human_T_beta",
            "TRB": "human_T_beta",
            "human_T_beta": "human_T_beta",
        }
        if chain_dict.get(chain) is not None:
            self.chain = chain_dict.get(chain)
        else:
            raise Exception(f"{chain} is not a valid chain")

        self.parallel = parallel
        self.n_cpus = n_cpus
        self.model_folder = Path(f"prioritcr\modules\olga\default_models\human_T_beta")
        self.model = self._construct_model()

    def _construct_model(self):
        """
        Constructs the OLGA model
        Returns the calculator
        chain (alpha or beta): the T-cell receptor chain of interest
        """
        import sys

        sys.path.append("modules/")
        from olga import load_model, generation_probability

        model_folder = self.model_folder
        params_file_name = Path(model_folder, "model_params.txt")
        marginals_file_name = Path(model_folder, "model_marginals.txt")
        V_anchor_pos_file = Path(model_folder, "V_gene_CDR3_anchors.csv")
        J_anchor_pos_file = Path(model_folder, "J_gene_CDR3_anchors.csv")
        alphabet_filename = None

        genomic_data = load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(
            params_file_name, V_anchor_pos_file, J_anchor_pos_file
        )
        self.genomic_data = genomic_data

        if self.chain == "human_T_alpha":
            generative_model = load_model.GenerativeModelVDJ()
            generative_model.load_and_process_igor_model(marginals_file_name)
            self.generative_model = generative_model
            return generation_probability.GenerationProbabilityVJ(
                generative_model, genomic_data, alphabet_filename
            )

        elif self.chain == "human_T_beta":
            generative_model = load_model.GenerativeModelVDJ()
            generative_model.load_and_process_igor_model(marginals_file_name)
            self.generative_model = generative_model
            return generation_probability.GenerationProbabilityVDJ(
                generative_model, genomic_data, alphabet_filename
            )

    def generate_sequences(
        self, n: int, multiprocessing: bool = False, n_cpus: int = 8
    ) -> List[str]:
        """
        Generates a list of CDR3 sequences using the OLGA model

        Parameters
        ----------
        n: int
            number of CDR3 sequences to generate

        Returns
        -------
        sequences : List of CDR3 sequences
        """
        if self.parallel:
            slices = [n // self.n_cpus] * (self.n_cpus - 1) + [
                n // self.n_cpus + n % self.n_cpus
            ]
            with Pool(self.n_cpus) as p:
                out = p.map(self._generate_seq_helper, slices)
                p.close()
                p.join()
            return list(chain.from_iterable(out))

        else:
            generator = SequenceGenerationVDJ(self.generative_model, self.genomic_data)
            out = []
            for i in range(n):
                out.append(generator.gen_rnd_prod_CDR3()[1])
            return out

    def _generate_seq_helper(self, slice):
        generator = SequenceGenerationVDJ(self.generative_model, self.genomic_data)
        return [generator.gen_rnd_prod_CDR3()[1] for i in range(slice)]

    def compute(self, sequence: str, v: str = None, j: str = None) -> float:
        """
        Computes the generation probability for a given CDR3 amino acid sequence

        Parameters
        ----------
        sequence : str
            amino acid sequence of CDR3 loop
        v : str, optional
            V segment formatted according to IMGT convention (e.g. 'TRBV10-3*0')
        j : str, optional
            J segment formatted according to IMGT convention (e.g. 'TRBJ1-5*01')

        Returns
        -------
        generation probability : float

        """
        return self.model.compute_aa_CDR3_pgen(sequence, v, j)

    def compute_regex(self, regex_sequence: str, v: str = None, j: str = None) -> float:
        """
        Compute generation probabilities of regex patterns

        Parameters
        ----------
        regex_sequences : List[str]
            a list of CDR3 sequences
        v : str, optional
            V segment formatted according to IMGT convention (e.g. 'TRBV10-3*0')
        j : str, optional
            J segment formatted according to IMGT convention (e.g. 'TRBJ1-5*01')

        Returns
        generation_probabilities : List[float]
            list of generation probabilities
        """
        return self.model.compute_regex_CDR3_template_pgen(regex_sequence, v, j)

    def compute_hamming_1(self, sequence: str, v: str = None, j: str = None) -> float:
        return self.model.compute_hamming_dist_1_pgen(sequence, v, j)

    def compute_multiple(self, sequence: pd.Series):
        """
        Compute multiple generation probabilites from vectors

        Parameters
        ----------
        sequence : pd.Series[str]
            amino acid sequence of CDR3 loop

        Returns
        -------
        prob : series
            series of generation probabilities
        """

        if self.parallel:
            with multiprocessing.Pool(self.n_cpus) as pool:
                probabilities = pool.map(self.compute, sequence)

        else:
            probabilities = np.vectorize(self.compute)(sequence)

        return pd.Series(probabilities)

    def compute_multiple_vj(
        self, sequence: pd.Series, v: pd.Series = None, j: pd.Series = None
    ):
        """
        Compute multiple generation probabilites from vectors

        Parameters
        ----------
        sequence : pd.Series[str]
            amino acid sequence of CDR3 loop
        v : pd.Series[str], optional
            V segment formatted according to IMGT convention (e.g. 'TRBV10-3')
        j : pd.Series[str], optional
            J segment formatted according to IMGT convention (e.g. 'TRBJ1-5')

        Returns
        -------
        prob : series
            series of generation probabilities
        """
        if self.parallel:
            with multiprocessing.Pool(self.n_cpus) as pool:
                if v is not None:
                    probabilities = pool.starmap(self.compute, zip(sequence, v, j))
                else:
                    probabilities = pool.map(self.compute, sequence)

        else:
            probabilities = sequence.apply(self.compute)

        return pd.Series(probabilities)
