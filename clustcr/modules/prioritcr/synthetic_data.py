from .analysis import Cluster
from .tools import timed
from .generation_probability import PGen
from collections import Counter


class SyntheticRepertoire:
    """
    Synthetic repertoire, generated using OLGA. This
    class can be queried for the presence of certain sequences
    or patterns, or compared against a real repertoire.
    """

    def __init__(self, size: int, multiprocessing: bool = True) -> None:
        self.multiprocessing = multiprocessing
        self.sequences = self._generate_sequences(size)
        self.seq_counter = Counter(self.sequences)

    def __repr__(self) -> str:
        return f"Synthetic repertoire (size:{len(self.sequences)})"

    @timed
    def _generate_sequences(self, size):
        return PGen().generate_sequences(size, multiprocessing=self.multiprocessing)

    def count_sequence(self, sequence: str) -> int:
        return self.seq_counter.get(sequence, 0)

    def count_sequences_from_cluster(self, cluster: Cluster):
        return sum([self.seq_counter.get(s, 0) for s in cluster.sequences])


class MultiSyntheticRepertoire(SyntheticRepertoire):
    """
    Synthetic repertoire consisting of n times a repertoire
    of a specific size.
    """

    def __init__(self, size: int, n: int) -> None:
        self.n = n
        super().__init__(n * size)

    def __repr__(self) -> str:
        return f"Multisynthetic repertoire (size:{self.n} x {len(self.sequences)})"

    def count_sequence(self, sequence: str) -> int:
        return self.seq_counter.get(sequence, 0) / self.n

    def count_sequences_from_cluster(self, cluster: Cluster):
        return sum([self.seq_counter.get(s, 0) for s in cluster.sequences]) / self.n
