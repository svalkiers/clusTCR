import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.stats import zscore
from typing import Callable, List, Tuple
from sklearn.base import BaseEstimator, TransformerMixin

properties = {
    "basicity": {
        "A": 206.4,
        "B": 210.7,
        "C": 206.2,
        "D": 208.6,
        "E": 215.6,
        "F": 212.1,
        "G": 202.7,
        "H": 223.7,
        "I": 210.8,
        "K": 221.8,
        "L": 209.6,
        "M": 213.3,
        "N": 212.8,
        "P": 214.4,
        "Q": 214.2,
        "R": 237.0,
        "S": 207.6,
        "T": 211.7,
        "V": 208.7,
        "W": 216.1,
        "X": 210.2,
        "Y": 213.1,
        "Z": 214.9,
    },
    "hydrophobicity": {
        "A": 0.16,
        "B": -3.14,
        "C": 2.50,
        "D": -2.49,
        "E": -1.50,
        "F": 5.00,
        "G": -3.31,
        "H": -4.63,
        "I": 4.41,
        "K": -5.00,
        "L": 4.76,
        "M": 3.23,
        "N": -3.79,
        "P": -4.92,
        "Q": -2.76,
        "R": -2.77,
        "S": -2.85,
        "T": -1.08,
        "V": 3.02,
        "W": 4.88,
        "X": 4.59,
        "Y": 2.00,
        "Z": -2.13,
    },
    "helicity": {
        "A": 1.24,
        "B": 0.92,
        "C": 0.79,
        "D": 0.89,
        "E": 0.85,
        "F": 1.26,
        "G": 1.15,
        "H": 0.97,
        "I": 1.29,
        "K": 0.88,
        "L": 1.28,
        "M": 1.22,
        "N": 0.94,
        "P": 0.57,
        "Q": 0.96,
        "R": 0.95,
        "S": 1.00,
        "T": 1.09,
        "V": 1.27,
        "W": 1.07,
        "X": 1.29,
        "Y": 1.11,
        "Z": 0.91,
    },
    "mutation_stability": {
        "A": 13,
        "C": 52,
        "D": 11,
        "E": 12,
        "F": 32,
        "G": 27,
        "H": 15,
        "I": 10,
        "K": 24,
        "L": 34,
        "M": 6,
        "N": 6,
        "P": 20,
        "Q": 10,
        "R": 17,
        "S": 10,
        "T": 11,
        "V": 17,
        "W": 55,
        "Y": 31,
    },
}


class Resampler(BaseEstimator, TransformerMixin):
    """
    Converts an array of lists with different lengths to an array of
    lists with the same lenght by resampling the lists, effectively
    stretching or squeezing the information but retaining the patterns
    which are present.

    Parameters
    ----------
    length : int, optional
        length to pad to, default = 30
    kind : str, optional
        kind of interpolation algorithm to use in the scipy `interp1d()` function
    parallel: bool, optional
        uses multiple cores to parallelize calculations, default = True
    n_cpus: int, optional
        when using parallel: number of processes to run at the same time, default = 4
    """

    def __init__(
        self,
        length: int = 30,
        kind: str = "linear",
        parallel: bool = True,
        n_cpus: int = 4,
    ) -> None:

        self.length = length
        self.kind = kind
        self.index_out = np.linspace(0, 1, length)
        self.parallel = parallel
        self.n_cpus = n_cpus

    def fit(self, X, y=None):
        """
        Fit transformer.

        Returns
        -------
        self : object
            fitted transformer
        """
        return self

    def _resample(self, array_in: np.array) -> np.array:
        """
        Resamples an input array to match the length specified when initiating the class

        Parameters
        ----------
        array_in: np.array
            Input array

        Returns
        -------
        array_out: np.array
            Output array with the length of ´length_out´
        """
        length_in = len(array_in)
        index_in = np.linspace(0, 1, length_in)
        return interp1d(index_in, array_in, kind=self.kind)(self.index_out)

    def transform(self, X):
        """
        Transform X by resampling the lists of values to equal length.

        Parameters
        ----------
        X : array-like object of lists
            the lists to resample

        Returns:
        --------
        X_out : array of lists
            transformed input, each list is resampled to the required length
        """
        if type(X) == pd.DataFrame:
            conversion = X.applymap(self._resample)
        else:
            conversion = X.apply(self._resample)
        return conversion


# todo : add multiprocessing
class NumericalEncoder(BaseEstimator, TransformerMixin):
    """
    Numerically encodes amino acid sequences based on the physico-chemical
    properties of the amino acids

    Parameters
    ----------
    property : {'basicity', 'hydrophobicity', 'helicity', 'mutation_stability'}
        the amino acid property to encode,
        additional properties can be added via the physicochemical_properties.py file
    """

    def __init__(self, property):

        p = properties.get(property)
        self.property = property
        self.property_dict = self._z_normalize(p)

    @staticmethod
    def _z_normalize(col: pd.Series) -> pd.Series:
        return zscore(pd.Series(col), axis=None, ddof=1).to_dict()

    def fit(self, X, y=None) -> None:
        """
        Fit transformer.

        Returns
        -------
        self : object
            fitted transformer
        """
        return self

    def _encode(self, sequence) -> list:
        """
        Encode a sequence of amino acids to a numeric value according to a
        specific property.
        """
        return [self.property_dict.get(a) for a in sequence]

    def transform(self, X: pd.Series) -> pd.Series:
        """
        Transform X through encoding amino acid properties as values.

        Parameters
        ----------
        X : array-like
            the CDR3 AA sequences to encode.

        Returns:
        --------
        X_out : array of lists
            transformed input, each CDR3 is encoded to a list of values per AA
        """
        return X.apply(self._encode)


class NumericalEncoderAll(BaseEstimator, TransformerMixin):
    """
    Uses the NumericalEncoder class to encode amino acid sequences
    based on all the physico-chemical properties of the amino acids

    Parameters
    ----------
    properties : list, optional
        List of amino acid properties to encode, default = all
    """

    def __init__(self, properties: list = None):

        if properties:
            self.properties = properties
        else:
            self.properties = list(p.keys())

        self.fndict = {}
        for prop in self.properties:
            self.fndict[prop] = NumericalEncoder(property=prop)

    def fit(self, X, y=None):
        """
        Fit transformer.

        Returns
        -------
        self : object
            fitted transformer
        """
        return self

    def transform(self, X):
        """
        Transform X through encoding amino acid properties as values.

        Parameters
        ----------
        X : array-like
            the CDR3 AA sequences to encode.

        Returns:
        --------
        X_out : array of lists
            transformed input, each CDR3 is encoded to a list of values per AA
        """
        res = {}
        for prop in self.properties:
            res[prop] = self.fndict.get(prop).fit_transform(X)
        return pd.DataFrame(res)


class NumericalEncoderAll(BaseEstimator, TransformerMixin):
    """
    Uses the NumericalEncoder class to encode amino acid sequences
    based on all the physico-chemical properties of the amino acids

    Parameters
    ----------
    properties : list, optional
        List of amino acid properties to encode, default = all
    """

    def __init__(self, properties: list = None):

        if properties:
            self.properties = properties
        else:
            self.properties = list(p.keys())

        self.fndict = {}
        for prop in self.properties:
            self.fndict[prop] = NumericalEncoder(property=prop)

    def fit(self, X, y=None):
        """
        Fit transformer.

        Returns
        -------
        self : object
            fitted transformer
        """
        return self

    def transform(self, X):
        """
        Transform X through encoding amino acid properties as values.

        Parameters
        ----------
        X : array-like
            the CDR3 AA sequences to encode.

        Returns:
        --------
        X_out : array of lists
            transformed input, each CDR3 is encoded to a list of values per AA
        """
        res = {}
        for prop in self.properties:
            res[prop] = self.fndict.get(prop).fit_transform(X)
        return pd.DataFrame(res)


class DFFeatureUnion(BaseEstimator, TransformerMixin):
    """
    A `FeatureUnion` class specifically adopted for DataFrames.

    Parameters
    ----------
    transformer_list : [('identifying_string, Function()), ...]
        list of transformations to apply, according to the FeatureUnion conventions
    """

    def __init__(self, transformer_list: List[Tuple[str, Callable]]) -> None:
        self.transformer_list = transformer_list

    def fit(self, X, y=None):
        for (name, t) in self.transformer_list:
            t.fit(X, y)
        return self

    def transform(self, X):

        res = [t.transform(X) for _, t in self.transformer_list]
        out = []

        for r, tup in zip(res, self.transformer_list):
            tname, t = tup
            if isinstance(r, (list, pd.Series)):
                r = pd.Series(r, name=tname)
            elif isinstance(r, np.ndarray):
                if len(r.shape) == 1:
                    r = pd.Series(r, name=tname)
                else:
                    colnames = [f"{tname}_{i}" for i in range(r.shape[1])]
                    r = pd.DataFrame(r, names=colnames)
            elif isinstance(r, pd.DataFrame):
                # colnames = [f"{tname}_{i}" for i in range(r.shape[1])]
                # r.columns = colnames
                pass

            else:
                raise ValueError(f"{type(r)} is not a recognized data type (yet)")

            out.append(r)

        return pd.concat(out, axis=1)


def expand_all(df):
    dfs = {}
    for col in df:
        length = len(df[col][0])
        dfs[col] = pd.DataFrame(
            df[col].to_list(), columns=[f"{col}_{i}" for i in range(length)]
        )
    return pd.concat(dfs.values(), axis=1)
