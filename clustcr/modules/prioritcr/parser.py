import pandas as pd
from pathlib import Path
from abc import ABC, abstractmethod
import re


class ParserABC(ABC):
    """
    Parser class:
    ------------
    output dataframe: v_segm, j_segm, CDR3
    """

    def __init__(
        self, min_read_count: int = 0, remove_unproductive: bool = False
    ) -> None:
        self.min_read_count = min_read_count
        self.remove_unproductive = remove_unproductive

    @abstractmethod
    def parse(self, filename: str) -> pd.DataFrame:
        pass

    def _remove_unproductive_sequences(self, df):
        """
        removes sequences containing internal termination codons
        """
        return df[df["CDR3"].str.contains("\_|\*", regex=True) == False]


class VdjdbParser(ParserABC):
    def __init__(self, min_quality_score: int = 0) -> None:
        self.pattern = re.compile(r"((TR[AB])([VJ])(\d+)-?(\d?)?)\*?(\d+)?")
        self.min_qual_score = min_quality_score
        super().__init__()

    def parse(self, filename: str) -> pd.DataFrame:
        file = Path(filename)
        df = pd.read_csv(file, sep="\t")
        df.columns = df.columns.str.replace(".", "_", regex=False)
        df = df.query(
            f'gene=="TRB" & species=="HomoSapiens" & vdjdb_score >= {self.min_qual_score}'
        )

        out_df = pd.DataFrame.from_dict(
            {
                "v_segm": df["v_segm"].apply(self._parse_segment),
                "j_segm": df["j_segm"].apply(self._parse_segment),
                "CDR3": df["cdr3"],
                "epitope": df["antigen_epitope"],
            }
        )
        out_df.reset_index(drop=True, inplace=True)
        return out_df

    def _parse_segment(self, segment: str) -> str:
        if not isinstance(segment, str):
            return None
        m = re.search(self.pattern, segment)
        if m:
            return m.group(1)
        else:
            return None


class MixcrParser(ParserABC):
    def parse(self, filename: str) -> pd.DataFrame:
        file = Path(filename)
        usecols = self._determine_usecols(file)
        colnames = [
            "CDR3",
            "CDR3_NT",
            "read_count",
            "read_proportion",
            "v_segm",
            "j_segm",
        ]
        df = pd.read_csv(file, sep="\t", usecols=usecols).rename(
            dict(zip(usecols, colnames)), axis=1
        )[colnames]

        if self.min_read_count:
            df = df.reset_index(drop=True).query(f"read_count >= {self.min_read_count}")

        if self.remove_unproductive:
            df = self._remove_unproductive_sequences(df)

        return df

    def _determine_usecols(self, file):

        with open(file, "r") as f:
            first_line = f.readline()

        if "clonalSequence" in first_line:
            return [
                "aaSeqCDR3",
                "clonalSequence",
                "cloneCount",
                "cloneFraction",
                "allVHitsWithScore",
                "allJHitsWithScore",
            ]

        elif "targetSequences" in first_line:
            return [
                "aaSeqCDR3",
                "targetSequences",
                "cloneCount",
                "cloneFraction",
                "allVHitsWithScore",
                "allJHitsWithScore",
            ]


class OlgaParser(ParserABC):
    def parse(self, filename: str) -> pd.DataFrame:
        file = Path(filename)
        df = pd.read_csv(file, sep="\t", names=["AA", "CDR3", "v_segm", "j_segm"])
        return df[["v_segm", "j_segm", "CDR3"]]


class PogorelyyParser(ParserABC):
    def __init__(
        self, min_read_count: int = 0, remove_unproductive: bool = True
    ) -> None:
        super().__init__()
        self.min_read_count = min_read_count
        self.remove_unproductive = remove_unproductive

    def parse(self, filename: str) -> pd.DataFrame:
        file = Path(filename)
        df = pd.read_csv(
            file,
            sep="\t",
            header=0,
            names=[
                "rank",
                "read_count",
                "read_proportion",
                "CDR3_NT",
                "CDR3",
                "v_segm",
                "j_segm",
            ],
        )
        df = df.reset_index(drop=True).query(f"read_count >= {self.min_read_count}")

        if self.remove_unproductive:
            df = self._remove_unproductive_sequences(df)

        return df


class ImmunoseqParser(ParserABC):
    def parse(self, filename: str):
        usecols = [
            "count (templates/reads)",
            "frequencyCount (%)",
            "nucleotide",
            "aminoAcid",
            "vGeneName",
            "jGeneName",
        ]
        file = Path(filename)
        df = pd.read_csv(file, sep="\t", usecols=usecols).rename(
            columns=dict(
                zip(
                    usecols,
                    [
                        "read_count",
                        "read_proportion",
                        "CDR3_NT",
                        "CDR3",
                        "v_segm",
                        "j_segm",
                    ],
                )
            )
        )

        if self.remove_unproductive:
            df = self._remove_unproductive_sequences(df)
            df = df.dropna(subset=["CDR3"])

        if self.min_read_count:
            df = df.query(f"read_count > {self.min_read_count}")

        return df
