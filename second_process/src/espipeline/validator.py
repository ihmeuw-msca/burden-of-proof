"""
Validation tests
"""
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

from espipeline.filemanager import FileManager


class Validator:
    """
    Validation tests for each risk outcome pair
    """

    name = "unknown"
    data_files = {}

    def __init__(self, filemanager: FileManager):
        self.filemanager = filemanager
        required_files = [
            "risk_cause_metadata.csv",
            "study_data.csv",
            "output_data.csv"
        ]
        for file in required_files:
            if file not in self.data_files:
                raise KeyError(f"data_files must include {file}.")

    def __repr__(self) -> str:
        return f"Validator(name={self.name})"

    def _run_tests(self, pairs: List[str]) -> pd.DataFrame:
        tests = [name.replace("test_", "")
                 for name in dir(self) if name.startswith("test_")]

        results = []
        for pair in pairs:
            result = pd.DataFrame({"pair": pair}, index=[0])
            for test in tests:
                test_fun = getattr(self, f"test_{test}")
                try:
                    test_fun(self.filemanager[pair])
                    result[test] = True
                except AssertionError:
                    result[test] = False
            results.append(result)

        df_results = pd.concat(results)
        df_results.reset_index(inplace=True, drop=True)
        df_results.to_csv(self.filemanager.o_path /
                          f"{self.name}_validation.csv", index=False)
        return df_results

    def run_tests(self, pairs: List[str] = None) -> pd.DataFrame:
        """
        Run validation tests
        """
        if pairs is None:
            if self.name not in self.filemanager.pair_types:
                raise ValueError(f"Does not have pair type {self.name}.")
            pairs = self.filemanager.pairs_by_type[self.name]
        return self._run_tests(pairs)

    def test_files_exist(self, path: Path):
        """
        Test if all required files exist
        """
        files = self.data_files.keys()
        assert all((path / file).exists() for file in files)

    def test_study_data_columns_exist(self, path: Path):
        """
        Test important columns exist in study_data.csv
        """
        df = pd.read_csv(path / "study_data.csv")
        columns = pd.Index(self.data_files["study_data.csv"])
        assert np.all(columns.isin(df.columns))
        assert not np.any(np.isnan(df[columns]))
        assert not np.any(np.isinf(df[columns]))

    @staticmethod
    def test_study_data_is_outlier(path: Path):
        """
        Test is_outlier column has only value 0 or 1
        """
        df = pd.read_csv(path / "study_data.csv")
        assert np.all(df.is_outlier.isin([0, 1]))

    def test_output_data_columns_exist(self, path: Path):
        """
        Test all important columns exist in output_data.csv
        """
        df = pd.read_csv(path / "output_data.csv")
        columns = pd.Index(self.data_files["output_data.csv"])
        assert np.all(columns.isin(df.columns))
        assert not np.any(np.isnan(df[columns]))
        assert not np.any(np.isinf(df[columns]))

    @staticmethod
    def test_output_data_cause(path: Path):
        """
        Test cause lower less or equal to upper in output_data.csv
        and exp of log cause equal to linear cause
        """
        df = pd.read_csv(path / "output_data.csv")
        assert np.all(df.outer_log_cause_lower <= df.inner_log_cause_lower)
        assert np.all(df.inner_log_cause_lower <= df.log_cause)
        assert np.all(df.log_cause <= df.inner_log_cause_upper)
        assert np.all(df.inner_log_cause_upper <= df.outer_log_cause_upper)
        for name in ["outer", "inner"]:
            assert np.allclose(df[[f"{name}_log_cause_lower",
                                   f"{name}_log_cause_upper"]].mean(axis=1),
                               df.log_cause)
        for prefix in ["outer", "inner"]:
            for suffix in ["lower", "upper"]:
                assert np.allclose(np.exp(df[f"{prefix}_log_cause_{suffix}"]),
                                   df[f"{prefix}_linear_cause_{suffix}"])
        np.allclose(np.exp(df.log_cause), df.linear_cause)


class ContinuousValidator(Validator):
    """
    Validator for continuous risk outcome
    """

    name = "continuous"
    data_files = {
        "risk_cause_metadata.csv": [
            "rei_id",
            "cause_id",
            "risk_type",
            "score",
            "pub_bias",
            "risk_lower",
            "risk_upper",
            "risk_unit"
        ],
        "study_data.csv": [
            "rei_id",
            "cause_id",
            "study_id",
            "log_rr",
            "log_rr_se",
            "is_outlier",
            "ref_risk_lower",
            "ref_risk",
            "ref_risk_upper",
            "alt_risk_lower",
            "alt_risk",
            "alt_risk_upper",
            "log_ref_cause",
            "linear_ref_cause",
            "log_alt_cause_lower",
            "log_alt_cause",
            "log_alt_cause_upper",
            "linear_alt_cause_lower",
            "linear_alt_cause",
            "linear_alt_cause_upper"
        ],
        "output_data.csv": [
            "rei_id",
            "cause_id",
            "age_group_id",
            "sex_id",
            "location_id",
            "outer_log_cause_lower",
            "inner_log_cause_lower",
            "log_cause",
            "inner_log_cause_upper",
            "outer_log_cause_upper",
            "risk",
            "outer_linear_cause_lower",
            "inner_linear_cause_lower",
            "linear_cause",
            "inner_linear_cause_upper",
            "outer_linear_cause_upper"
        ],
        "signal_model.pkl": None,
        "linear_model.pkl": None
    }

    @staticmethod
    def test_study_data_risk(path: Path):
        """
        Test risk lower less or equal than upper in study_data.csv
        """
        df = pd.read_csv(path / "study_data.csv")
        assert np.all(df.ref_risk_lower <= df.ref_risk_upper)
        assert np.all(df.alt_risk_lower <= df.alt_risk_upper)
        for name in ["ref", "alt"]:
            assert np.allclose(df[[f"{name}_risk_lower",
                                   f"{name}_risk_upper"]].mean(axis=1),
                               df[f"{name}_risk"])

    @staticmethod
    def test_study_data_cause(path: Path):
        """
        Test cause lower less or equal than upper in study_data.csv
        and the exponential of log cause match linear cause
        """
        df = pd.read_csv(path / "study_data.csv")
        assert np.all(df.log_alt_cause_lower <= df.log_alt_cause_upper)
        assert np.allclose(df[["log_alt_cause_lower",
                               "log_alt_cause_upper"]].mean(axis=1),
                           df.log_alt_cause)
        for name in ["ref_cause",
                     "alt_cause_lower",
                     "alt_cause",
                     "alt_cause_upper"]:
            assert np.allclose(np.exp(df[f"log_{name}"]), df[f"linear_{name}"])


class DichotomousValidator(Validator):
    """
    Validator for dichotomous risk outcome
    """

    name = "dichotomous"
    data_files = {
        "risk_cause_metadata.csv": [
            "rei_id",
            "cause_id",
            "risk_type",
            "score",
            "pub_bias"
        ],
        "study_data.csv": [
            "rei_id",
            "cause_id",
            "study_id",
            "log_rr",
            "log_rr_se",
            "is_outlier"
        ],
        "output_data.csv": [
            "rei_id",
            "cause_id",
            "age_group_id",
            "sex_id",
            "location_id",
            "outer_log_cause_lower",
            "inner_log_cause_lower",
            "log_cause",
            "inner_log_cause_upper",
            "outer_log_cause_upper"
        ],
        "model.pkl": None
    }
