"""
Main GBD Evience Score Pipeline
"""
import shutil
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple, Union, ValuesView

import dill
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mrtool import MRBRT, LinearCovModel, MRBeRT, MRData
from numpy import ndarray
from pandas import DataFrame

from espipeline.process import Process
from espipeline.utils import get_beta_info, get_gamma_info, get_pub_bias


class PublicationNoTrimmingProcess(Process):
    """
    Sensitivity analysis for trimming
    """
    def __init__(self, i_path: Path, o_path: Path):
        super().__init__(i_path, o_path)
        self.o_path = self.o_path.parents[0]
        df_risk = pd.read_csv("ids_rei_2021-05-12.csv")
        df_cause = pd.read_csv("ids_cause_2021-05-12.csv")
        
        rei_name = dict(zip(df_risk.rei_id, df_risk.rei_name))
        cause_name = dict(zip(df_cause.cause_id, df_cause.cause_name))
        acause = dict(zip(df_cause.cause_id, df_cause.acause))
        
        self.rei_name = rei_name[self.risk_cause_metadata.rei_id.values[0]]
        self.rei = self.rei_name.replace("High ", "").replace("Low ", "").replace(" low", "").replace(" high", "").lower()
        
        self.acause = acause[self.risk_cause_metadata.cause_id.values[0]]
        self.cause_name = cause_name[self.risk_cause_metadata.cause_id.values[0]]
        
        self.new_pair_model = PairModel(
            self.name,
            self.signal_model,
            self.new_linear_model,
            exp_limits=(self.output_data.risk.min(),
                        self.output_data.risk.max())
        )
        
    def run(self):
        self.get_no_trimming_model()
        self.plot_models()
   
    def get_no_trimming_model(self):
        self.nt_signal_model = deepcopy(self.signal_model)
        self.nt_signal_model.inlier_pct = 1.0
        for model in self.nt_signal_model.sub_models:
            model.inlier_pct = 1.0
        self.nt_signal_model.fit_model(inner_max_iter=1000, inner_print_level=5)
        
        self.nt_pair_model = PairModel(
            self.name,
            self.nt_signal_model,
            self.linear_model,
            exp_limits=(self.output_data.risk.min(),
                        self.output_data.risk.max())
        )   
        
        # get data for new linear model
        df_signal = self.nt_signal_model.data.to_df()
        df_signal["signal"] = self.nt_pair_model.get_data_signal()
        df_signal["re_signal"] = self.nt_pair_model.get_data_re_signal(use_new=True)
        df_signal.study_id = df_signal.study_id.astype(int)
        
        df_new_linear = df_signal.copy()
        for name in self.new_pair_model.bias_covs:
            df_new_linear[name] = 0.0

        data = MRData()
        data.load_df(
            df_new_linear,
            col_obs="obs",
            col_obs_se="obs_se",
            col_study_id="study_id",
            col_covs=["signal", "re_signal"] + self.new_pair_model.bias_covs
        )

        cov_models = [
            LinearCovModel("signal",
                           use_re=False,
                           prior_beta_uniform=[1.0, 1.0]),
            LinearCovModel("re_signal",
                           use_re=True,
                           prior_beta_uniform=[0.0, 0.0]),
            LinearCovModel("intercept",
                           use_re=True,
                           prior_beta_uniform=[0.0, 0.0])
        ]
        for bias_cov in self.new_pair_model.bias_covs:
            old_cov_model = self.linear_model.get_cov_model(bias_cov)
            cov_models.append(
                LinearCovModel(bias_cov,
                               use_re=False,
                               prior_beta_gaussian=old_cov_model.prior_beta_gaussian.copy())
            )
        self.nt_linear_model = MRBRT(data, cov_models=cov_models)
        self.nt_linear_model.fit_model(inner_print_level=5, inner_max_iter=1000)
        
        self.nt_pair_model = PairModel(
            self.name,
            self.nt_signal_model,
            self.nt_linear_model,
            exp_limits=(self.output_data.risk.min(),
                        self.output_data.risk.max())
        )
        
        # get new study data
        self.nt_study_data = self.study_data.copy()
        pred = self.nt_pair_model.get_pred(self.nt_study_data.ref_risk.values)
        self.nt_study_data["log_ref_cause"] = pred.log_cause.values
        self.nt_study_data["log_alt_cause"] = (
            self.nt_study_data.log_rr.values +
            self.nt_study_data.log_ref_cause.values
        )
        self.nt_study_data["log_alt_cause_lower"] = (
            self.nt_study_data.log_alt_cause.values -
            1.96*self.nt_study_data.log_rr_se.values
        )
        self.nt_study_data["log_alt_cause_lower"] = (
            self.nt_study_data.log_alt_cause.values +
            1.96*self.nt_study_data.log_rr_se.values
        )

        log_col_names = ["log_ref_cause",
                         "log_alt_cause_lower",
                         "log_alt_cause",
                         "log_alt_cause_upper"]
        for log_col_name in log_col_names:
            linear_col_name = log_col_name.replace("log", "linear")
            self.nt_study_data[linear_col_name] = np.exp(self.nt_study_data[log_col_name])

        # get new output data
        self.nt_output_data = self.output_data.copy()
        pred = self.nt_pair_model.get_pred(self.output_data.risk.values)
        log_col_names = ["outer_log_cause_lower",
                         "inner_log_cause_lower",
                         "log_cause",
                         "inner_log_cause_upper",
                         "outer_log_cause_upper"]
        for log_col_name in log_col_names:
            linear_col_name = log_col_name.replace("log", "linear")
            self.nt_output_data[log_col_name] = pred[log_col_name]
            self.nt_output_data[linear_col_name] = np.exp(pred[log_col_name])

        # get new risk cause metadata
        self.nt_risk_cause_metadata = self.risk_cause_metadata.copy()
        self.nt_risk_cause_metadata.score = self.nt_pair_model.get_score(pred)
        self.nt_risk_cause_metadata.pub_bias = self.nt_pair_model.get_pub_bias()
        self.nt_risk_cause_metadata.risk_lower = self.nt_pair_model.exp_bounds[0]
        self.nt_risk_cause_metadata.risk_upper = self.nt_pair_model.exp_bounds[1]
        self.nt_risk_cause_metadata.to_csv(self.i_path / "nt_risk_cause_metadata.csv", index=False)
        
    def plot_models(self):
        # plot the comparison
        _, ax = plt.subplots(1, 2, figsize=(16, 3.5), sharey=True)
        # plot data
        ax[0].scatter(
            self.new_study_data.alt_risk,
            self.new_study_data.log_alt_cause,
            s=5.0/self.new_study_data.log_rr_se,
            color="gray",
            alpha=0.5
        )
        outlier_index = self.new_study_data.is_outlier == 1
        ax[0].scatter(
            self.new_study_data.alt_risk[outlier_index],
            self.new_study_data.log_alt_cause[outlier_index],
            s=5.0/self.new_study_data.log_rr_se[outlier_index],
            color="red",
            alpha=0.5,
            marker="x"
        )
        ax[1].scatter(
            self.nt_study_data.alt_risk,
            self.nt_study_data.log_alt_cause,
            s=5.0/self.nt_study_data.log_rr_se,
            color="gray",
            alpha=0.5
        )
        
        # plot curves
#         cov_model = self.signal_model.sub_models[0].cov_models[0]
#         risk = self.new_output_data.risk.values
#         covs = {name: np.repeat(risk.min(), risk.size) for name in cov_model.ref_cov}
#         covs.update({name: risk for name in cov_model.alt_cov})
#         signal = self.signal_model.predict(MRData(covs=covs))
        
#         wt_beta_info = get_beta_info(self.new_linear_model)
#         nt_beta_info = get_beta_info(self.nt_linear_model)
        
#         wt_curve = signal*wt_beta_info[0]
#         nt_curve = signal*nt_beta_info[0]
        
        colors = ["#008080", "#FFA500"]
        labels = ["with trimming", "no trimming"]
        
        for i, output_data in enumerate([self.new_output_data, self.nt_output_data]):
            ax[i].plot(output_data.risk, output_data.log_cause, color=colors[i], label=labels[i])
            ax[i].fill_between(output_data.risk,
                               output_data.inner_log_cause_lower,
                               output_data.inner_log_cause_upper,
                               color=colors[i],
                               alpha=0.2)
            ax[i].fill_between(output_data.risk,
                               output_data.outer_log_cause_lower,
                               output_data.outer_log_cause_upper,
                               color=colors[i],
                               alpha=0.2)
            log_cause_bound = output_data.outer_log_cause_lower if output_data.log_cause.mean() >= 0 else output_data.outer_log_cause_upper
            ax[i].plot(output_data.risk,
                       log_cause_bound,
                       color="red",
                       linewidth=1.2, alpha=0.9)
            ax[i].plot(output_data.risk,
                       log_cause_bound,
                       color="gray",
                       linewidth=1.2, alpha=0.5)
        
        # plot risk score bounds
        for b in [self.risk_cause_metadata.risk_lower.values[0],
                  self.risk_cause_metadata.risk_upper.values[0]]:
            for i in range(2):
                ax[i].axvline(b, linestyle="--", linewidth=1, color="gray")
        
        # plot 0 line
        title = (f"{self.rei_name} / {self.cause_name}")
        for i in range(2):
            ax[i].axhline(0.0, linestyle="-", linewidth=1, color="gray")
            ax[i].set_xlabel(f"{self.rei} ({self.risk_cause_metadata.risk_unit.values[0]})")
            ax[i].set_ylabel("log relative risk")
            ax[i].set_title(title, loc="left")
            ax[i].legend()
        
        plt.savefig(self.o_path / f"{self.name}_trimming.pdf", bbox_inches="tight")
        plt.close("all")


class PublicationFEProcess(Process):
    """
    Fixed effect model
    """
    def __init__(self, i_path: Path, o_path: Path):
        super().__init__(i_path, o_path)
        self.o_path = self.o_path.parents[0]
        df_risk = pd.read_csv("ids_rei_2021-05-12.csv")
        df_cause = pd.read_csv("ids_cause_2021-05-12.csv")
        
        rei_name = dict(zip(df_risk.rei_id, df_risk.rei_name))
        cause_name = dict(zip(df_cause.cause_id, df_cause.cause_name))
        acause = dict(zip(df_cause.cause_id, df_cause.acause))
        
        self.rei_name = rei_name[self.risk_cause_metadata.rei_id.values[0]]
        self.rei = self.rei_name.replace("High ", "").replace("Low ", "").replace(" low", "").replace(" high", "").lower()
        
        self.acause = acause[self.risk_cause_metadata.cause_id.values[0]]
        self.cause_name = cause_name[self.risk_cause_metadata.cause_id.values[0]]
   
    def run(self):
        # fit linear model without random effects
        data = self.new_linear_model.data

        cov_models = [LinearCovModel("signal", use_re=False, prior_beta_uniform=[0.0, 5.0])]
        for bias_cov in self.new_linear_model.cov_names:
            if bias_cov in ["signal", "re_signal"]:
                continue
            old_cov_model = self.new_linear_model.get_cov_model(bias_cov)
            cov_models.append(
                LinearCovModel(bias_cov,
                               use_re=False,
                               prior_beta_gaussian=old_cov_model.prior_beta_gaussian.copy())
            )
        linear_model = MRBRT(data, cov_models=cov_models)
        linear_model.fit_model(inner_max_iter=500)
        
        # plot the comparison
        _, ax = plt.subplots(figsize=(8, 5))
        # plot data
        ax.scatter(
            self.new_study_data.alt_risk,
            self.new_study_data.log_alt_cause,
            s=5.0/self.new_study_data.log_rr_se,
            color="gray",
            alpha=0.5
        )
        outlier_index = self.new_study_data.is_outlier == 1
        ax.scatter(
            self.new_study_data.alt_risk[outlier_index],
            self.new_study_data.log_alt_cause[outlier_index],
            s=5.0/self.new_study_data.log_rr_se[outlier_index],
            color="red",
            alpha=0.5,
            marker="x"
        )
        
        # plot curves
        cov_model = self.signal_model.sub_models[0].cov_models[0]
        risk = self.new_output_data.risk.values
        covs = {name: np.repeat(risk.min(), risk.size) for name in cov_model.ref_cov}
        covs.update({name: risk for name in cov_model.alt_cov})
        signal = self.signal_model.predict(MRData(covs=covs))
        
        re_beta_info = get_beta_info(self.new_linear_model)
        fe_beta_info = get_beta_info(linear_model)
        
        re_curve = signal*re_beta_info[0]
        fe_curve = signal*fe_beta_info[0]
        
        ax.plot(risk, re_curve, color="#008080", label="mixed effects")
        ax.plot(risk, fe_curve, color="#FFA500", label="fixed effects")
        
        # plot risk score bounds
        for b in [self.risk_cause_metadata.risk_lower.values[0],
                  self.risk_cause_metadata.risk_upper.values[0]]:
            ax.axvline(b, linestyle="--", linewidth=1, color="gray")
        
        # plot 0 line
        ax.axhline(0.0, linestyle="-", linewidth=1, color="gray")

        # add unit to the xaxis
        ax.set_xlabel(f"{self.rei} ({self.risk_cause_metadata.risk_unit.values[0]})")
        ax.set_ylabel("log relative risk")

        # title
        title = (f"{self.rei_name} / {self.cause_name}")
        ax.set_title(title, loc="left")
        
        # legend
        ax.legend()
        
        plt.savefig(self.o_path / f"{self.name}_compare.pdf", bbox_inches="tight")
        plt.close("all")


class PublicationSimProcess(Process):
    """
    Simulation for the publication, show the stochastic error not induce heterogeniety
    """
    num_sim = 100
    
    def __init__(self, i_path: Path, o_path: Path):
        super().__init__(i_path, o_path)
        self.o_path = self.o_path.parents[0]
    
    def run(self) -> pd.DataFrame:
        data = self.new_linear_model.data
        cov_models = self.new_linear_model.cov_models
        gamma = []
        gamma_sd = []
        for i in range(self.num_sim):
            # simulate data
            data = deepcopy(self.new_linear_model.data)
            data.obs = data.covs["signal"] + np.random.randn(data.num_obs)*data.obs_se
            # copy cov_models
            cov_models = deepcopy(self.new_linear_model.cov_models)
            # create model
            linear_model = MRBRT(data, cov_models)

            linear_model.fit_model(inner_max_iter=1000)
            gamma_info = get_gamma_info(linear_model)
            gamma.append(gamma_info[0])
            gamma_sd.append(gamma_info[1])
        return pd.DataFrame({
            "pair": self.name,
            "gamma": gamma,
            "gamma_sd": gamma_sd
        })


class PublicationScoreProcess(Process):
    """
    Produce different variation of the scores
    """
    def __init__(self, i_path: Path, o_path: Path):
        super().__init__(i_path, o_path)
        self.o_path = self.o_path.parents[0]
        self.new_pair_model = PairModel(
            self.name,
            self.signal_model,
            self.new_linear_model,
            exp_limits=(self.output_data.risk.min(),
                        self.output_data.risk.max())
        )
        
    def run(self) -> pd.DataFrame:
        pred = self.new_pair_model.get_pred(self.output_data.risk.values)
        score_15_85 = self.new_pair_model.get_score(pred)
        # compute score for 05 95
        self.new_pair_model.exp_bounds = (
            np.quantile(self.new_pair_model.ref_mid, 0.05),
            np.quantile(self.new_pair_model.alt_mid, 0.95)
        )
        score_05_95 = self.new_pair_model.get_score(pred)
        # compute score for 10 90
        self.new_pair_model.exp_bounds = (
            np.quantile(self.new_pair_model.ref_mid, 0.10),
            np.quantile(self.new_pair_model.alt_mid, 0.90)
        )
        score_10_90 = self.new_pair_model.get_score(pred)
        return pd.DataFrame({
            "pair": self.name,
            "score_15_85": score_15_85,
            "score_05_95": score_05_95,
            "score_10_90": score_10_90
        }, index=[0])
        

class PulicationFigureProcess(Process):
    """
    Produce figures for publication
    """
    
    def __init__(self, i_path: Path, o_path: Path):
        super().__init__(i_path, o_path)
        self.o_path = self.o_path.parents[0]
        self.new_pair_model = PairModel(
            self.name,
            self.signal_model,
            self.new_linear_model,
            exp_limits=(self.output_data.risk.min(),
                        self.output_data.risk.max())
        )
        df_risk = pd.read_csv("ids_rei_2021-05-12.csv")
        df_cause = pd.read_csv("ids_cause_2021-05-12.csv")
        
        rei_name = dict(zip(df_risk.rei_id, df_risk.rei_name))
        cause_name = dict(zip(df_cause.cause_id, df_cause.cause_name))
        acause = dict(zip(df_cause.cause_id, df_cause.acause))
        
        self.rei_name = rei_name[self.risk_cause_metadata.rei_id.values[0]]
        self.rei = self.rei_name.replace("High ", "").replace("Low ", "").replace(" low", "").replace(" high", "").lower()
        
        self.acause = acause[self.risk_cause_metadata.cause_id.values[0]]
        self.cause_name = cause_name[self.risk_cause_metadata.cause_id.values[0]]

    def run(self):
        # plot figures
        plot_model_publication(self.new_study_data,
                               self.new_output_data,
                               self.new_risk_cause_metadata,
                               self.new_pair_model,
                               self.rei,
                               self.rei_name,
                               self.acause,
                               self.cause_name)
        plt.savefig(self.o_path / f"{self.name}.pdf", bbox_inches="tight")
        plt.close("all")
        

class ScoreChangeProcess(Process):
    """
    Change the way of how to evaluating the score.
    """

    def __init__(self, i_path: Path, o_path: Path):
        super().__init__(i_path, o_path)

        self.new_pair_model = PairModel(
            self.name,
            self.signal_model,
            self.new_linear_model,
            exp_limits=(self.output_data.risk.min(),
                        self.output_data.risk.max())
        )

    def run(self):
        # change results
        pred = self.new_pair_model.get_pred(self.output_data.risk.values)
        self.new_risk_cause_metadata.score = self.new_pair_model.get_score(pred)
        self.new_risk_cause_metadata.risk_lower = self.new_pair_model.exp_bounds[0]
        self.new_risk_cause_metadata.risk_upper = self.new_pair_model.exp_bounds[1]
        # save results
        self.new_risk_cause_metadata.to_csv(
            self.o_path / "new_risk_cause_metadata.csv",
            index=False
        )
        # plot figures
        plot_model(self.new_study_data,
                   self.new_output_data,
                   self.new_risk_cause_metadata,
                   self.new_pair_model)
        plt.savefig(self.o_path / "new_model_figure.pdf", bbox_inches="tight")
        plt.close("all")


class PostContinuousProcess(Process):
    """
    Post process for GBD 2020 continuous risk
    """

    def __init__(self, i_path: Path, o_path: Path):
        super().__init__(i_path, o_path)

        # extract the model specs
        self.pair_model = PairModel(
            self.name,
            self.signal_model,
            self.linear_model,
            exp_limits=(self.output_data.risk.min(),
                        self.output_data.risk.max())
        )

        # place holders
        for attr in ["new_linear_model",
                     "new_pair_model",
                     "new_output_data",
                     "new_risk_cause_metadata",
                     "pair_model_result",
                     "pair_model_prediction",
                     "score",
                     "pub_bias"]:
            setattr(self, attr, None)

    def run(self):
        super().run()
        self.get_new_linear_model()
        self.get_new_data()
        self.save_results()
        self.plot_models()

    def get_new_linear_model(self):
        # get data for new linear model
        df_signal = self.signal_model.data.to_df()
        df_signal["signal"] = self.pair_model.get_data_signal()
        df_signal["re_signal"] = self.pair_model.get_data_re_signal(use_new=True)
        df_signal.study_id = df_signal.study_id.astype(int)

        # get data frame for linear model
        df_linear = self.linear_model.data.to_df()
        df_linear.study_id = df_linear.study_id.astype(int)

        # match signal
        keys = np.unique(df_linear[["obs", "obs_se", "study_id"]].values, axis=0)
        df_new_linear = []
        for key in keys:
            sub_df_signal = df_signal[
                (df_signal.obs == key[0]) &
                (df_signal.obs_se == key[1]) &
                (df_signal.study_id == key[2])
            ].sort_values("signal")
            sub_df_linear = df_linear[
                (df_linear.obs == key[0]) &
                (df_linear.obs_se == key[1]) &
                (df_linear.study_id == key[2])
            ].sort_values("signal")
            for col in ["signal", "re_signal"]:
                sub_df_linear[col] = sub_df_signal[col].values[:sub_df_linear.shape[0]]
            df_new_linear.append(sub_df_linear)
        df_new_linear = pd.concat(df_new_linear)
#         print(f"df_signal: {df_signal.shape[0]}, df_new_linear: {df_new_linear.shape[0]}")
#         df_new_linear = df_signal.copy()
#         for name in self.pair_model.bias_covs:
#             df_new_linear[name] = 0.0

        data = MRData()
        data.load_df(
            df_new_linear,
            col_obs="obs",
            col_obs_se="obs_se",
            col_study_id="study_id",
            col_covs=["signal", "re_signal"] + self.pair_model.bias_covs
        )

        cov_models = [
            LinearCovModel("signal",
                           use_re=False,
                           prior_beta_uniform=[0.0, 1.0 if self.name == "sbp_cvd_ihd" else 5.0]),
            LinearCovModel("re_signal",
                           use_re=True,
                           prior_beta_uniform=[0.0, 0.0]),
            LinearCovModel("intercept",
                           use_re=True,
                           prior_beta_uniform=[0.0, 0.0])
        ]
        for bias_cov in self.pair_model.bias_covs:
            old_cov_model = self.linear_model.get_cov_model(bias_cov)
            cov_models.append(
                LinearCovModel(bias_cov,
                               use_re=False,
                               prior_beta_gaussian=old_cov_model.prior_beta_gaussian.copy())
            )
        self.new_linear_model = MRBRT(data, cov_models=cov_models)
        self.new_linear_model.fit_model(inner_print_level=5, inner_max_iter=1000)
        self.new_pair_model = PairModel(
            self.name,
            self.signal_model,
            self.new_linear_model,
            exp_limits=(self.output_data.risk.min(),
                        self.output_data.risk.max())
        )

    def get_new_data(self):
        # get new study data
        self.new_study_data = self.study_data.copy()
        pred = self.new_pair_model.get_pred(self.new_study_data.ref_risk.values)
        self.new_study_data["log_ref_cause"] = pred.log_cause.values
        self.new_study_data["log_alt_cause"] = (
            self.new_study_data.log_rr.values +
            self.new_study_data.log_ref_cause.values
        )
        self.new_study_data["log_alt_cause_lower"] = (
            self.new_study_data.log_alt_cause.values -
            1.96*self.new_study_data.log_rr_se.values
        )
        self.new_study_data["log_alt_cause_lower"] = (
            self.new_study_data.log_alt_cause.values +
            1.96*self.new_study_data.log_rr_se.values
        )

        log_col_names = ["log_ref_cause",
                         "log_alt_cause_lower",
                         "log_alt_cause",
                         "log_alt_cause_upper"]
        for log_col_name in log_col_names:
            linear_col_name = log_col_name.replace("log", "linear")
            self.new_study_data[linear_col_name] = np.exp(self.new_study_data[log_col_name])

        # get new output data
        self.new_output_data = self.output_data.copy()
        pred = self.new_pair_model.get_pred(self.output_data.risk.values)
        log_col_names = ["outer_log_cause_lower",
                         "inner_log_cause_lower",
                         "log_cause",
                         "inner_log_cause_upper",
                         "outer_log_cause_upper"]
        for log_col_name in log_col_names:
            linear_col_name = log_col_name.replace("log", "linear")
            self.new_output_data[log_col_name] = pred[log_col_name]
            self.new_output_data[linear_col_name] = np.exp(pred[log_col_name])

        # get new risk cause metadata
        self.new_risk_cause_metadata = self.risk_cause_metadata.copy()
        self.new_risk_cause_metadata.score = self.new_pair_model.get_score(pred)
        self.new_risk_cause_metadata.pub_bias = self.new_pair_model.get_pub_bias()
        self.new_risk_cause_metadata.risk_lower = self.new_pair_model.exp_bounds[0]
        self.new_risk_cause_metadata.risk_upper = self.new_pair_model.exp_bounds[1]

    def save_results(self):
        # copy old results
        for file_name in ["signal_model.pkl",
                          "linear_model.pkl",
                          "study_data.csv",
                          "output_data.csv",
                          "risk_cause_metadata.csv"]:
            shutil.copy(self.i_path / file_name, self.o_path / file_name)

        # save new results
        with open(self.o_path / "new_linear_model.pkl", "wb") as f:
            dill.dump(self.new_linear_model, f)

        for df_name in ["new_study_data",
                        "new_output_data",
                        "new_risk_cause_metadata"]:
            file_name = f"{df_name}.csv"
            getattr(self, df_name).to_csv(self.o_path / file_name, index=False)

    def plot_models(self):
        # plot original model
        plot_model(self.study_data,
                   self.output_data,
                   self.risk_cause_metadata,
                   self.pair_model)
        plt.savefig(self.o_path / "model_figure.pdf", bbox_inches="tight")

        # plot new model
        plot_model(self.new_study_data,
                   self.new_output_data,
                   self.new_risk_cause_metadata,
                   self.new_pair_model)
        plt.savefig(self.o_path / "new_model_figure.pdf", bbox_inches="tight")
        plt.close("all")

    def __repr__(self) -> str:
        return f"{type(self).__name__}(name={self.name})"


@dataclass
class PairModel:
    name: str
    signal_model: Union[MRBRT, MRBeRT]
    linear_model: MRBRT

    ref_covs: List[str] = field(init=False)
    alt_covs: List[str] = field(init=False)
    bias_covs: List[str] = field(init=False)

    exp_limits: Tuple[float] = field(default=())
    exp_bounds: Tuple[float] = field(init=False)

    j_shaped: bool = field(init=False)

    beta_info: Tuple[float] = field(init=False)
    gamma_info: Tuple[float] = field(init=False)
    residual_info: Tuple[ndarray] = field(init=False)
    outlier_index: ndarray = field(init=False)

    tmrel_info: Tuple[float] = field(init=False)

    def __post_init__(self):
        if isinstance(self.signal_model, MRBeRT):
            cov_model = self.signal_model.sub_models[0].cov_models[0]
        else:
            cov_model = self.signal_model.cov_models[0]
        self.ref_covs = cov_model.ref_cov
        self.alt_covs = cov_model.alt_cov

        self.bias_covs = [
            cov_model.name
            for cov_model in self.linear_model.cov_models
            if cov_model.name not in ["signal", "re_signal", "linear"]
        ]

        alt_exp = self.signal_model.data.get_covs(self.alt_covs)
        ref_exp = self.signal_model.data.get_covs(self.ref_covs)
        exp = np.hstack([alt_exp.ravel(), ref_exp.ravel()])
        alt_mid = alt_exp.mean(axis=1)
        ref_mid = ref_exp.mean(axis=1)
        self.alt_mid = alt_mid
        self.ref_mid = ref_mid

        if len(self.exp_limits) == 0:
            self.exp_limits = (exp.min(), exp.max())
        self.exp_bounds = (np.quantile(ref_mid, 0.15),
                           np.quantile(alt_mid, 0.85))

        self.j_shaped = (cov_model.prior_spline_monotonicity is None)

        self.beta_info = get_beta_info(self.linear_model)
        self.gamma_info = get_gamma_info(self.linear_model)

        signal = self.get_data_signal()
        re_signal = self.get_data_re_signal()
        residual = self.signal_model.data.obs - signal*self.beta_info[0]
        residual_sd = np.sqrt(self.signal_model.data.obs_se**2 +
                              re_signal**2*self.gamma_info[0])
        self.residual_info = (residual, residual_sd)
        if isinstance(self.signal_model, MRBeRT):
            trimming_weights = self.signal_model.get_w_soln()
        else:
            trimming_weights = self.signal_model.w_soln
        self.outlier_index = trimming_weights < 0.1

        # get tmrel infomation
        risk = np.linspace(*self.exp_limits, 100)
        signal = self.get_signal(risk)
        if self.j_shaped:
            self.tmrel_info = (risk[np.argmin(signal)], signal.min())
        else:
            indices = np.cumprod(np.abs(signal) < 1e-10).astype(bool)
            if indices.any():
                tmrel_index = np.arange(risk.size)[indices].max()
            else:
                tmrel_index = 0
            self.tmrel_info = (risk[tmrel_index], signal[tmrel_index])
        self.exp_bounds = (
            max(self.tmrel_info[0], self.exp_bounds[0]),
            self.exp_bounds[1]
        )
        if self.exp_bounds[0] > self.exp_bounds[1]:
            self.j_shaped = False
            indices = np.cumprod(np.abs(signal) < 1e-10).astype(bool)
            if indices.any():
                tmrel_index = np.arange(risk.size)[indices].max()
            else:
                tmrel_index = 0
            self.tmrel_info = (risk[tmrel_index], signal[tmrel_index])
            self.exp_bounds = (max(np.quantile(ref_mid, 0.15),
                                   self.tmrel_info[0]),
                               self.exp_bounds[1])

    def __repr__(self) -> str:
        return f"{type(self).__name__}(name={self.name})"

    @property
    def is_new(self) -> bool:
        return "re_signal" in self.linear_model.cov_names

    @property
    def is_locked(self) -> bool:
        # return "metab_bmi_adult" in self.name
        return False

    def get_data_signal(self) -> ndarray:
        return self.signal_model.predict(self.signal_model.data)

    def get_data_re_signal(self, use_new: bool = False) -> ndarray:
        # alt_mid = self.signal_model.data.get_covs(self.alt_covs).mean(axis=1)
        # ref_mid = self.signal_model.data.get_covs(self.ref_covs).mean(axis=1)
        # if (not self.is_locked) and (self.is_new or use_new):
        #     ref_mid.fill(self.exp_limits[0])
        # re_signal = alt_mid - ref_mid
        covs = deepcopy(self.signal_model.data.covs)
        if (not self.is_locked) and (self.is_new or use_new):
            for cov_name in self.ref_covs:
                covs[cov_name].fill(self.exp_limits[0])
        data = MRData(covs=covs)
        re_signal = self.signal_model.predict(data)

        return re_signal

    def get_signal(self, risk: ndarray) -> ndarray:
        covs = {}
        for cov_name in self.alt_covs:
            covs[cov_name] = risk
        for cov_name in self.ref_covs:
            covs[cov_name] = np.repeat(self.exp_limits[0], risk.size)
        data = MRData(covs=covs)
        return self.signal_model.predict(data)

    def get_re_signal(self, risk: ndarray) -> ndarray:
        # re_signal = risk - self.exp_limits[0]
        re_signal = self.get_signal(risk)
        return re_signal

    def get_pred(self, risk: ndarray) -> DataFrame:
        signal = self.get_signal(risk)
        re_signal = self.get_re_signal(risk)

        signal -= self.tmrel_info[1]
        re_signal -= self.tmrel_info[1]

        beta, beta_sd = self.beta_info
        gamma, gamma_sd = self.gamma_info
        outer_log_cause_sd = np.sqrt(
            signal**2*beta_sd**2 + re_signal**2*(gamma + 2*gamma_sd)
        )
        outer_log_cause_lower = signal*beta - 1.645*outer_log_cause_sd
        inner_log_cause_lower = signal*(beta - 1.645*beta_sd)
        log_cause = signal*beta
        inner_log_cause_upper = signal*(beta + 1.645*beta_sd)
        outer_log_cause_upper = signal*beta + 1.645*outer_log_cause_sd

        return pd.DataFrame({
            "risk": risk,
            "outer_log_cause_lower": outer_log_cause_lower,
            "inner_log_cause_lower": inner_log_cause_lower,
            "log_cause": log_cause,
            "inner_log_cause_upper": inner_log_cause_upper,
            "outer_log_cause_upper": outer_log_cause_upper
        })

    def get_score(self, pred: DataFrame) -> float:
        risk = pred.risk.values
        log_cause = pred.log_cause.values
        outer_log_cause_lower = pred.outer_log_cause_lower.values
        outer_log_cause_upper = pred.outer_log_cause_upper.values

        score_index = ((risk >= self.exp_bounds[0]) &
                       (risk <= self.exp_bounds[1]))
        score_sign = np.sign(log_cause[score_index].mean())
        return min(outer_log_cause_lower[score_index].mean()*score_sign,
                   outer_log_cause_upper[score_index].mean()*score_sign)

    def get_pub_bias(self) -> int:
        return get_pub_bias(*self.residual_info)


def plot_model(study_data: DataFrame,
               output_data: DataFrame,
               risk_cause_metadata: DataFrame,
               pair_model: PairModel,
               ax: List[plt.Axes] = None):
    if ax is None:
        _, ax = plt.subplots(1, 2, figsize=(16, 5))

    # plot data
    ax[0].scatter(
        study_data.alt_risk,
        study_data.log_alt_cause,
        s=5.0/study_data.log_rr_se,
        color="gray",
        alpha=0.5
    )
    outlier_index = study_data.is_outlier == 1
    ax[0].scatter(
        study_data.alt_risk[outlier_index],
        study_data.log_alt_cause[outlier_index],
        s=5.0/study_data.log_rr_se[outlier_index],
        color="red",
        alpha=0.5,
        marker="x"
    )

    # plot prediction
    ax[0].plot(output_data.risk, output_data.log_cause,
               color="#008080", linewidth=1)

    # plot uncertainties
    ax[0].fill_between(output_data.risk,
                       output_data.inner_log_cause_lower,
                       output_data.inner_log_cause_upper,
                       color="#008080",
                       alpha=0.2)
    ax[0].fill_between(output_data.risk,
                       output_data.outer_log_cause_lower,
                       output_data.outer_log_cause_upper,
                       color="#008080",
                       alpha=0.2)

    # plot bounds
    for b in [risk_cause_metadata.risk_lower.values[0],
              risk_cause_metadata.risk_upper.values[0]]:
        ax[0].axvline(b, linestyle="--", linewidth=1, color="k")

    # plot 0 line
    ax[0].axhline(0.0, linestyle="-", linewidth=1, color="k")

    # add unit to the xaxis
    ax[0].set_xlabel(risk_cause_metadata.risk_unit.values[0])

    # title
    title = (f"name={pair_model.name}, "
             f"score={risk_cause_metadata.score.values[0]: .3f}")
    ax[0].set_title(title, loc="left")

    # plot residual
    residual, residual_sd = pair_model.residual_info
    ax[1].set_ylim(residual_sd.max(), 0.0)
    ax[1].scatter(residual, residual_sd,
                  color="gray", alpha=0.4)
    ax[1].scatter(residual[pair_model.outlier_index],
                  residual_sd[pair_model.outlier_index],
                  color="red", alpha=0.4, marker="x")
    ax[1].fill_betweenx(
        [0.0, residual_sd.max()],
        [0.0, -1.96*residual_sd.max()],
        [0.0, 1.96*residual_sd.max()],
        color="#B0E0E6", alpha=0.4
    )
    ax[1].plot([0.0, -1.96*residual_sd.max()],
               [0.0, residual_sd.max()],
               linewidth=1, color="#87CEFA")
    ax[1].plot([0.0, 1.96*residual_sd.max()],
               [0.0, residual_sd.max()],
               linewidth=1, color="#87CEFA")
    ax[1].axvline(0.0, color="k", linewidth=1, linestyle="--")
    ax[1].set_xlabel("residual")
    ax[1].set_ylabel("residual sd")

    
def plot_model_publication(study_data: DataFrame,
                           output_data: DataFrame,
                           risk_cause_metadata: DataFrame,
                           pair_model: PairModel,
                           rei: str,
                           rei_name: str,
                           acause: str,
                           cause_name: str,
                           ax: List[plt.Axes] = None):
    if ax is None:
        _, ax = plt.subplots(1, 2, figsize=(16, 5))

    # plot data
    ax[0].scatter(
        study_data.alt_risk,
        study_data.log_alt_cause,
        s=5.0/study_data.log_rr_se,
        color="gray",
        alpha=0.5
    )
    outlier_index = study_data.is_outlier == 1
    ax[0].scatter(
        study_data.alt_risk[outlier_index],
        study_data.log_alt_cause[outlier_index],
        s=5.0/study_data.log_rr_se[outlier_index],
        color="red",
        alpha=0.5,
        marker="x"
    )

    # plot prediction
    ax[0].plot(output_data.risk, output_data.log_cause,
               color="#008080", linewidth=1.2)

    # plot uncertainties
    ax[0].fill_between(output_data.risk,
                       output_data.inner_log_cause_lower,
                       output_data.inner_log_cause_upper,
                       color="#008080",
                       alpha=0.2)
    ax[0].fill_between(output_data.risk,
                       output_data.outer_log_cause_lower,
                       output_data.outer_log_cause_upper,
                       color="#008080",
                       alpha=0.2)
    log_cause_bound = output_data.outer_log_cause_lower if output_data.log_cause.mean() >= 0 else output_data.outer_log_cause_upper
    ax[0].plot(output_data.risk,
               log_cause_bound,
               color="red",
               linewidth=1.2, alpha=0.9)
    ax[0].plot(output_data.risk,
               log_cause_bound,
               color="gray",
               linewidth=1.2, alpha=0.5)
    lb = output_data.risk.min()
    ub = study_data.alt_risk.max()
    ax[0].set_xlim(lb - (ub - lb)*0.05, ub + (ub - lb)*0.1)

    # plot bounds
    for b in [risk_cause_metadata.risk_lower.values[0],
              risk_cause_metadata.risk_upper.values[0]]:
        ax[0].axvline(b, linestyle="--", linewidth=1, color="gray")

    # plot 0 line
    ax[0].axhline(0.0, linestyle="-", linewidth=1, color="gray")

    # add unit to the xaxis
    ax[0].set_xlabel(f"{rei} ({risk_cause_metadata.risk_unit.values[0]})")
    ax[0].set_ylabel("log relative risk")

    # title
    title = (f"{rei_name} / {cause_name}")
    ax[0].set_title(title, loc="left")

    # plot residual
    residual, residual_sd = pair_model.residual_info
    residual_sd_max = residual_sd.max()*1.1
    ax[1].set_ylim(residual_sd_max, 0.0)
    ax[1].scatter(residual, residual_sd,
                  color="gray", alpha=0.4)
    ax[1].scatter(residual[pair_model.outlier_index],
                  residual_sd[pair_model.outlier_index],
                  color="red", alpha=0.4, marker="x")
    ax[1].fill_betweenx(
        [0.0, residual_sd_max],
        [0.0, -1.96*residual_sd_max],
        [0.0, 1.96*residual_sd_max],
        color="#B0E0E6", alpha=0.4
    )
    ax[1].plot([0.0, -1.96*residual_sd_max],
               [0.0, residual_sd_max],
               linewidth=1, color="#87CEFA")
    ax[1].plot([0.0, 1.96*residual_sd_max],
               [0.0, residual_sd_max],
               linewidth=1, color="#87CEFA")
    ax[1].axvline(0.0, color="gray", linewidth=1, linestyle="--")
    ax[1].set_xlabel("residual")
    ax[1].set_ylabel("residual sd")