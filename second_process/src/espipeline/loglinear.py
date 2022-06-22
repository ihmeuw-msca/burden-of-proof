"""
Main GBD Evience Score Pipeline
"""
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame

from mrtool import MRBRT

from espipeline.process import Process


class PostLogLinearProcess(Process):
    """
    Post process for GBD 2020 loglinear risk
    """

    def plot_models(self):
        # plot original model
        plot_model(self.study_data,
                   self.output_data,
                   self.risk_cause_metadata,
                   self.model,
                   self.name)
        plt.savefig(self.o_path / "model_figure.pdf", bbox_inches="tight")
        plt.close("all")

    def run(self):
        super().run()
        self.plot_models()


def get_data_signal(model: MRBRT) -> np.ndarray:
    for cov_model in model.cov_models:
        if len(cov_model.ref_cov) != 0:
            break
    alt_cov = model.data.get_covs(cov_model.alt_cov).mean(axis=1)
    ref_cov = model.data.get_covs(cov_model.ref_cov).mean(axis=1)
    return alt_cov - ref_cov


def plot_model(study_data: DataFrame,
               output_data: DataFrame,
               risk_cause_metadata: DataFrame,
               model: MRBRT,
               name: str,
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
    title = (f"name={name}, "
             f"score={risk_cause_metadata.score.values[0]: .3f}")
    ax[0].set_title(title, loc="left")

    # plot residual
    residual = model.data.obs - model.predict(model.data)
    residual_sd = np.sqrt(model.data.obs_se**2 + get_data_signal(model)**2*model.gamma_soln[0])
    outlier_index = model.w_soln < 0.1
    ax[1].set_ylim(residual_sd.max(), 0.0)
    ax[1].scatter(residual, residual_sd,
                  color="gray", alpha=0.4)
    ax[1].scatter(residual[outlier_index],
                  residual_sd[outlier_index],
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
