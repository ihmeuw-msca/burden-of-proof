"""
Main GBD Evience Score Pipeline
"""
import matplotlib.pyplot as plt

from espipeline.process import Process


class PostDichotomousProcess(Process):
    """
    Post process for GBD 2020 dichotomous risk
    """

    def plot_model(self):
        _, ax = plt.subplots(figsize=(8, 5))

        # plot data
        ax.scatter(self.study_data.log_rr,
                   self.study_data.log_rr_se,
                   color="grey", alpha=0.4)
        outlier_index = self.study_data.is_outlier == 1
        ax.scatter(self.study_data.log_rr[outlier_index],
                   self.study_data.log_rr_se[outlier_index],
                   color="grey", alpha=0.4)

        # plot funnel
        beta = self.model.fe_soln["intercept"][0]
        se_max = self.study_data.log_rr_se.max()
        ax.fill_betweenx(
            [0.0, se_max],
            [beta, beta - 1.96*se_max],
            [beta, beta + 1.96*se_max],
            color="#B0E0E6", alpha=0.4
        )
        ax.plot([beta, beta - 1.96*se_max],
                [0.0, se_max],
                linewidth=1, color="#87CEFA")
        ax.plot([beta, beta + 1.96*se_max],
                [0.0, se_max],
                linewidth=1, color="#87CEFA")
        ax.set_ylim([se_max, 0.0])

        # plot vertical lines
        ax.axvline(0.0, color="k")
        ax.axvline(beta, color="#008080")
        ax.fill_betweenx([0.0, se_max],
                         [self.output_data.outer_log_cause_lower.values[0]]*2,
                         [self.output_data.outer_log_cause_upper.values[0]]*2,
                         color="#008080", alpha=0.2)
        ax.fill_betweenx([0.0, se_max],
                         [self.output_data.inner_log_cause_lower.values[0]]*2,
                         [self.output_data.inner_log_cause_upper.values[0]]*2,
                         color="#008080", alpha=0.2)

        # set title and labels
        title = (f"name={self.name}, "
                 f"score={self.risk_cause_metadata.score.values[0]: .3f}")
        ax.set_title(title, loc="left")
        ax.set_xlabel("log_rr")
        ax.set_ylabel("log_rr_se")

        plt.savefig(self.o_path / "model_figure.pdf", bbox_inches="tight")
        plt.close("all")

    def run(self):
        super().run()
        self.plot_model()
