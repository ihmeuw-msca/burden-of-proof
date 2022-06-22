"""
Utility Functions
"""
from itertools import chain
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from mrtool import MRBRT
from mrtool.core.other_sampling import (extract_simple_lme_hessian,
                                        extract_simple_lme_specs)
from numpy import ndarray
from scipy.stats import norm


def get_fe_hessian(model: MRBRT) -> ndarray:
    specs = extract_simple_lme_specs(model)
    return extract_simple_lme_hessian(specs)


def get_re_fisher(model: MRBRT) -> ndarray:
    lt = model.lt
    return lt.get_gamma_fisher(lt.gamma)


def get_beta_info(model: MRBRT, name: str = "signal") -> Tuple[float]:
    # get beta solution
    cov_index = model.cov_names.index(name)
    beta = model.beta_soln[cov_index]
    beta_hessian = get_fe_hessian(model)
    beta_sd = 1.0/np.sqrt(np.diag(beta_hessian))[cov_index]
    return (beta, beta_sd)


def get_gamma_info(model: MRBRT) -> Tuple[float]:
    # get gamma solution
    gamma = model.gamma_soln[0]
    gamma_fisher = get_re_fisher(model)
    gamma_sd = 1.0/np.sqrt(gamma_fisher[0, 0])
    return (gamma, gamma_sd)


def get_pval(mean, sd, one_sided: bool = False) -> float:
    zscore = np.abs(mean/sd)
    if one_sided:
        pval = 1 - norm.cdf(zscore)
    else:
        pval = 2*(1 - norm.cdf(zscore))
    return pval


def egger_regression(residual, residual_sd, one_sided: bool = True) -> Dict[str, float]:
    weighted_residual = residual/residual_sd
    r_mean = weighted_residual.mean()
    r_sd = 1/np.sqrt(weighted_residual.size)
    r_pval = get_pval(r_mean, r_sd, one_sided=one_sided)
    return {
        "mean": r_mean,
        "sd": r_sd,
        "pval": r_pval
    }


def get_pub_bias(*args, **kwargs) -> int:
    result = egger_regression(*args, **kwargs)
    return int(result["pval"] < 0.05)


def list_all_files(path: Path) -> List[Path]:
    if path.is_file():
        return [path]
    return chain.from_iterable(
        list_all_files(sub_path) for sub_path in path.iterdir()
    )
