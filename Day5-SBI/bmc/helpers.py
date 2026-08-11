import numpy as np
import bayesflow as bf

from ssms.basic_simulators.simulator import simulator as ssm_simulator
from ssms.config import model_config as ssms_model_config


# Model configurations
lba_cfg = ssms_model_config["lba3"]
param_names_lba = lba_cfg["params"]
lower_lba = np.array(lba_cfg["param_bounds"][0])
upper_lba = np.array(lba_cfg["param_bounds"][1])

lca_cfg = ssms_model_config["lca_3"]
param_names_lca = lca_cfg["params"]
lower_lca = np.array(lca_cfg["param_bounds"][0])
upper_lca = np.array(lca_cfg["param_bounds"][1])

lca_no_bias_cfg = ssms_model_config["lca_no_z_3"]
param_names_lca_no_bias = lca_no_bias_cfg["params"]
lower_lca_no_bias = np.array(lca_no_bias_cfg["param_bounds"][0])
upper_lca_no_bias = np.array(lca_no_bias_cfg["param_bounds"][1])


def lba_prior():
    return {
        name: np.random.uniform(lo, hi)
        for name, lo, hi in zip(param_names_lba, lower_lba, upper_lba)
    }


def lba_sim(A, b, v0, v1, v2, num_trials=100):
    """Simulate multiple (RT, choice) trials for LBA."""
    result = ssm_simulator(
        theta={"A": A, "b": b, "v0": v0, "v1": v1, "v2": v2},
        model="lba3",
        n_samples=num_trials,
        delta_t=0.001,
    )
    obs = np.array([result["rts"][:, 0], result["choices"][:, 0]]).T.astype(
        "float32", copy=False
    )
    return {"obs": obs}


def lca_prior():
    return {
        name: np.random.uniform(lo, hi)
        for name, lo, hi in zip(param_names_lca, lower_lca, upper_lca)
    }


def lca_sim(v0, v1, v2, a, z0, z1, z2, g, b, t, num_trials=100):
    """Simulate multiple (RT, choice) trials for LCA with bias."""
    result = ssm_simulator(
        theta={
            "v0": v0,
            "v1": v1,
            "v2": v2,
            "a": a,
            "z0": z0,
            "z1": z1,
            "z2": z2,
            "g": g,
            "b": b,
            "t": t,
        },
        model="lca_3",
        n_samples=num_trials,
        delta_t=0.001,
    )
    obs = np.array([result["rts"][:, 0], result["choices"][:, 0]]).T.astype(
        "float32", copy=False
    )
    return {"obs": obs}


def lca_no_bias_prior():
    return {
        name: np.random.uniform(lo, hi)
        for name, lo, hi in zip(
            param_names_lca_no_bias, lower_lca_no_bias, upper_lca_no_bias
        )
    }


def lca_no_bias_sim(v0, v1, v2, a, g, b, t, num_trials=100):
    """Simulate multiple (RT, choice) trials for LCA without starting-point bias."""
    result = ssm_simulator(
        theta={"v0": v0, "v1": v1, "v2": v2, "a": a, "g": g, "b": b, "t": t},
        model="lca_no_z_3",
        n_samples=num_trials,
        delta_t=0.001,
    )
    obs = np.array([result["rts"][:, 0], result["choices"][:, 0]]).T.astype(
        "float32", copy=False
    )
    return {"obs": obs}


def build_meta_model(num_trials=100):
    def _lba_sim_fixed(A, b, v0, v1, v2):
        return lba_sim(A, b, v0, v1, v2, num_trials=num_trials)

    def _lca_sim_fixed(v0, v1, v2, a, z0, z1, z2, g, b, t):
        return lca_sim(v0, v1, v2, a, z0, z1, z2, g, b, t, num_trials=num_trials)

    def _lca_no_bias_sim_fixed(v0, v1, v2, a, g, b, t):
        return lca_no_bias_sim(v0, v1, v2, a, g, b, t, num_trials=num_trials)

    simulators = [
        bf.make_simulator([lba_prior, _lba_sim_fixed]),
        bf.make_simulator([lca_prior, _lca_sim_fixed]),
        bf.make_simulator([lca_no_bias_prior, _lca_no_bias_sim_fixed]),
    ]
    return bf.simulators.ModelComparisonSimulator(simulators)
