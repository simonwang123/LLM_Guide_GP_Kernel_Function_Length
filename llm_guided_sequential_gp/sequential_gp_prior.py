"""Soft log-length-scale prior used by the sequential GP MAP fits."""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from sklearn.exceptions import ConvergenceWarning


CATEGORY_LOG_SHIFT = {"short": -math.log(2.0), "base": 0.0, "long": math.log(2.0)}


@dataclass(frozen=True)
class PriorConfig:
    """Pilot settings fixed across rule candidates; the LLM selects categories."""

    base_factor: float = 0.5
    log_std: float = math.log(2.0)
    length_scale_lower: float = 0.01
    length_scale_upper: float = 100.0
    noise_lower: float = 1e-4
    noise_upper: float = 1e-1
    noise_initial: float = 1e-3

    def validate(self, dimension: int) -> None:
        if type(dimension) is not int or dimension < 1:
            raise ValueError("prior dimension must be a positive integer")
        values = (self.base_factor, self.log_std, self.length_scale_lower,
                  self.length_scale_upper, self.noise_lower, self.noise_upper,
                  self.noise_initial)
        if any(type(value) not in (int, float) or not math.isfinite(value) or value <= 0
               for value in values):
            raise ValueError("prior settings must be finite and positive")
        if self.length_scale_lower >= self.length_scale_upper:
            raise ValueError("length-scale optimizer bounds must increase")
        if self.noise_lower >= self.noise_upper:
            raise ValueError("noise optimizer bounds must increase")
        if not self.noise_lower < self.noise_initial < self.noise_upper:
            raise ValueError("noise initial value must lie inside optimizer bounds")
        base = self.base_factor * math.sqrt(dimension)
        if not (self.length_scale_lower < base / 2
                and base * 2 < self.length_scale_upper):
            raise ValueError("all prior centers must lie inside the optimizer bounds")


def log_centers(categories: dict[str, str], variables: tuple[str, ...],
                config: PriorConfig) -> np.ndarray:
    config.validate(len(variables))
    if not isinstance(categories, dict) or set(categories) != set(variables):
        raise ValueError("one prior category is required for each variable")
    if any(type(value) is not str or value not in CATEGORY_LOG_SHIFT
           for value in categories.values()):
        raise ValueError("prior categories must be short, base, or long")
    baseline = math.log(config.base_factor) + 0.5 * math.log(len(variables))
    return np.asarray([baseline + CATEGORY_LOG_SHIFT[categories[name]]
                       for name in variables], dtype=float)


def lengthscale_theta_indices(kernel: object, dimension: int) -> np.ndarray:
    """Find Matérn ARD entries in sklearn's log-parameter vector."""
    offset = 0
    for parameter in kernel.hyperparameters:
        if parameter.fixed:
            continue
        if parameter.name == "k1__k2__length_scale":
            if parameter.n_elements != dimension:
                raise ValueError("Matérn kernel does not have the expected ARD dimension")
            return np.arange(offset, offset + dimension, dtype=int)
        offset += parameter.n_elements
    raise ValueError("kernel has no trainable ARD Matérn length scales")


class LengthscaleMAPOptimizer:
    """Minimize negative GP log marginal likelihood plus a Gaussian log prior."""

    def __init__(self, centers: np.ndarray, log_std: float, indices: np.ndarray):
        self.centers = np.asarray(centers, dtype=float)
        self.log_std = float(log_std)
        self.indices = np.asarray(indices, dtype=int)
        if (self.centers.ndim != 1 or self.indices.shape != self.centers.shape
                or self.log_std <= 0 or not math.isfinite(self.log_std)
                or not np.isfinite(self.centers).all()
                or np.any(self.indices < 0) or len(np.unique(self.indices)) != len(self.indices)):
            raise ValueError("invalid length-scale MAP optimizer settings")

    def penalty(self, theta: np.ndarray) -> float:
        deviation = (np.asarray(theta)[self.indices] - self.centers) / self.log_std
        return float(0.5 * (deviation @ deviation))

    def __call__(self, obj_func, initial_theta: np.ndarray,
                 bounds: np.ndarray) -> tuple[np.ndarray, float]:
        def objective(theta: np.ndarray) -> tuple[float, np.ndarray]:
            value, gradient = obj_func(theta, eval_gradient=True)
            deviation = theta[self.indices] - self.centers
            derivative = np.asarray(gradient, dtype=float).copy()
            derivative[self.indices] += deviation / self.log_std**2
            penalty = 0.5 * (deviation @ deviation) / self.log_std**2
            return float(value + penalty), derivative

        result = minimize(objective, initial_theta, method="L-BFGS-B", jac=True,
                          bounds=bounds)
        if not np.isfinite(result.fun) or not np.isfinite(result.x).all():
            raise ValueError("MAP optimizer produced a non-finite result")
        if not result.success:
            warnings.warn(f"MAP optimizer did not converge: {result.message}",
                          ConvergenceWarning, stacklevel=2)
        return np.asarray(result.x, dtype=float), float(result.fun)
