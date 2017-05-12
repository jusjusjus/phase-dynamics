
import logging
import numpy as np
from .futils.futils import futils


logger = logging.getLogger(name="phd.core")

period = np.inf

def set_period(p):
    """Set the periodicity of periodic variables."""
    global period
    assert p > 0.0, "Period must be larger than 0 (gave %.3f)." % (p)
    period = p
    futils.set_period(p)

set_period(2.0 * np.pi)


def mod(x):
    return np.mod(x, period)


def unmod(x):
    """Unwraps a phase variable."""

    x = np.asarray(x, dtype=np.float64)
    return futils.unmod(x)


def threshold_data(x, threshold, n_min=1):
    """Computes list of super-threshold slices.
    
    Each slice has length > n_min, and within each slice x > threshold.
    """

    assert n_min > 0
    x = np.asarray(x, dtype=np.float64)

    segments, num_segments = futils.threshold_data(x, threshold, n_min)
    slices = [
            slice(seg[0], seg[1])
            for seg in segments[:num_segments]
    ]
    return slices


def poincare_times(x, x0=0.0, interp=True):
    x     = unmod(x)-x0
    x_mod = mod(x)
    idx, ti, n_idx = futils.poincare_times_no_false_returns(x_mod, x)
    if n_idx == 0:
        logger.info("no crossings found.")
        return np.array([]), np.array([])
    idx = idx[:n_idx] # -1 because fortran
    ti  = ti[:n_idx]

    if interp == True:
        idx, ti = futils.poincare_times_interpolate(idx, ti, x_mod)

    idx -= 1 # -1 because fortran
    return idx, ti[:-1]

