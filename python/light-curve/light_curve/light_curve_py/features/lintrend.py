import numpy as np

from ._base import BaseFeature
from ._lstsq import least_squares


class LinearTrend(BaseFeature):
    def _eval(self, t, m, sigma=None):
        n = len(t)

        if n == 2:
            return (m[1] - m[0]) / (t[1] - t[0]), 0, 0

        slope, chi2 = least_squares(t, m, None)

        red_chi2 = chi2 / (n - 2)
        sxx = np.var(t, ddof=n - 1)
        return np.array([slope, np.sqrt(red_chi2 / sxx), red_chi2])

    @property
    def size(self):
        return 3


__all__ = ("LinearTrend",)
