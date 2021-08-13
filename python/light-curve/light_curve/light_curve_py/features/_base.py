import numpy as np

from abc import ABC, abstractmethod


class BaseFeature(ABC):
    @staticmethod
    def _normalize_input(t, m, sigma, sorted):
        t = np.asarray(t)
        m = np.asarray(m)
        if sigma is not None:
            sigma = np.asarray(sigma)

        if sorted is None:
            if np.any(np.diff(t) <= 0):
                raise ValueError("t must be sorted")
        elif not sorted:
            idx = np.argsort(t)
            t = t[idx]
            m = m[idx]
            if sigma is not None:
                sigma = sigma[idx]

        return t, m, sigma

    def __call__(self, t, m, sigma=None, sorted=None, fill_value=None):
        t, m, sigma = self._normalize_input(t, m, sigma, sorted)
        return self._eval(t, m, sigma)

    @abstractmethod
    def _eval(self, t, m, sigma=None):
        pass
