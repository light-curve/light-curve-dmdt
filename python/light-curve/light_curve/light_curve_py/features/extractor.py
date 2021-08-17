from dataclasses import dataclass
from typing import Collection, Union

import numpy as np

from ._base import BaseFeature
from light_curve.light_curve_ext import Extractor as _RustExtractor, _FeatureEvaluator as _RustBaseFeature


@dataclass()
class _PyExtractor(BaseFeature):
    features: Collection[Union[BaseFeature, _RustBaseFeature]] = ()

    def _eval(self, t, m, sigma=None):
        raise NotImplementedError("_eval is missed for _PyExtractor")

    def _eval_and_fill(self, t, m, sigma, fill_value):
        return np.concatenate(
            [np.atleast_1d(feature._eval_and_fill(t, m, sigma, fill_value)) for feature in self.features]
        )

    @property
    def size(self):
        return sum(feature.size for feature in self.features)


class Extractor:
    def __new__(cls, features: Collection[Union[BaseFeature, _RustBaseFeature]]):
        if len(features) > 0 and all(isinstance(feature, _RustBaseFeature) for feature in features):
            return _RustExtractor(features)
        else:
            return _PyExtractor(features)


__all__ = ("Extractor",)
