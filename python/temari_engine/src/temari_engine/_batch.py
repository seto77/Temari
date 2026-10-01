"""Shared validation helpers for the NumPy batch readers."""

import numbers

import numpy as np


def _is_bool(value):
    return isinstance(value, (bool, np.bool_))


def as_float64_1d(values, name, scalar_guard):
    """Convert a one-dimensional real array-like using the scalar API's guard.

    Numeric NumPy arrays take the vectorized path. Object arrays are checked one
    item at a time so bools and non-real objects cannot sneak through coercion.
    """
    if isinstance(values, (list, tuple)) and not isinstance(values, np.ndarray):
        for value in values:
            if _is_bool(value):
                raise TypeError("%s は実数 (bool 不可): %r" % (name, value))
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError("%s は1次元の実数配列" % name) from exc
    if raw.ndim != 1:
        raise TypeError("%s は1次元配列: shape=%r" % (name, raw.shape))
    if raw.dtype.kind in "iuf":
        out = np.asarray(raw, dtype=np.float64)
        bad = np.flatnonzero(~np.isfinite(out))
        if bad.size:
            scalar_guard(float(out[int(bad[0])]))
        return out
    if raw.dtype.kind != "O":
        raise TypeError("%s は実数配列 (bool・文字列・複素数は不可)" % name)

    out = np.empty(raw.shape, dtype=np.float64)
    for i, value in enumerate(raw):
        if _is_bool(value) or not isinstance(value, numbers.Real):
            raise TypeError("%s は実数 (bool 不可): %r" % (name, value))
        out[i] = scalar_guard(value)
    return out
