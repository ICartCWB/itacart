"""Tests for :mod:`itacart.geometry`."""

from __future__ import annotations

import itacart  # noqa: F401

# No tests are portable from itacart_core: it has no coverage for this module.
# The 406 figure that used to stand here was the size of its whole suite, not a
# count of anything reusable. Measured in F6 with
# `grep -rlE "neighbor|grid_disk|adjacen" itacart_core/`, which returned nothing.
