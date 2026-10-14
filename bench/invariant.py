# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""The bench scripts' one assertion: Python's `assert` vanishes under -O; this never does
(lessons log, 8 Oct 2026). Import it with `sys.path` pointing at bench/."""


def invariant(condition, message):
    """Halts on a broken invariant, saying which."""
    if not condition:
        raise AssertionError(message)
