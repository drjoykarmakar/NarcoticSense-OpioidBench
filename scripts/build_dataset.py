#!/usr/bin/env python
"""Build a versioned normalized dataset from a local EI-MS source file."""

import sys

from opioidbench.cli import main

raise SystemExit(main(["data", "build", *sys.argv[1:]]))
