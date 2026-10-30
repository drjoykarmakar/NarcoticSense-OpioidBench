#!/usr/bin/env python
"""Thin wrapper for importing a legally obtained source file."""

import sys

from opioidbench.cli import main

raise SystemExit(main(["data", "import", *sys.argv[1:]]))
