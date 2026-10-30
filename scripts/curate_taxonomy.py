#!/usr/bin/env python
"""Curate exact evidence-backed opioid taxonomy labels."""

import sys

from opioidbench.cli import main

raise SystemExit(main(["metadata", "curate-taxonomy", *sys.argv[1:]]))
