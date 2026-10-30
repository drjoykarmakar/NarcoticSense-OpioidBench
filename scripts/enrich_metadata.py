#!/usr/bin/env python
"""Enrich compound metadata from cache or explicitly enabled PubChem access."""

import sys

from opioidbench.cli import main

raise SystemExit(main(["metadata", "enrich", *sys.argv[1:]]))
