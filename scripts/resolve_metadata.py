#!/usr/bin/env python
"""Resolve chemical identities with local ChEBI and PubChem fallback."""

import sys

from opioidbench.cli import main

raise SystemExit(main(["metadata", "resolve", *sys.argv[1:]]))
