.PHONY: setup sources download demo dataset enrich resolve taxonomy validate test lint reproduce

setup:
	python -m pip install -e '.[dev]'

sources:
	opioidbench data sources

download:
	@echo "Automatic SWGDRUG download is disabled by source policy. Obtain a legally permitted archive and use: opioidbench data import --source swgdrug_eims --input /path/to/archive.zip --destination data/raw/swgdrug_3.14 --retrieval-date YYYY-MM-DD"

demo:
	opioidbench demo --fixture tests/fixtures/synthetic_eims.jdx --output results/demo

dataset:
	@test -n "$(INPUT)" || (echo "Usage: make dataset INPUT=/path/to/library.hpj" && exit 2)
	opioidbench data build --input "$(INPUT)" --source swgdrug_eims --output data/processed/opioidbench-eims-v0.1.0 $(if $(SOURCE_ARCHIVE),--source-archive "$(SOURCE_ARCHIVE)",) $(if $(COMPOUND_LIST),--compound-list "$(COMPOUND_LIST)",) $(if $(EXACT_MASS_LIST),--exact-mass-list "$(EXACT_MASS_LIST)",) $(if $(AGILENT_ARCHIVE),--agilent-archive "$(AGILENT_ARCHIVE)",) $(if $(NIST_ARCHIVE),--nist-archive "$(NIST_ARCHIVE)",)

enrich:
	opioidbench metadata enrich --input data/processed/opioidbench-eims-v0.1.0/compounds.csv --output data/processed/opioidbench-eims-v0.1.0/enrichment --cache data/external/pubchem_cache

resolve:
	@test -n "$(CHEBI_SDF)" || (echo "Usage: make resolve CHEBI_SDF=/path/to/chebi.sdf.gz" && exit 2)
	opioidbench metadata resolve --input data/processed/opioidbench-eims-v0.1.0/compounds.csv --chebi-sdf "$(CHEBI_SDF)" --output data/processed/opioidbench-eims-v0.1.0/resolution --cache data/external/pubchem_cache

taxonomy:
	opioidbench metadata curate-taxonomy --input data/processed/opioidbench-eims-v0.1.0/resolution/resolved_compounds.csv --evidence metadata/opioid_evidence.csv --exclusions metadata/taxonomy_exclusions.csv --output metadata/opioid_taxonomy.csv --review-output results/private/opioid_taxonomy_review.csv

validate:
	opioidbench data validate --input tests/fixtures/synthetic_eims.jdx --source synthetic_fixture

test:
	python -m pytest

lint:
	ruff check .
	ruff format --check .

reproduce: test lint demo
	@echo "Synthetic CI reproduction complete. The v0.4.0 real-source pilot requires legally obtained local data."
