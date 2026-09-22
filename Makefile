.PHONY: bootstrap test test-qe audit pseudos clean-scratch

bootstrap:        ## one-time local setup (venv + QE env)
	bash scripts/bootstrap.sh

test:             ## fast unit tests (no QE)
	.venv/bin/python -m pytest

test-qe:          ## includes QE integration tests
	.venv/bin/python -m pytest -m "qe"

audit:            ## environment audit -> results/provenance/environment.json
	bash scripts/check_environment.sh

pseudos:          ## fetch + validate pseudopotentials, write manifest
	.venv/bin/python scripts/fetch_pseudopotentials.py

clean-scratch:    ## delete QE scratch only (keeps inputs/results)
	find calculations -type d -name out -prune -exec rm -rf {} +
	rm -rf work/* .tmp/*
