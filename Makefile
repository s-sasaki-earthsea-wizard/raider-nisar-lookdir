.DEFAULT_GOAL := help
PY := .venv/bin/python

.PHONY: help setup verify verify-local clean

help:  ## Show this help
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

setup:  ## Create the isolated venv and install deps with uv
	uv venv --python 3.12 .venv
	uv pip install --python $(PY) -r pyproject.toml

verify:  ## Stream NISAR metadata, re-derive look side, write reports/lookdir_verification.md
	$(PY) verify_lookdir.py

verify-local:  ## Same, but read pre-downloaded granules from data/ (DIR=... to override)
	$(PY) verify_lookdir.py --local $(or $(DIR),data)

clean:  ## Remove the venv and caches (keeps reports/)
	rm -rf .venv __pycache__
