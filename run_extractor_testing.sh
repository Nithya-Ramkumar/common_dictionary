#!/usr/bin/env bash
# Activate the virtual environment and run the chemistry extractor in testing mode
 
cd "$(dirname "$0")"
source .venv/bin/activate
cd src
python3 -m domains.chemistry.main --env testing --expected-dir test_config1 