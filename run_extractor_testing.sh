#!/usr/bin/env bash
# Run the chemistry extractor in testing mode using conda run

cd "$(dirname "$0")"
cd src
conda run -n chem_env python3 -m domains.chemistry.main --env testing