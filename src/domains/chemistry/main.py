import argparse
import sys
from domains.chemistry.extractor import ExtractorOrchestrator
from config.env_loader import EnvironmentLoader
from config_reconciliation import run_config_reconciliation

# ===============================================
# **CHEMISTRY EXTRACTOR MAIN ENTRYPOINT USAGE**
#
# How to specify the expected_dir (config directory for validation):
#
# 1. Command-Line Argument (Recommended):
#    python3 main.py --env testing --expected-dir test_config1
#    - This uses 'test_config1' as the directory for config validation.
#    - The --expected-dir argument takes precedence over the environment variable.
#
# 2. Environment Variable:
#    export EXPECTED_CONFIG_DIR=test_config1
#    python3 main.py --env testing
#    - If --expected-dir is not provided, the script will use the EXPECTED_CONFIG_DIR environment variable.
#
# Precedence:
#    --expected-dir (CLI argument) > EXPECTED_CONFIG_DIR (env var) > script default
#
# Summary Table:
# | How to Specify         | Example Command                                              | Takes Precedence? |
# |------------------------|-------------------------------------------------------------|-------------------|
# | Command-line argument  | python3 main.py --env testing --expected-dir test_config1   | Highest           |
# | Environment variable   | export EXPECTED_CONFIG_DIR=test_config1 && python3 main.py  | Lower             |
#
# Only one method is needed. If both are provided, the CLI argument is used.
# ===============================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--env', type=str, help='Environment to use (overrides COMMON_DICT_ENV)')
    parser.add_argument('--expected-dir', type=str, help='Expected config directory for verification (overrides EXPECTED_CONFIG_DIR)')
    args, unknown = parser.parse_known_args()
    try:
        run_config_reconciliation(domain='chemistry', environment=args.env, expected_dir=args.expected_dir)
    except Exception as e:
        print(f"[ERROR] Config reconciliation failed: {e}")
        sys.exit(1)
    env_loader = EnvironmentLoader(domain='chemistry', environment=args.env)
    orchestrator = ExtractorOrchestrator(domain='chemistry', environment=args.env)
    orchestrator.run()

if __name__ == '__main__':
    main() 