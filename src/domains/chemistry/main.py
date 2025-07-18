import argparse
import sys
import os
from config.env_loader import EnvironmentLoader
from logging_util import setup_logging_from_env

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--env', type=str, help='Environment to use (overrides COMMON_DICT_ENV)')
    parser.add_argument('--expected-dir', type=str, help='Expected config directory for verification (overrides EXPECTED_CONFIG_DIR)')
    args, unknown = parser.parse_known_args()

    # Set up environment and logging early
    env_loader = EnvironmentLoader(domain='chemistry', environment=args.env)
    setup_logging_from_env(env_loader)

    from config_reconciliation import run_config_reconciliation
    from domains.chemistry.extractor import ExtractorOrchestrator

    try:
        run_config_reconciliation(domain='chemistry', environment=args.env, expected_dir=args.expected_dir)
    except Exception as e:
        import logging
        logging.getLogger("main").error(f"[ERROR] Config reconciliation failed: {e}")
        sys.exit(1)

    orchestrator = ExtractorOrchestrator(domain='chemistry', environment=args.env)
    orchestrator.run()

if __name__ == '__main__':
    main() 