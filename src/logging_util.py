import logging
import os

def setup_logging_from_env(env_loader):
    log_settings = env_loader.get_logging_settings()
    loglevel = log_settings.get('level', 'DEBUG').upper()  # Default to DEBUG
    logformat = log_settings.get('format', '[%(levelname)s] %(name)s: %(message)s')
    logfile = log_settings.get('file', None)
    error_logfile = log_settings.get('error_log_path', None)

    # Remove any existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Create formatter
    formatter = logging.Formatter(logformat)

    # Console handler (set to global loglevel)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(getattr(logging, loglevel, logging.DEBUG))
    console_handler.setFormatter(formatter)

    handlers = [console_handler]

    # File handler (set to global loglevel)
    if logfile:
        os.makedirs(os.path.dirname(logfile), exist_ok=True)
        file_handler = logging.FileHandler(logfile)
        file_handler.setLevel(getattr(logging, loglevel, logging.DEBUG))
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)

    # Error file handler (optional, set to ERROR)
    if error_logfile and error_logfile != logfile:
        os.makedirs(os.path.dirname(error_logfile), exist_ok=True)
        error_handler = logging.FileHandler(error_logfile)
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(formatter)
        handlers.append(error_handler)

    logging.basicConfig(
        level=getattr(logging, loglevel, logging.DEBUG),
        handlers=handlers
    )

    # Per-module debug
    for module in ['pubchem', 'rdkit', 'reaxys', 'config', 'extractor', 'env', 'output_generator', 'base_source', 'chebi', 'test_rdkit_manual']:
        if env_loader.get_debug_flag(module):
            logging.getLogger(module).setLevel(logging.DEBUG) 