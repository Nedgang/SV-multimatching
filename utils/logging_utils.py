##########
# IMPORT #
##########
import logging
import sys


############
# FUNCTION #
############
def setup_logging(verbose: bool = False, log_file: str = None) -> None:
    """
    Configure the application logging system.

    Parameters
    ----------
    verbose :
        If True, display detailed log messages on stderr.
        Otherwise, only informational messages are displayed.
    log_file :
        If provided, log messages are also written to this file.
    """
    level = logging.DEBUG if verbose else logging.INFO
    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s", datefmt="%H:%M:%S"
    )
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logging.basicConfig(level=level, handlers=[console_handler], force=True)
    if log_file is not None:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logging.getLogger().addHandler(file_handler)

