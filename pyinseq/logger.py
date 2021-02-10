#!/usr/bin/env python3
"""

Custom logger for pyinseq

"""
import sys
import tqdm
import logging

# Custom formatter
FORMATTER = logging.Formatter(
    fmt="%(asctime)s - %(levelname)s - %(module)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M",
)

def setup_filehandler(logger, logfile):
    fh = logging.FileHandler(logfile, 'a')
    fh.setFormatter(FORMATTER)
    logger.addHandler(fh)
    return

## ONLY LOGGER IN PYINSEQ ##
logger = logging.getLogger("pyinseq")
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler(stream=sys.stdout)
stream_handler.setFormatter(FORMATTER)
stream_handler.setLevel(logging.INFO)
logger.addHandler(stream_handler)

class TqdmBarLogger:
    """

    Custom progress bar that works with Logger

    """
    def __init__(self, logger: logging.Logger):
        self.logger = logger

    def info(self, msg):
        """ Removes logger from display and logs message"""
        self.p_bar.close()
        self.logger.info(msg)
        self.p_bar.disable = False
        # Displays it
        self.p_bar.refresh()

    def update(self):
        self.p_bar.update()

    def set_progress_bar(self, num_iterations: int, desc: str, unit: str):
        self.p_bar = tqdm.tqdm(total=num_iterations, desc=desc, unit=unit, leave=False, position=0)

    def clear_progress_bar(self):
        self.p_bar = None

# Only tqdm_logger
tqdm_logger = TqdmBarLogger(logger)









