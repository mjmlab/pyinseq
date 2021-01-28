#!/usr/bin/env python3
"""
Holds classes for pyinseq logging.

"""
import tqdm
import logging


# This controls the stdout logging.
formatter = logging.Formatter(
    fmt="%(asctime)s - %(levelname)s - %(module)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M",
)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(module)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M",
)


class TqdmBarLogger:
    """

    Custom progress bar that works with Logger

    """

    def __init__(
        self, logger: logging.Logger, num_reads: int, description: str, pos: int = 0
    ):
        self.logger = logger
        self.p_bar = tqdm.tqdm(
            total=num_reads, desc=description, unit="reads", leave=False, position=pos
        )

    def info(self, msg):
        """ Removes logger from display and logs message"""
        self.p_bar.close()
        self.logger.info(msg)
        self.p_bar.disable = False
        # Displays it
        self.p_bar.refresh()

    def update(self):
        """ Calls update on progress bar """
        self.p_bar.update()
        return


# Create logger just once and use in 'pyinseq' module
logger = logging.getLogger("pyinseq")
