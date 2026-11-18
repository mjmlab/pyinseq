#!/usr/bin/env python3

"""

Custom logger for pyinseq

"""

import io
import sys
import tqdm
import logging


class CustomLoggers:

    FORMATTER = logging.Formatter(
        fmt="{asctime} - {levelname} - {module} - {message}",
        datefmt="%Y-%m-%d %H:%M",
        style="{",
    )

    def __init__(self):
        # Set io's
        self.logger_io = io.StringIO()
        self.snake_io = io.StringIO()
        # Available at import time; handlers are attached by setup_logger()
        self.logger = logging.getLogger("pyinseq")

    def add_logfile_paths(self, settings):
        self.log = settings.log
        self.summary_log = settings.summary_log

    def summarize_step(self):
        """Dumps io buffer into a log summary file"""
        with open(self.summary_log, "a+") as f:
            snake_message = self.snake_io.getvalue()
            pyinseq_message = self.logger_io.getvalue()
            f.write(snake_message)
            f.write("(PYINSEQ INFO)\n" + pyinseq_message)
        # Flush io's
        self.snake_io.flush()
        self.logger_io.flush()

    def setup_logger(self, formatter=None):
        """Sets up logger with sys.stdout handler"""
        self.logger = logging.getLogger("pyinseq")
        self.logger.setLevel(logging.INFO)
        # Replace the stdout handler from any previous call so that repeated
        # setup (e.g. running the CLI several times in one process) does not
        # duplicate output and always writes to the current sys.stdout
        previous = getattr(self, "_stdout_handler", None)
        if previous is not None:
            self.logger.removeHandler(previous)
        stream_handler = logging.StreamHandler(stream=sys.stdout)
        if formatter:
            stream_handler.setFormatter(formatter)
        else:
            stream_handler.setFormatter(CustomLoggers.FORMATTER)
        stream_handler.setLevel(logging.INFO)
        self.logger.addHandler(stream_handler)
        self._stdout_handler = stream_handler
        return


## ONLY LOGGER IN PYINSEQ ##
pyinseq_logger = CustomLoggers()


class TqdmBarLogger:
    """

    Custom progress bar that works with Logger

    """

    def __init__(
        self, logger: logging.Logger, num_iterations: int, desc: str, unit: str
    ):
        self.logger = logger
        self.p_bar = tqdm.tqdm(
            total=num_iterations,
            desc=desc,
            unit=unit,
            leave=False,
            position=0,
            colour="green",
        )

    def info(self, msg):
        """Removes logger from display and logs message"""
        self.p_bar.close()
        self.logger.info(msg)
        # Set disable to False, brings bar back
        self.p_bar.disable = False

    def update(self):
        self.p_bar.update()

    def close(self):
        self.p_bar = None


def add_fileHandler(logger, logfile, formatter=None):
    """Adds a filehandler to logger"""
    fh = logging.FileHandler(logfile, "a")
    if formatter:
        fh.setFormatter(formatter)
    else:
        fh.setFormatter(CustomLoggers.FORMATTER)
    logger.addHandler(fh)
    return


def add_streamHandler(logger, stream, formatter=False):
    """Adds a stream handler to logger"""
    stream_handler = logging.StreamHandler(stream=stream)
    if formatter:
        stream_handler.setFormatter(CustomLoggers.FORMATTER)
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)
    return
