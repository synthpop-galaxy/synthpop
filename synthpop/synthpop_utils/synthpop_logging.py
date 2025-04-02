"""
This module consist the logging class of the SynthPopFramework.
It mainly works like a standard python logger.
But can change the logging location and provided function
to create sections in the logfile .
"""

__all__ = ["SynthpopLogger", "logger", "log_basic_statistics"]
__author__ = "J. Klüter"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2023-02-10"

import logging
import shutil
from datetime import datetime
import os
import tempfile
import numpy as np

try:
    from constants import SYNTHPOP_DIR
except (ImportError, ValueError):
    from ..constants import SYNTHPOP_DIR

LENGTH = 75


class SynthpopLogger(logging.Logger):
    def __init__(
            self,
            name, level=logging.INFO,
            debug_file=os.path.join(SYNTHPOP_DIR, "synthpop.logfile.log")
            ):
        super().__init__(name, level)

        self.level = logging.DEBUG
        self.stream_level = level - 5
        self.file_level = level - 5

        # Create a stream handler with default INFO level
        self.stream_logger = logging.StreamHandler()
        self.stream_logger.setLevel(level)

        # create debug handler
        self.debugger = logging.FileHandler(debug_file)
        self.debugger.setLevel(logging.DEBUG)

        # Create the logging formatter
        self.file_formatter = logging.Formatter('%(relativeCreated)d - %(message)s')
        self.stream_formatter = logging.Formatter(' %(relativeCreated)d - %(message)s')
        self.debug_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

        # Add the stream handler to the logger object
        self.stream_logger.setFormatter(self.stream_formatter)
        self.addHandler(self.stream_logger)
        # add debugger to
        self.debugger.setFormatter(self.debug_formatter)
        self.addHandler(self.debugger)

        # Set a flag to indicate whether file logging is enabled or not
        self.file_logging_enabled = False
        self.temp_file = None
        self.current_file = None
        self.filelogger = None

    def setup_file_logging(self, stream_level, file_level=None):
        """
        sets up the file logging.

        Parameters
        ----------
        stream_level : int
            set level for stream
        file_level: int or  None
            set level for stream
            if None it will be stream_level - 5

        Returns
        -------

        """
        # create a temporarily log file.
        self.file_level = stream_level - 5 if file_level is None else file_level
        self.stream_level = stream_level
        self.current_file = tempfile.NamedTemporaryFile(mode='w+', delete=True)

        # Set the logging levels for the stream handlers
        self.stream_logger.setLevel(stream_level)

        # Create file handler.
        self.filelogger = logging.FileHandler(self.current_file.name)
        self.filelogger.setLevel(self.file_level)

        # Add the handler to the logger object
        self.filelogger.setFormatter(self.file_formatter)
        self.debugger.setFormatter(self.debug_formatter)

        self.addHandler(self.filelogger)

        # Set the flag to indicate that file logging is now enabled
        self.file_logging_enabled = True

        # log date and time
        now = datetime.now()
        self.info(f'Execution Date: {now.strftime("%d-%m-%Y %H:%M:%S")}')

    def create_info_section(self, msg):
        if len(msg) > LENGTH - 6:
            updated_msg = f"## {msg} ##"
        else:
            ll = LENGTH - len(msg) - 2
            updated_msg = f"{'#' * (ll // 2)} {msg} {'#' * ((ll + 1) // 2)}"

        if self.file_logging_enabled:
            self.filelogger.stream.write(f"\n\n{updated_msg}\n")
        if self.stream_level <= 25:
            self.stream_logger.stream.write(f"\n\n{updated_msg}\n")

    def create_info_subsection(self, msg, level=25):

        if len(msg) > LENGTH - 2:
            updated_msg = f"# {msg} --"
        else:
            ll = LENGTH - len(msg) - 3
            updated_msg = f"# {msg} {'-' * ll}"

        if self.file_logging_enabled:
            self.filelogger.stream.write(f"\n\n{updated_msg}\n")
        if self.stream_level <= level:
            self.stream_logger.stream.write(f"\n\n{updated_msg}\n")

    def save_log_file(self, file_path):
        self.current_file.seek(0)
        with open(file_path, 'w') as f:
            shutil.copyfileobj(self.current_file, f)

    def update_location(self, filename, no_log_file=False):
        """
        Copy the logfile header to a new location
        and continues logging at the new location

        Parameters
        ----------
        filename : str
            new logging location
        no_log_file : bool
            if True. replace the logging with a tmp file to prevent file logging
        """
        if no_log_file:
            if self.temp_file is not None:
                self.temp_file.close()
            self.temp_file = tempfile.NamedTemporaryFile(mode='w', delete=True)
            filename = self.temp_file.name
        else:
            dirname= os.path.dirname(filename)
            os.makedirs(dirname, exist_ok=True)
            self.save_log_file(filename)
        if self.filelogger in self.handlers:
            self.removeHandler(self.filelogger)
        self.filelogger = logging.FileHandler(filename)
        # set up new file handler with new location

        self.filelogger.setLevel(self.stream_level)
        self.filelogger.setFormatter(self.file_formatter)
        self.addHandler(self.filelogger)

    def cleanup(self):
        if self.file_logging_enabled:
            self.current_file.close()

            if self.temp_file is not None:
                self.temp_file.close()
            if self.file_logging_enabled:
                self.remove_file_handler()
            self.file_logging_enabled = False

    def flush(self):
        for handler in self.handlers:
            handler.flush()

    def remove_file_handler(self):
        if self.filelogger:
            self.removeHandler(self.filelogger)
            self.filelogger.close()
            self.filelogger = None

    def log2file(self, level, msg, *args, **kwargs):
        """
        Log a message only to the file handler of the specified logger.

        Parameters
        ----------
        level : int
            The logging level of the message to log.
        msg : str
            The message to log.
        *args : positional arguments
            Additional positional arguments to include in the log message.
        **kwargs : keyword arguments
            Additional keyword arguments to include in the log message.

        Returns
        -------
        None
        """
        self.debug(f'Logging to file: {msg}')
        if self.file_logging_enabled:
            rec = logging.LogRecord(self.name, level, '', level, msg, args, None)
            self.filelogger.emit(rec)

    def log2stream(self, level, msg, *args, **kwargs):
        """
        Log a message only to the stream handler of the specified logger.

        Parameters
        ----------
        level : int
            The logging level of the message to log.
        msg : str
            The message to log.
        *args : positional arguments
            Additional positional arguments to include in the log message.
        **kwargs : keyword arguments
            Additional keyword arguments to include in the log message.

        Returns
        -------
        None
        """
        self.debug(f'Logging to file: {msg}')
        rec = logging.LogRecord(self.name, level, '', level, msg, args, None)
        self.stream_logger.emit(rec)

logger = SynthpopLogger('synthpop_logging')

def log_basic_statistics(df, var_name, criteria=None):
    logger.log(15, '# Basic Statistics:')
    if criteria is not None:
        df = df.loc[criteria]
    if len(df) == 0:
        # Check if there are any stars
        logger.log(15, f'# No Stars from population {var_name}:')
        return

    logger.log(15, f'{var_name} = [')
    for col in df.columns:
        msg = f'    {{"name":{col}, '
        msg += f'"mean": {df.loc[:, col].mean(skipna=True):.4f}, '
        msg += f'"min": {df.loc[:, col].min(skipna=True):.4f}, '
        msg += f'"max": {df.loc[:, col].max(skipna=True):.4f}, '
        msg += f'"std": {df.loc[:, col].std(skipna=True):.4f}}},'
        logger.log(15,  msg)
    logger.log(15,'    ]')


