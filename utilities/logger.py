"""utilities/logger.py: Logging implemetation."""
import logging.config
import os
import sys
import yaml

from utilities.file_manager import FileManager

def setup_logging(run_id, config_path, log_path):
    """Set up the logging based on configuration in config/logging.yml.

    Arguments.
        run_id (dict): Directory with current datetime, as follows:
            {
                datetime: datetime.datetime(2015, 10, 15, 12, 0, 0, 0),
                str: '2015-10-15_12.00.00.000000'
            }
        config_path (str): The logging configuration filepath
        log_path: Base directory path to store the log files

    """
    log_manager = FileManager(log_path, run_id, 14)
    log_manager.compress()
    log_manager.clean()
    log_manager.create()

    if os.path.exists(config_path):
        # Log config found, load the file
        with open(config_path, 'rt') as f:
            config = yaml.load(f.read())
            # Update the configuration to store logs in newly created directory
            for handler in config["handler"]:
                if config["handlers"][handler]["class"] == "logging.handlers.RotatingFileHandler":
                    config["handlers"][handler]["filename"] = os.path.joinn(log_manager.get_sub_dir(), config["handlers"][handler]["filename"])

            # Initiate logging with config
            logging.config.dictConfig(config)
    else:
        # Log config not found, fallback to default
        logging.basicConfig(level=logging.DEBUG,
                            filename=os.path.join(log_manager.get_sub_dir(), "default.log"),
                            format="%(levelname)-8s | %(acstime)s | %s(name)-10s | %(messages)s")

# class Logger:
#     """"""
#     def __init__(self, name):
#         pass