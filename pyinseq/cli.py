"""

Main entry-point for pyinseq

"""

import sys
import traceback

# Module imports
from pyinseq.logger import pyinseq_logger
from pyinseq.parsers import get_args
from pyinseq.settings import SETTINGS_CONSTRUCTORS
from pyinseq.utils import (
    write_config_file,
    execute,
    create_experiment_directories,
    pickle_object,
    get_config_dict,
)


def run_tests():
    """Run the pyinseq test suite that ships with the package."""
    try:
        from pyinseq.tests.test_utils import run_pytest
    except ImportError:
        print(
            "pytest is required for --test. Install it with: pip install 'pyinseq[test]'",
            file=sys.stderr,
        )
        return 1
    return run_pytest()


def main():
    # Setup up universal logger
    pyinseq_logger.setup_logger()
    logger = pyinseq_logger.logger

    logger.info("Process command line arguments")
    parser, args = get_args()
    # Print help if no args given
    if not sys.argv[1:]:
        parser.print_help()
        sys.exit(0)

    if args.test:
        sys.exit(run_tests())

    if args.config:
        logger.info(f"Config file provided as: {args.config}")
        config_file = args.config
    else:
        config_file = write_config_file(
            args, format=args.config_format, default=args.get_default_config
        )

    config_dict = get_config_dict(config_file)
    command = config_dict["command"]

    constructor = SETTINGS_CONSTRUCTORS[command]
    settings = constructor(config_file)
    create_experiment_directories(settings)

    if settings.command != "genomeprep":
        settings.dump_sample_dict_to_yml()

    # Pickle settings
    pickle_object(settings, settings.settings_pickle)
    logger.info(f"Workflow selected: {settings.command}")
    logger.info(f"***** BEGIN {settings.command.upper()} WORKFLOW *****")

    returncode = 1
    try:
        returncode = execute(settings.snakemake_cmd)
    except Exception:
        logger.exception("Exception while executing snakemake subprocess")
        traceback.print_exc()

    # Cleanup settings and combine logging output
    settings.clean_up()
    if returncode != 0:
        logger.error(
            f"***** {settings.command.upper()} WORKFLOW FAILED "
            f"(snakemake exit code {returncode}) *****"
        )
        sys.exit(returncode)
    logger.info(f"***** {settings.command.upper()} WORKFLOW FINISHED *****")


if __name__ == "__main__":
    main()
