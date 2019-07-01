#
# Module to call from command line.

import sys
import argparse
from utilities.add_tools import add_parent_tool, add_tool, add_subtool

parser = argparse.ArgumentParser(description='Initialize metadata for a tool, script, or workflow.')


def main():
    raise NotImplementedError

if __name__ == "__main__":
    args = parser.parse_args()
    main()
    sys.exit(0)