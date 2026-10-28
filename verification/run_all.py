"""Run the repository's finite scientific checks."""
import argparse
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from verify_core import run_checks as core_checks
from verify_wml import run_checks as wml_checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path.")
    args = parser.parse_args()
    report = {
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "core": core_checks(),
        "wml": wml_checks(),
        "status": "passed",
        "scope": "Finite identity and convention checks; general claims rely on the proofs.",
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
