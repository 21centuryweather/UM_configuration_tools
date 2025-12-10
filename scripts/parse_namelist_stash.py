"""
Parse UM stash namelist sections into CSV.

Usage:
    python parse_stash.py input.txt output.csv
"""

import re
import csv
import sys
from pathlib import Path


def parse_stash(input_file, output_file):

    with open(input_file, "r") as f:
        lines = [line.strip() for line in f]

    # Skip until first [namelist:umstash_streq or [!namelist:umstash_streq
    start_idx = next(
        (i for i, line in enumerate(lines) if re.match(r"^\[!?namelist:umstash_streq", line)),
        None,
    )
    if start_idx is None:
        raise ValueError("No namelist sections found in file.")

    lines = lines[start_idx:]

    # Identify section start and end indices
    section_starts = [i for i, line in enumerate(lines) if re.match(r"^\[!?namelist:umstash_streq", line)]
    section_ends = section_starts[1:] + [len(lines)]

    records = []

    for start, end in zip(section_starts, section_ends):
        block = lines[start:end]
        header = block[0]

        incl = not header.startswith("[!namelist")
        match = re.search(r"\(([^)]+)\)", header)
        stash_id = match.group(1) if match else None

        record = {"id": stash_id, "incl": incl}

        # Extract key=value pairs
        for line in block[1:]:
            if "=" in line:
                key, val = line.split("=", 1)
                key = key.strip()
                val = val.strip().strip("'")
                record[key] = val

        records.append(record)

    # Collect all unique keys to use as CSV columns
    all_keys = sorted({k for r in records for k in r.keys()})
    fieldnames = all_keys

    # Write CSV
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"✅ Parsed {len(records)} stash sections and saved to '{output_file}'.")


def main():
    if len(sys.argv) != 3:
        print("Usage: python parse_namelist_stash.py <input.txt> <output.csv>")
        sys.exit(1)

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    if not input_path.exists():
        print(f"❌ Input file not found: {input_path}")
        sys.exit(1)

    parse_stash(input_path, output_path)


if __name__ == "__main__":
    main()
