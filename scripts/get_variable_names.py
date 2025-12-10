"""
Extend parsed UM stash CSV with skos:notation and merge with reference stash table.

Usage:
    python join_stash.py parsed.csv stash.csv output.csv
"""

import sys
import pandas as pd


def main():
    if len(sys.argv) != 4:
        print("Usage: python get_variable_names.py <parsed.csv> <stash.csv> <output.csv>")
        sys.exit(1)

    parsed_file, ref_file, output_file = sys.argv[1:]

    # Load parsed data (from parse_stash.py)
    df = pd.read_csv(parsed_file)

    # Convert isec/item to numeric and fill missing values
    df["isec"] = pd.to_numeric(df["isec"], errors="coerce").fillna(0).astype(int)
    df["item"] = pd.to_numeric(df["item"], errors="coerce").fillna(0).astype(int)

    # Create skos:notation column
    df["skos:notation"] = df.apply(
        lambda r: f"m01s{r['isec']:02d}i{r['item']:03d}", axis=1
    )

    # Load reference stash table
    ref = pd.read_csv(ref_file)

    # Merge on skos:notation
    merged = df.merge(ref[["skos:notation", "rdfs:label"]], on="skos:notation", how="left")
    merged = merged.rename(columns={"rdfs:label": "name"})
    merged = merged.rename(columns={"skos:notation": "stash_code"})


    # Write output
    merged.to_csv(output_file, index=False)
    print(f"✅ Merged file written to {output_file}")


if __name__ == "__main__":
    main()
