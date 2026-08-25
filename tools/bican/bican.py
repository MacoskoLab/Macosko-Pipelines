"""
python tools/bican/bican.py <BCL> <RNAindex> uses bican pipeline log spreadsheet to identify the row and print the info in the following format
python tools/bican/bican.py <BCL> prints the info for all matching rows

Output is printed to stdout (with a "bican all < /dev/null" hint between rows,
for copy-pasting into a shell that has already `source bican.sh`-d). The same
info is also written to <BCL>.txt as plain data: one block of `KEY=VALUE`
lines per matching row, blocks separated by a blank line, no hint lines. That
file is meant to be replayed with bican.sh, e.g.:
  ./bican.sh <BCL>.txt all --dryrun

Auth: uses your ADC to impersonate a service account (no key file needed).
  Run once beforehand:  gcloud auth application-default login
  You need the "Service Account Token Creator" role on that SA.
Pass --config <profile.json> (e.g. bican.json) to override the "sa"/"sheet_key" defaults below.
"""

import sys
import json
import argparse
import gspread
import pandas as pd
from google.auth import default as google_auth_default
from google.auth import impersonated_credentials

BUCKET = "gs://fc-79019ce4-0d1c-4750-b109-458c7c5c4e68"

DEFAULT_PROFILE = {
    "sa": "pipelines@velina-208320.iam.gserviceaccount.com",
    "sheet_key": "1-xdrCPnTKmrnaPlew7ToOqjPRv0rsSUP5GOrAKiZK64",
}

# Parse arguments
parser = argparse.ArgumentParser(description='BICAN pipeline log lookup')
parser.add_argument("bcl", type=str)
parser.add_argument("rna_index", type=str, nargs="?", default=None)
parser.add_argument("--config", type=str, default=None, help="Path to JSON profile overriding the built-in BICAN defaults (e.g. bican.json)")
parser.add_argument("--sa", type=str, default=None, help="Service account to impersonate (short name or full email) for Google Sheets/Drive access; overrides the profile 'sa'")
args = parser.parse_args()

# Load the profile (built-in defaults, overridden by --config)
profile = dict(DEFAULT_PROFILE)
if args.config is not None:
    with open(args.config) as f:
        profile.update(json.load(f))  # unknown (submit/results-only) keys are harmless
if args.sa is not None:
    profile["sa"] = args.sa

bcl = args.bcl.strip("/ \t\n\r")
rna_index = args.rna_index.strip("/ \t\n\r") if args.rna_index else None

# Load the worksheet using impersonated service account credentials
# (avoids needing a service account key file - uses ADC + impersonation)
source_creds, gcp_project = google_auth_default(scopes=["https://www.googleapis.com/auth/cloud-platform"])
if "@" in profile["sa"]:
    sa_email = profile["sa"]
else:
    assert gcp_project, "Could not infer project from ADC; pass a full SA email via --sa or profile 'sa'"
    sa_email = f"{profile['sa']}@{gcp_project}.iam.gserviceaccount.com"
sheets_creds = impersonated_credentials.Credentials(
    source_credentials=source_creds,
    target_principal=sa_email,
    target_scopes=[
        "https://www.googleapis.com/auth/spreadsheets.readonly",
        "https://www.googleapis.com/auth/drive.readonly",
    ],
    lifetime=3600,
)
sh = gspread.authorize(sheets_creds).open_by_key(profile["sheet_key"])
ws = sh.worksheet("Slide-tags")
df = pd.DataFrame(ws.get_all_records(value_render_option=gspread.utils.ValueRenderOption.formatted))
df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))

# Subset to BCL
df = df[df["BCL"] == bcl]
assert len(df.index) > 0, f"No rows found for BCL '{bcl}'"

# Subset to RNAIndex if provided (unless "all")
if rna_index is not None and rna_index.lower() != "all":
    df = df[df["RNAIndex"] == rna_index]
    assert len(df.index) > 0, f"No rows found for BCL '{bcl}' and RNAIndex '{rna_index}'"

# Column name mapping: spreadsheet column -> output variable name
COL_MAP = {
    "Sample": "sample",
    "BCL": "bcl",
    "RNAIndex": "rnaidx",
    "Puck": "puck",
    "PuckAlt": "puckalt",
    "SBIndex": "spidx",
    "fastqpost": "fastqpost",
    "spfastqpath": "spfastqpath",
    "fastqprefix": "spfastqprefix",
    "gex": "gex",
    "gexpath": "gexpath",
}

def quote_if_needed(key, value):
    """Wrap value in quotes if it's empty or contains spaces/special chars."""
    if key in ("gexpath", "fastqpost", "spfastqpath"):
        return f'"{value}"'
    if " " in value or "," in value:
        return f'"{value}"'
    return value

def row_lines(row):
    lines = []
    for col, var in COL_MAP.items():
        value = row.get(col, "")
        value = "" if pd.isna(value) else str(value)
        value = value.replace("\n", "").replace("\r", "")
        if col == "Puck" and value == "":
            alt = row.get("PuckAlt", "")
            value = "" if pd.isna(alt) else str(alt).replace("\n", "").replace("\r", "")
        lines.append(f"{var}={quote_if_needed(var, value)}")
    lines.append(f"gcp_bucket={BUCKET}")
    return lines

blocks = [row_lines(row) for _, row in df.iterrows()]

# Print to stdout with copy-paste hints between rows
for i, lines in enumerate(blocks):
    if i > 0:
        print()
        print("bican all < /dev/null")
        print()
    print("\n".join(lines))
print()
print("bican all < /dev/null")
print()

# Write pure data (no hint lines) for bican.sh to replay
outfile = f"{bcl}.txt"
with open(outfile, "w") as out:
    out.write("\n\n".join("\n".join(lines) for lines in blocks) + "\n")
print(f"Wrote output to {outfile}", file=sys.stderr)
