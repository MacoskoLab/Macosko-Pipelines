import os
import pandas as pd
from functools import reduce
from collections import Counter
from google.cloud import storage
BUCKET = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
bucket = storage.Client().bucket(BUCKET)
bucket.reload()

# Load the blobs
fastqs = [blob.name for blob in bucket.list_blobs(prefix=f"fastqs") if blob.name.endswith(".fastq.gz")]
gex = [blob.name for blob in bucket.list_blobs(prefix=f"gene-expression") if blob.name.endswith("/web_summary.html")]
recon = [blob.name for blob in bucket.list_blobs(prefix=f"recon") if blob.name.endswith("/summary.pdf")]
tags = [blob.name for blob in bucket.list_blobs(prefix=f"slide-tags") if blob.name.endswith("/obj.qs")]

# Split on BCL/Index
df1 = pd.DataFrame([(k[0],k[1],v) for k,v in Counter([(a.split("/")[1], a.split("/")[2].split("_S")[0]) for a in fastqs]).items()], columns=["BCL", "Index", "FASTQs"])
df2 = pd.DataFrame([(g.split("/")[1], g.split("/")[2], g) for g in gex], columns=["BCL", "Index", "GEX"])
df3 = pd.DataFrame([(t.split("/")[1], t.split("/")[2], t) for t in tags], columns=["BCL", "Index", "Tags"])
df4 = pd.DataFrame([(r.split("/")[1], r.split("/")[2], r) for r in recon], columns=["BCL", "Index", "Recon"])

# Assert only 1 file exists
assert not df1.duplicated(subset=["BCL", "Index"]).any()
assert not df2.duplicated(subset=["BCL", "Index"]).any()
assert not df3.duplicated(subset=["BCL", "Index"]).any()
assert not df4.duplicated(subset=["BCL", "Index"]).any()

# Merge into df
df = df1.merge(df2, on=["BCL","Index"]).merge(df3, on=["BCL","Index"]).merge(df4, on=["BCL", "Index"])
df = reduce(lambda left, right: pd.merge(left, right, on=["BCL", "Index"], how="outer"), [df1, df2, df3, df4])
df = df.sort_values(by=['BCL', 'Index'])

# Format links
def blob_to_link(blob, name):
    if isinstance(blob, str):
        return f'=HYPERLINK("https://storage.cloud.google.com/{BUCKET}/{blob}", "{name}")'
    else:
        return pd.NA

df["GEX"] = df["GEX"].apply(lambda g: blob_to_link(g, "web_summary.html"))
df["Recon"] = df["Recon"].apply(lambda r: blob_to_link(r, "summary.pdf"))

################################################################################

import gspread
from google.auth import default
from gspread_dataframe import set_with_dataframe

# Upload sheet
gc = gspread.authorize(default()[0])
sh = gc.open_by_url("https://docs.google.com/spreadsheets/d/1NOaWXARQiSA6fquOtcouQPREPN4buYIf13tq_F6D9As")
ws = sh.worksheet("Results")
set_with_dataframe(ws, df)

# Add horizontal lines
requests = []
for i in range(1, len(df)):
    if df.loc[i, 'BCL'] != df.loc[i-1, 'BCL']:
        requests.append({
            "updateBorders": {
                "range": {"sheetId": ws.id, "startRowIndex": i, "endRowIndex": i+1},
                "bottom": {"style": "SOLID"}
            }
        })

# Auto-resize all columns
requests.append({
    "autoResizeDimensions": {
        "dimensions": {"sheetId": ws.id, "dimension": "COLUMNS", "startIndex": 0, "endIndex": ws.col_count}
    }
})

# Submit changes
ws.spreadsheet.batch_update({"requests": requests})
