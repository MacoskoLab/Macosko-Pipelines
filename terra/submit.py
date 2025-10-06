import sys
import math
import gspread
import argparse
import pandas as pd
import firecloud.api as fapi
from google.auth import default
from google.cloud import storage
from gspread_dataframe import get_as_dataframe

wnamespace = "testmybroad"
workspace = "Macosko-Pipelines"
cnamespace = "macosko-pipelines"

def get_args():
    parser = argparse.ArgumentParser(description='Terra job submission script')
    parser.add_argument("workflow", type=str)
    parser.add_argument("bcl", type=str)
    parser.add_argument("index", type=str)
    parser.add_argument("--dryrun", action='store_true')
    args = parser.parse_args()
    return args

args = get_args()
workflow = args.workflow.lower()     ; print(f"workflow: {workflow}")
bcl = args.bcl.strip("/ \t\n\r")     ; print(f"     bcl: {bcl}")
index = args.index.strip("/ \t\n\r") ; print(f"   index: {index}")
dryrun = args.dryrun                 ; print(f"  dryrun: {dryrun}")

assert workflow in ["cellranger-count", "slide-tags", "recon", "reconstruction"]
assert not any(c.isspace() for c in bcl), f"remove whitespace from bcl ({bcl})"
assert not any(c.isspace() for c in index), f"remove whitespace from index ({index})"

# Load bucket
BUCKET = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
bucket = storage.Client().bucket(BUCKET)
bucket.reload()


# Load the worksheet, select the columns
sh = gspread.authorize(default()[0]).open_by_key("1NOaWXARQiSA6fquOtcouQPREPN4buYIf13tq_F6D9As")
if workflow in ["cellranger-count", "slide-tags"]:
    df = get_as_dataframe(sh.worksheet("Slide-tags"))
    cols = ["BCL", "Reference", "RNAIndex", "SBIndex", "Puck", "params"]
elif workflow in ["recon", "reconstruction"]:
    df = get_as_dataframe(sh.worksheet("Recon"))
    cols = ["BCL", "Index", "bc1", "bc2", "params"]


# Clean the worksheet
assert all(col in df.columns for col in cols)
df = df[cols].apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))

print(f"Total rows found: {len(df.index)}")


# Subset the worsheet to the BCL
df = df[df["BCL"] == bcl]
print(f"BCL rows found: {len(df.index)}")
assert len(df.index) > 0, f"BCL rows not found ({bcl})"


# Subset the worksheet to the index
if workflow in ["cellranger-count", "slide-tags"]:
    idx_col = "RNAIndex"
elif workflow in ["recon", "reconstruction"]:
    idx_col = "Index"
    
def assert_unique_column(series):
    assert series.notna().all(), f"Column has NA values:\n{series}"
    assert (~series.str.strip().eq("")).all(), f"Column has empty values:\n{series}"
    assert series.is_unique, f"Column has repeated values:\n{series}"

assert_unique_column(df[idx_col])
df = df if index.lower() == "all" else df[df[idx_col] == index]
print(f"Index rows found: {len(df.index)}")
assert len(df.index) >= 1, f"No index rows found ({index})"


# Assert necessary supplementary files exist
def assert_full_column(series):
    assert series.notna().all(), "Column has NA values"
    assert (~series.str.strip().eq("")).all(), "Column has empty values"

if workflow == "cellranger-count":
    # Assert the input references exist
    assert_full_column(df["Reference"])
    ref_blobs = bucket.list_blobs(prefix=f"references")
    refs = {blob.name.split("/")[1] for blob in ref_blobs if blob.name.endswith("reference.json")}
    assert df["Reference"].isin(refs).all(), f"Reference {set(df['Reference'])-refs} does not exist in the bucket"

    # Assert the output gene-expression folder does not exist
    count_blobs = bucket.list_blobs(prefix=f"gene-expression/{bcl}")
    counts = {blob.name.split("/")[2] for blob in count_blobs}
    assert not df["RNAIndex"].isin(counts).any(), f"Output {set(df['RNAIndex'])&counts} already exists in the bucket"

elif workflow == "slide-tags":
    # Assert RNA input exists
    count_blobs = bucket.list_blobs(prefix=f"gene-expression/{bcl}")
    counts = {blob.name.split("/")[2] for blob in count_blobs}
    assert df["RNAIndex"].isin(counts).all(), f"GEX for {set(df['RNAIndex'])-counts} does not exist in the bucket"

    # Assert puck file exists and is unique
    assert_full_column(df["Puck"])
    recon_pucks = [blob.name for blob in bucket.list_blobs(prefix=f"recon") if blob.name.endswith("/Puck.csv")]
    insitu_pucks = [blob.name for blob in bucket.list_blobs(prefix=f"pucks") if blob.name.endswith(".csv")]
    all_pucks = recon_pucks + insitu_pucks
    
    def pucks_to_URIs(pucks):
        URIs = []
        for puck in [s.strip() for s in pucks.split(',')]:
            matches = [s for s in all_pucks if puck in s]
            assert len(matches) == 1, f"{len(matches)} pucks found for '{puck}': {matches}"
            URIs.append("gs://"+BUCKET+"/"+matches[0])
        return URIs

    pucks = df["Puck"].apply(pucks_to_URIs)


# Compute memory requirements
fastq_blobs = bucket.list_blobs(prefix=f"fastqs/{bcl}")
fastqs = [(blob.name, blob.size) for blob in fastq_blobs if blob.name.endswith(".fastq.gz")]
getfastqsizes = lambda inds: [math.ceil(sum(s for n,s in fastqs if "/"+i+"_S" in n) / 1e9) for i in inds]

if workflow == "cellranger-count":
    # Compute the FASTQ sizes
    mem_GBs = getfastqsizes(df["RNAIndex"])
    mem_GBs = [math.ceil(5*mem+20) for mem in mem_GBs]

elif workflow == "slide-tags":
    # Compute the FASTQ sizes
    mem_GBs_fastq = getfastqsizes(df["SBIndex"])
    mem_GBs_fastq = [math.ceil(2*mem) for mem in mem_GBs_fastq]

    # Compute the SBcounts.h5 size
    tags_blobs = bucket.list_blobs(prefix=f"slide-tags/{bcl}")
    tags = [(blob.name, blob.size) for blob in tags_blobs if blob.name.endswith("/SBcounts.h5")]
    mem_GBs_mat = [max((s for n,s in tags if "/"+i+"/" in n), default=0) / 1e9 for i in df["RNAIndex"]]
    mem_GBs_mat = [math.ceil(20*mem) for mem in mem_GBs_mat]
    
    # Take the max (TODO)
    mem_GBs = [max(x,y) for x,y in zip(mem_GBs_fastq, mem_GBs_mat)]
    
elif workflow in ["recon", "reconstruction"]:
    # Compte the FASTQ sizes
    mem_GBs_fastq = getfastqsizes(df["Index"])
    mem_GBs_fastq = [math.ceil(2*mem) for mem in mem_GBs_fastq]
    print('mem_GBs_fastq', mem_GBs_fastq)
    
    # Compute the intermediate matrix sizes
    recon_blobs = bucket.list_blobs(prefix=f"recon/{bcl}")
    recons = [(blob.name, blob.size) for blob in recon_blobs if blob.name.endswith("/knn2.npz")]
    mem_GBs_mat = [max((s for n,s in recons if "/"+i+"/" in n), default=0) / 1e9 for i in df["Index"]]
    mem_GBs_mat = [math.ceil(25*mem) for mem in mem_GBs_mat]

    # Take the max (TODO)
    mem_GBs = [max(x,y) for x,y in zip(mem_GBs_fastq, mem_GBs_mat)]

assert all(m > 0 for m in mem_GBs), f"Incomplete memory estimation: {mem_GBs} (missing input files)"
mem_GBs = [math.ceil(max(mem, 64)) for mem in mem_GBs]
print(f"Memory (GB): {mem_GBs}")


# Compute the names, assert the jobs are not already running
job_names = ["_".join([workflow, idx, bcl]) for idx in df[idx_col]]
subs = fapi.list_submissions(wnamespace, workspace).json()
subs = [sub for sub in subs if sub["status"] not in ["Done","Aborted"]]
running_job_names = [sub["userComment"] for sub in subs if "userComment" in sub]
assert not set(job_names) & set(running_job_names), "Jobs already running!"


# Print + exit if testing
if dryrun:
    print(df)
    print("Dry run complete, no errors found - exiting...")
    sys.exit(0)
# TODO: recon won't error if missing fastqs when reselecting barcodes
# TODO: cached bc1/bc2 will likely differ from the sheet
# TODO: check data types

### Terra Submission ###########################################################

def submit(config, user_comment=""):
    # Validate the configuration
    res = fapi.validate_config(wnamespace, workspace, cnamespace, config).json()
    assert res["extraInputs"] == [], f"ERROR: extra input: \n{res['extraInputs']}"
    assert res["invalidInputs"] == {}, f"ERROR: invalid input: \n{res['invalidInputs']}"
    assert res["invalidOutputs"] == {}, f"ERROR: invalid output: \n{res['invalidOutputs']}"
    assert res["missingInputs"] == [], f"ERROR: missing input: \n{res['missingInputs']}"
    
    # Submit the job
    res = fapi.create_submission(wnamespace, workspace, cnamespace, config, user_comment=user_comment, use_callcache=False)
    if res.status_code != 201:
        print("Full error response:", res.json())
    assert res.status_code == 201, str(res.status_code) + ": " + res.json()['message']
    print(f"Submitted {config} {user_comment}")

def run_cellranger_count(bcl, index, reference, mem_GB, disk_GB, params=None, user_comment=""):
    body = fapi.get_workspace_config(wnamespace, workspace, cnamespace, "cellranger-count").json()
    body["inputs"]["cellranger_count.bcl"] = f'"{bcl}"'
    body["inputs"]["cellranger_count.index"] = f'"{index}"'
    body["inputs"]["cellranger_count.reference"] = f'"{reference}"'
    body["inputs"]["cellranger_count.mem_GB"] = f'{mem_GB}'
    body["inputs"]["cellranger_count.disk_GB"] = f'{disk_GB}'
    body["inputs"]["cellranger_count.params"] = f'"{params}"' if pd.notna(params) else f''
    body["inputs"]["cellranger_count.docker"] = f''
    res = fapi.update_workspace_config(wnamespace, workspace, cnamespace, "cellranger-count", body)
    assert res.status_code == 200, res.json()['message']
    
    submit("cellranger-count", user_comment)
    return True

def run_reconstruction(bcl, index, mem_GB, disk_GB, bc1=None, bc2=None, lanes=None, params=None, user_comment=""):
    body = fapi.get_workspace_config(wnamespace, workspace, cnamespace, "reconstruction").json()  
    body["inputs"]["reconstruction.bcl"] = f'"{bcl}"'
    body["inputs"]["reconstruction.index"] = f'"{index}"'
    body["inputs"]["reconstruction.mem_GB"] = f'{mem_GB}'
    body["inputs"]["reconstruction.disk_GB"] = f'{disk_GB}'
    body["inputs"]["reconstruction.bc1"] = f'{bc1}' if pd.notna(bc1) else f''
    body["inputs"]["reconstruction.bc2"] = f'{bc2}' if pd.notna(bc2) else f''
    body["inputs"]["reconstruction.lanes"] = f'{lanes}' if pd.notna(lanes) else f''
    body["inputs"]["reconstruction.params"] = f'"{params}"' if pd.notna(params) else f''
    body["inputs"]["reconstruction.docker"] = f''
    res = fapi.update_workspace_config(wnamespace, workspace, cnamespace, "reconstruction", body)
    assert res.status_code == 200, res.json()['message']
    
    submit("reconstruction", user_comment)
    return True

def run_slidetags(bcl, rna_index, sb_index, puck_paths, mem_GB, disk_GB, params=None, user_comment=""):
    body = fapi.get_workspace_config(wnamespace, workspace, cnamespace, "slide-tags").json()
    body["inputs"]["slide_tags.bcl"] = f'"{bcl}"'
    body["inputs"]["slide_tags.rna_index"] = f'"{rna_index}"'
    body["inputs"]["slide_tags.sb_index"] = f'"{sb_index}"'
    body["inputs"]["slide_tags.puck_paths"] = "[" + ", ".join(f'"{gs}"' for gs in puck_paths) + "]"
    body["inputs"]["slide_tags.mem_GB"] = f'{mem_GB}'
    body["inputs"]["slide_tags.disk_GB"] = f'{disk_GB}'
    body["inputs"]["slide_tags.params"] = f'"{params}"' if pd.notna(params) else f''
    body["inputs"]["slide_tags.docker"] = f''
    res = fapi.update_workspace_config(wnamespace, workspace, cnamespace, "slide-tags", body)
    assert res.status_code == 200, res.json()['message']

    submit("slide-tags", user_comment)
    return True

if workflow == "cellranger-count":
    for r, m, j in zip(df.itertuples(index=False), mem_GBs, job_names):
        run_cellranger_count(r.BCL, r.RNAIndex, r.Reference, 100, m, None, j)
        
elif workflow == "slide-tags":
    for r, p, m, j in zip(df.itertuples(index=False), pucks, mem_GBs, job_names):
        run_slidetags(r.BCL, r.RNAIndex, r.SBIndex, p, m, m, r.params, j)
elif workflow in ["recon", "reconstruction"]:
    for r, m, j in zip(df.itertuples(index=False), mem_GBs, job_names):
        assert r.Index.count("-") <= 1
        idx, _, lanes = r.Index.partition("-") ; 
        run_reconstruction(r.BCL, idx, m, m, r.bc1, r.bc2, lanes or None, r.params, j)

### Terra Commands #############################################################

# List all submissions
subs = fapi.list_submissions(wnamespace, workspace).json()
subs = [sub for sub in subs if sub["status"] not in ["Done","Aborted","Aborting"]]
print(f"Currently running submissions: {len(subs)}")

# Abort all submissions
# ids = [sub["submissionId"] for sub in subs]
# [fapi.abort_submission(wnamespace, workspace, submission_id) for submission_id in ids]
