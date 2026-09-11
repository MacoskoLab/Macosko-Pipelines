import sys
import json
import math
import gspread
import argparse
import pandas as pd
import firecloud.api as fapi
from google.cloud import storage
from google.auth import default as google_auth_default
from google.auth import impersonated_credentials
from gspread_dataframe import get_as_dataframe

# Built-in profile defaults reproduce the Macosko-Pipelines workspace behavior.
# Pass --config <profile.json> to override these values (e.g. bican.json) and
# serve as a drop-in for another workspace.
DEFAULT_PROFILE = {
    "sa": "pipelines",
    "wnamespace": "testmybroad",
    "workspace": "Macosko-Pipelines",
    "bucket": "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36",
    "sheet_key": "1NOaWXARQiSA6fquOtcouQPREPN4buYIf13tq_F6D9As",
    "fastq_match": "_S",
    "singlecell_concat": True,
    "use_sbbcl": False,
    "interactive_puck": False,
    "subfolder_aware_cache": False,
    "methods": {
        "cellranger-count": {"namespace": "macosko-pipelines", "config": "cellranger-count", "extra_inputs": []},
        "reconstruction":   {"namespace": "macosko-pipelines", "config": "reconstruction",   "extra_inputs": ["selection"]},
        "slide-tags":       {"namespace": "macosko-pipelines", "config": "slide-tags",       "extra_inputs": []},
    },
}

def get_args():
    parser = argparse.ArgumentParser(description='Terra job submission script')
    parser.add_argument("workflow", type=str)
    parser.add_argument("bcl", type=str)
    parser.add_argument("index", type=str)
    parser.add_argument("--config", type=str, default=None, help="Path to JSON profile overriding the built-in Macosko defaults")
    parser.add_argument("--sa", type=str, default=None, help="Service account to impersonate (short name or full email) for Google Sheets access; overrides the profile 'sa'")
    parser.add_argument("--dryrun", action='store_true')
    parser.add_argument("--mem", type=int, default=None, help="Override memory (GB) for all jobs")
    parser.add_argument("--branch", type=str, default=None, help="Git ref to pull pipeline scripts from: branch/tag name or commit SHA (slide-tags/recon only; default: main)")
    parser.add_argument("--pr", type=int, default=None, help="Pull request number to pull pipeline scripts from; overrides --branch (slide-tags/recon only)")
    parser.add_argument("--subfolder", type=str, default=None, help="Output subfolder name override (slide-tags/recon only)")
    parser.add_argument("--selection", type=str, default=None, help="Submit curated bead selections instead of the full puck (recon only): a selection name, or 'all' for every selection found. Selections are created by tools/puck-select.py")
    parser.add_argument("--bucket", type=str, default=None, help="Override GCS bucket name (slide-tags only)")
    parser.add_argument("--tag", type=str, default=None, help="Override docker image tag (slide-tags only; default: latest)")
    args = parser.parse_args()
    return args

args = get_args()

# Load the profile (built-in defaults, overridden by --config)
profile = json.loads(json.dumps(DEFAULT_PROFILE))  # deep copy
if args.config is not None:
    with open(args.config) as f:
        override = json.load(f)
    methods_override = override.pop("methods", {})
    profile.update(override)
    for name, method in methods_override.items():
        profile["methods"].setdefault(name, {}).update(method)
if args.sa is not None:
    profile["sa"] = args.sa

workflow = args.workflow.lower()     ; print(f"workflow: {workflow}")
bcl = args.bcl.strip("/ \t\n\r")     ; print(f"     bcl: {bcl}")
index = args.index.strip("/ \t\n\r") ; print(f"   index: {index}")
dryrun = args.dryrun                 ; print(f"  dryrun: {dryrun}")
mem_override = args.mem              ; print(f"     mem: {mem_override}")
branch = args.branch                 ; print(f"  branch: {branch}")
pr = args.pr                         ; print(f"      pr: {pr}")
subfolder = args.subfolder           ; print(f"  subfolder: {subfolder}")
selection = args.selection           ; print(f"  selection: {selection}")
bucket_override = args.bucket if args.bucket is not None else profile["bucket"] ; print(f"  bucket: {bucket_override}")
tag = args.tag                       ; print(f"     tag: {tag}")
print(f"workspace: {profile['wnamespace']}/{profile['workspace']}")

assert workflow in ["cellranger-count", "slide-tags", "recon", "reconstruction"]
assert not any(c.isspace() for c in bcl), f"remove whitespace from bcl ({bcl})"
assert not any(c.isspace() for c in index), f"remove whitespace from index ({index})"
assert selection is None or workflow in ["recon", "reconstruction"], "--selection is recon-only"

# Resolve the output subfolder exactly as the WDLs do, so every cache lookup below targets
# the same path the workflow will actually stat and reuse.
if subfolder:
    sub = subfolder
elif pr is not None:
    sub = f"pr-{pr}"
elif branch is not None and branch != "main":
    sub = branch
else:
    sub = None

# Load bucket (uses standard ADC with cloud-platform scope only)
BUCKET = profile["bucket"]
bucket = storage.Client().bucket(BUCKET)
bucket.reload()


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
if workflow in ["cellranger-count", "slide-tags"]:
    if workflow == "cellranger-count" and profile["singlecell_concat"]:
        df = pd.concat([get_as_dataframe(sh.worksheet("Slide-tags")),
                        get_as_dataframe(sh.worksheet("SingleCell"))], ignore_index=True)
    else:
        df = get_as_dataframe(sh.worksheet("Slide-tags"))
    cols = ["BCL", "Reference", "RNAIndex", "SBIndex", "Puck", "params"]
    if profile["use_sbbcl"]:
        if "SBBCL" not in df.columns:
            df["SBBCL"] = pd.NA
        cols = ["BCL", "Reference", "RNAIndex", "SBIndex", "SBBCL", "Puck", "params"]
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


# Expand recon rows into one job per bead selection
if workflow in ["recon", "reconstruction"]:
    def recon_base(idx):
        """GCS prefix reconstruction.wdl writes an index's outputs to (no trailing slash)."""
        d = idx[:-len("-12345678")] if idx.endswith("-12345678") else idx
        return f"recon/{bcl}/{d}" + (f"/{sub}" if sub else "")

    df["selection"] = pd.NA
    if selection:
        # A selection is named by the folder holding the selection.json that
        # tools/puck-select.py uploads; reconstruction.wdl reads the same file.
        found = {}
        sel_blobs = [b.name for b in bucket.list_blobs(prefix=f"recon/{bcl}")
                     if b.name.endswith("/selection.json")]
        for i in df[idx_col]:
            prefix = recon_base(i) + "/"
            for name in sel_blobs:
                if not name.startswith(prefix):
                    continue
                sel = name[len(prefix):-len("/selection.json")]
                # exactly one level down: not the base dir itself, not a nested subfolder
                if sel and "/" not in sel and selection in ("all", sel):
                    found.setdefault(i, []).append(sel)

        missing = [i for i in df[idx_col] if i not in found]
        assert not missing, (f"No selection.json found for {missing} under recon/{bcl}"
                             + (f" named '{selection}'" if selection != "all" else "")
                             + " - run tools/puck-select.py first")
        df = pd.DataFrame([{**r, "selection": s} for r in df.to_dict("records")
                                                 for s in sorted(found[r[idx_col]])])
        print(f"Selection rows found: {len(df.index)} ({sorted(set(df['selection']))})")


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
            if not profile["interactive_puck"]:
                assert len(matches) == 1, f"{len(matches)} pucks found for '{puck}': {matches}"
                URIs.append("gs://"+BUCKET+"/"+matches[0])
            elif len(matches) == 0:
                assert False, f"0 pucks found for '{puck}'"
            elif len(matches) == 1:
                URIs.append("gs://"+BUCKET+"/"+matches[0])
            else:
                print(f"\nMultiple pucks found for '{puck}':")
                for i, m in enumerate(matches):
                    print(f"  [{i}] {m}")
                choice = input(f"Select puck(s) (comma-separated, or 'all'): ").strip()
                if choice.lower() == "all":
                    URIs.extend("gs://"+BUCKET+"/"+m for m in matches)
                else:
                    for idx in choice.split(","):
                        idx = int(idx.strip())
                        assert 0 <= idx < len(matches), f"Invalid choice: {idx}"
                        URIs.append("gs://"+BUCKET+"/"+matches[idx])
        return URIs

    pucks = df["Puck"].apply(pucks_to_URIs)


# Compute memory requirements
fastq_blobs = bucket.list_blobs(prefix=f"fastqs/{bcl}")
fastqs = [(blob.name, blob.size) for blob in fastq_blobs if blob.name.endswith(".fastq.gz")]
getfastqsizes = lambda inds: [math.ceil(sum(s for n,s in fastqs if "/"+i+profile["fastq_match"] in n) / 1e9) for i in inds]

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
    if profile["subfolder_aware_cache"]:
        tags_suffix = f"/{sub}/SBcounts.h5" if sub else "/SBcounts.h5"
        mem_GBs_mat = [max((s for n,s in tags if n == f"slide-tags/{bcl}/{i}{tags_suffix}"), default=0) / 1e9 for i in df["RNAIndex"]]
    else:
        mem_GBs_mat = [max((s for n,s in tags if "/"+i+"/" in n), default=0) / 1e9 for i in df["RNAIndex"]]
    mem_GBs_mat = [math.ceil(20*mem) for mem in mem_GBs_mat]

    # Take the max (TODO)
    mem_GBs = [max(x,y) for x,y in zip(mem_GBs_fastq, mem_GBs_mat)]

elif workflow in ["recon", "reconstruction"]:
    # Compte the FASTQ sizes
    mem_GBs_fastq = getfastqsizes(df["Index"])
    mem_GBs_fastq = [math.ceil(2*mem) for mem in mem_GBs_fastq]

    # Compute the intermediate file sizes. Resolve the exact paths the WDL will stat rather
    # than substring matching, since a selection's knn2.npz also contains "/<Index>/" and
    # would otherwise be maxed in against the base run's.
    sizes = {blob.name: blob.size for blob in bucket.list_blobs(prefix=f"recon/{bcl}")}
    mem_GBs_mat = []
    for i, s in zip(df["Index"], df["selection"]):
        base = recon_base(i)
        work = f"{base}/{s}" if pd.notna(s) else base
        if sizes.get(f"{work}/knn2.npz", 0):
            mem_GBs_mat.append(math.ceil(25 * sizes[f"{work}/knn2.npz"] / 1e9))
        else:
            # No knn2.npz yet (first run, or knn.py died). Size off the matrix instead: an
            # observed knn2.npz is ~1.8x its matrix.csv.gz, so 25x knn ~= 45x matrix. A
            # selection's own matrix may not exist yet either, and is strictly smaller than
            # the base one it is cut from, so falling back to the base over-provisions safely.
            mat = sizes.get(f"{work}/matrix.csv.gz", 0) or sizes.get(f"{base}/matrix.csv.gz", 0)
            mem_GBs_mat.append(math.ceil(45 * mat / 1e9))

    # Take the max (TODO) - except a selection never reads the FASTQs
    mem_GBs = [y if pd.notna(s) else max(x, y)
               for x, y, s in zip(mem_GBs_fastq, mem_GBs_mat, df["selection"])]

assert all(m > 0 for m in mem_GBs), f"Incomplete memory estimation: {mem_GBs} (missing input files)"
mem_GBs = [math.ceil(max(mem, 64)) for mem in mem_GBs]
if mem_override is not None:
    mem_GBs = [mem_override] * len(mem_GBs)
print(f"Memory (GB): {mem_GBs}")


# Warn about CLI flags not supported by the active profile's method
_wf = "reconstruction" if workflow == "recon" else workflow
supported_extra = set(profile["methods"].get(_wf, {}).get("extra_inputs", []))
for flag_name, flag_val in [("branch", branch), ("pr", pr), ("subfolder", subfolder),
                            ("selection", selection), ("bucket", args.bucket), ("tag", tag)]:
    if flag_val is not None and flag_name not in supported_extra:
        print(f"WARNING: --{flag_name} is not supported by profile method '{_wf}' - ignoring")


# Compute the names, assert the jobs are not already running
wnamespace = profile["wnamespace"]
workspace = profile["workspace"]
# The selection belongs in the name: --selection all submits several jobs for one index at
# once, and they would otherwise share a userComment and collide in the in-flight guard below
if workflow in ["recon", "reconstruction"]:
    job_names = ["_".join([workflow, idx] + ([sel] if pd.notna(sel) else []) + [bcl])
                 for idx, sel in zip(df[idx_col], df["selection"])]
else:
    job_names = ["_".join([workflow, idx, bcl]) for idx in df[idx_col]]
resp = fapi.list_submissions(wnamespace, workspace)
if resp.status_code != 200:
    print(f"Terra API error {resp.status_code}: {resp.text[:500]}")
    sys.exit(1)
subs = resp.json()
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

def submit(config, ns, user_comment=""):
    # Validate the configuration
    res = fapi.validate_config(wnamespace, workspace, ns, config).json()
    assert res["extraInputs"] == [], f"ERROR: extra input: \n{res['extraInputs']}"
    assert res["invalidInputs"] == {}, f"ERROR: invalid input: \n{res['invalidInputs']}"
    assert res["invalidOutputs"] == {}, f"ERROR: invalid output: \n{res['invalidOutputs']}"
    assert res["missingInputs"] == [], f"ERROR: missing input: \n{res['missingInputs']}"

    # Submit the job
    res = fapi.create_submission(wnamespace, workspace, ns, config, user_comment=user_comment, use_callcache=False)
    assert res.status_code == 201, f"{res.status_code}: {res.json()['message']}"
    print(f"Submitted {config} {user_comment}")

def write_extra_inputs(body, extra_inputs, values):
    # Write only the inputs the active profile's method declares as supported
    for name in extra_inputs:
        if name not in values:
            continue
        val = values[name]
        # numeric inputs (pr) are written bare, strings are quoted
        if name == "pr":
            body["inputs"][f"slide_tags.{name}"] = f'{val}' if pd.notna(val) else f''
        else:
            body["inputs"][f"slide_tags.{name}"] = f'"{val}"' if pd.notna(val) else f''

def run_cellranger_count(method, bcl, index, reference, mem_GB, disk_GB, params=None, user_comment=""):
    ns, config = method["namespace"], method["config"]
    resp = fapi.get_workspace_config(wnamespace, workspace, ns, config)
    assert resp.status_code == 200, f"get_workspace_config({config}) failed ({resp.status_code}): {resp.text[:500]}"
    body = resp.json()
    body["inputs"]["cellranger_count.bcl"] = f'"{bcl}"'
    body["inputs"]["cellranger_count.index"] = f'"{index}"'
    body["inputs"]["cellranger_count.reference"] = f'"{reference}"'
    body["inputs"]["cellranger_count.mem_GB"] = f'{mem_GB}'
    body["inputs"]["cellranger_count.disk_GB"] = f'{disk_GB}'
    body["inputs"]["cellranger_count.params"] = f'"{params}"' if pd.notna(params) else f''
    body["inputs"]["cellranger_count.docker"] = f''
    res = fapi.update_workspace_config(wnamespace, workspace, ns, config, body)
    assert res.status_code == 200, res.json()['message']

    submit(config, ns, user_comment)
    return True

def run_reconstruction(method, bcl, index, mem_GB, disk_GB, bc1=None, bc2=None, lanes=None, params=None,
                       branch=None, pr=None, subfolder=None, selection=None, user_comment=""):
    ns, config = method["namespace"], method["config"]
    resp = fapi.get_workspace_config(wnamespace, workspace, ns, config)
    assert resp.status_code == 200, f"get_workspace_config({config}) failed ({resp.status_code}): {resp.text[:500]}"
    body = resp.json()
    body["inputs"]["reconstruction.bcl"] = f'"{bcl}"'
    body["inputs"]["reconstruction.index"] = f'"{index}"'
    body["inputs"]["reconstruction.mem_GB"] = f'{mem_GB}'
    body["inputs"]["reconstruction.disk_GB"] = f'{disk_GB}'
    body["inputs"]["reconstruction.bc1"] = f'{bc1}' if pd.notna(bc1) else f''
    body["inputs"]["reconstruction.bc2"] = f'{bc2}' if pd.notna(bc2) else f''
    body["inputs"]["reconstruction.lanes"] = f'{lanes}' if pd.notna(lanes) else f''
    body["inputs"]["reconstruction.params"] = f'"{params}"' if pd.notna(params) else f''
    for name in method["extra_inputs"]:
        val = {"branch": branch, "pr": pr, "subfolder": subfolder, "selection": selection}.get(name)
        if name == "pr":
            body["inputs"][f"reconstruction.{name}"] = f'{val}' if pd.notna(val) else f''
        else:
            body["inputs"][f"reconstruction.{name}"] = f'"{val}"' if pd.notna(val) else f''
    body["inputs"]["reconstruction.docker"] = f''
    res = fapi.update_workspace_config(wnamespace, workspace, ns, config, body)
    assert res.status_code == 200, res.json()['message']

    submit(config, ns, user_comment)
    return True

def run_slidetags(method, bcl, rna_index, sb_index, puck_paths, mem_GB, disk_GB, sb_bcl=None, params=None,
                  branch=None, pr=None, subfolder=None, bucket=None, tag=None, user_comment=""):
    ns, config = method["namespace"], method["config"]
    resp = fapi.get_workspace_config(wnamespace, workspace, ns, config)
    assert resp.status_code == 200, f"get_workspace_config({config}) failed ({resp.status_code}): {resp.text[:500]}"
    body = resp.json()
    body["inputs"]["slide_tags.bcl"] = f'"{bcl}"'
    body["inputs"]["slide_tags.rna_index"] = f'"{rna_index}"'
    body["inputs"]["slide_tags.sb_index"] = f'"{sb_index}"'
    body["inputs"]["slide_tags.puck_paths"] = "[" + ", ".join(f'"{gs}"' for gs in puck_paths) + "]"
    body["inputs"]["slide_tags.mem_GB"] = f'{mem_GB}'
    body["inputs"]["slide_tags.disk_GB"] = f'{disk_GB}'
    body["inputs"]["slide_tags.params"] = f'"{params}"' if pd.notna(params) else f''
    write_extra_inputs(body, method["extra_inputs"],
                       {"sb_bcl": sb_bcl, "branch": branch, "pr": pr,
                        "subfolder": subfolder, "bucket": bucket, "tag": tag})
    body["inputs"]["slide_tags.docker"] = f''
    res = fapi.update_workspace_config(wnamespace, workspace, ns, config, body)
    assert res.status_code == 200, res.json()['message']

    submit(config, ns, user_comment)
    return True

if workflow == "cellranger-count":
    method = profile["methods"]["cellranger-count"]
    for r, m, j in zip(df.itertuples(index=False), mem_GBs, job_names):
        run_cellranger_count(method, r.BCL, r.RNAIndex, r.Reference, 100, m, None, j)

elif workflow == "slide-tags":
    method = profile["methods"]["slide-tags"]
    for r, p, m, j in zip(df.itertuples(index=False), pucks, mem_GBs, job_names):
        sb_bcl = getattr(r, "SBBCL", None) if profile["use_sbbcl"] else None
        run_slidetags(method, r.BCL, r.RNAIndex, r.SBIndex, p, m, m, sb_bcl, r.params,
                      branch, pr, subfolder, bucket_override, tag, j)

elif workflow in ["recon", "reconstruction"]:
    method = profile["methods"]["reconstruction"]
    for r, m, j in zip(df.itertuples(index=False), mem_GBs, job_names):
        assert r.Index.count("-") <= 1
        idx, _, lanes = r.Index.partition("-") ;
        run_reconstruction(method, r.BCL, idx, m, m, r.bc1, r.bc2, lanes or None, r.params,
                           branch, pr, subfolder, r.selection, j)

### Terra Commands #############################################################

# List all submissions
subs = fapi.list_submissions(wnamespace, workspace).json()
subs = [sub for sub in subs if sub["status"] not in ["Done","Aborted","Aborting"]]
print(f"Currently running submissions: {len(subs)}")

# Abort all submissions
# ids = [sub["submissionId"] for sub in subs]
# [fapi.abort_submission(wnamespace, workspace, submission_id) for submission_id in ids]
