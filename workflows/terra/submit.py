import sys
import argparse
import pandas as pd
import firecloud.api as fapi
from google.cloud import storage
BUCKET = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"

def get_args():
    parser = argparse.ArgumentParser(description='Terra job submission script')
    parser.add_argument("workflow", type=str)
    parser.add_argument("bcl", type=str)
    parser.add_argument("index", type=str)
    parser.add_argument("--dryrun", action='store_true')
    args = parser.parse_args()
    return args
args = get_args()
workflow = args.workflow             ; print(f"workflow: {workflow}")
bcl = args.bcl.strip("/ \t\n\r")     ; print(f"     bcl: {bcl}")
index = args.index.strip("/ \t\n\r") ; print(f"   index: {index}")
dryrun = args.dryrun                 ; print(f"  dryrun: {dryrun}")

assert not any(c.isspace() for c in bcl), f"remove whitespace from bcl ({bcl})"
assert not any(c.isspace() for c in index), f"remove whitespace from index ({index})"

# cache=TRUE

# GCP helpers:
bucket = storage.Client().bucket(BUCKET)
bucket.reload() # Check bucket access

# Trim path if entire URL is specified
# def trim_path(path):
#     if path.startswith("gs://"):
#         path = path[len("gs://"):]
#     if path.startswith(BUCKET):
#         path = path[len(BUCKET):]
#     if path.startswith("/"):
#         path = path[len("/"):]
#     return path

# Check if a file exists in the bucket
# def file_exists(path):
#     blob = bucket.blob(trim_path(path))
#     return blob.exists()

# Check if a directory exists in the bucket
# def dir_exists(path):
#     if not path.endswith("/"):
#         path += "/"
    
#     blobs = bucket.list_blobs(prefix=trim_path(path), delimiter="/", max_results=1)
#     return any(True for _ in blobs)

'''
# The internal DNS within a terra instance redirects googleapis.com -> restricted.google.com and bans the request
# Don't mess with /etc/resolv.conf, it seems to have side effects
# Modify the python to use dns_server=8.8.8.8 without affecting the rest (hopefully)
import dns.resolver
import socket
def custom_getaddrinfo(host, port, family=0, type=0, proto=0, flags=0, dns_server='8.8.8.8'):
    resolver = dns.resolver.Resolver()
    resolver.nameservers = [dns_server]
    try:
        answer_ipv4 = resolver.resolve(host, 'A')
    except (dns.resolver.NoAnswer, dns.resolver.NXDOMAIN, dns.resolver.Timeout):
        answer_ipv4 = []
    try:
        answer_ipv6 = resolver.resolve(host, 'AAAA')
    except (dns.resolver.NoAnswer, dns.resolver.NXDOMAIN, dns.resolver.Timeout):
        answer_ipv6 = []
    if not answer_ipv4 and not answer_ipv6:
        raise socket.gaierror(socket.EAI_NONAME, 'Name or service not known')
    addresses = [(socket.AF_INET, a.to_text()) for a in answer_ipv4] + [(socket.AF_INET6, a.to_text()) for a in answer_ipv6]
    results = []
    for af, addr in addresses:
        if family != 0 and family != af:
            continue
        if type == socket.SOCK_STREAM:
            proto = socket.IPPROTO_TCP
        elif type == socket.SOCK_DGRAM:
            proto = socket.IPPROTO_UDP
        else:
            raise ValueError("Unknown socket type")
        results.append((af, type, proto, '', (addr, port)))
    if not results:
        raise socket.gaierror(socket.EAI_FAMILY, 'Address family not supported')
    return results
# Replace the default getaddrinfo with custom_getaddrinfo
socket.getaddrinfo = custom_getaddrinfo
'''

'''
workflow="cellranger-count"
bcl="250527_SL-EXB_0558_B22Y57YLT4"
index="all"

workflow="recon"
bcl="250521_SL-EXA_0515_A22YGT5LT4"
index="all"

workflow="slide-tags"
bcl="250527_SL-EXB_0558_B22Y57YLT4"
index="all"
'''

# Load the worksheet
if workflow in ["cellranger-count", "slide-tags"]:
    df = pd.read_csv('https://docs.google.com/spreadsheets/d/1NOaWXARQiSA6fquOtcouQPREPN4buYIf13tq_F6D9As/export?format=csv&gid=0')
    cols = ["BCL", "Reference", "RNAIndex", "SBIndex", "Puck"]
elif workflow in ["recon", "reconstruction"]:
    df = pd.read_csv('https://docs.google.com/spreadsheets/d/1NOaWXARQiSA6fquOtcouQPREPN4buYIf13tq_F6D9As/export?format=csv&gid=229934426')
    cols = ["BCL", "Index", "bc1", "bc2"]
else:
    sys.exit(f"Unknown workflow: {workflow}")


# Clean the worksheet
assert all(col in df.columns for col in cols)
df = df[cols].apply(lambda col: col.map(str.strip) if pd.api.types.is_string_dtype(col) else col)
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
else:
    sys.exit(f"Unknown workflow: {workflow}")
    
def assert_unique_column(series):
    assert series.notna().all(), f"Column has NA values:\n{series}"
    assert (~series.str.strip().eq("")).all(), f"Column has empty values:\n{series}"
    assert series.is_unique, f"Column has repeated values:\n{series}"

assert_unique_column(df[idx_col])
df = df if index == "all" else df[df[idx_col] == index]
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
    assert_unique_column(df["RNAIndex"])
    count_blobs = bucket.list_blobs(prefix=f"gene-expression/{bcl}")
    counts = {blob.name.split("/")[2] for blob in count_blobs}
    assert not df["RNAIndex"].isin(counts).any(), f"Output {set(df['RNAIndex'])&counts} already exists in the bucket"

elif workflow == "slide-tags":
    # Assert RNA input exists
    assert_unique_column(df["RNAIndex"])
    count_blobs = bucket.list_blobs(prefix=f"gene-expression/{bcl}")
    counts = {blob.name.split("/")[2] for blob in count_blobs}
    assert df["RNAIndex"].isin(counts).all(), f"GEX for {set(df['RNAIndex'])-counts} does not exist in the bucket"

    # Assert SB input exists
    assert_full_column(df["SBIndex"])
    bucket.list_blobs(prefix=f"fastqs/{bcl}")

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

    pucks = apply.pucks_to_URIs(df["Puck"])
    assert False

elif workflow in ["recon", "reconstruction"]:
    assert_unique_column(df["Index"])
else:
    sys.exit(f"Unknown workflow: {workflow}")


# Compute memory requirements
def compute_fastq_sizes(series):
    assert_full_column(series)
    fastq_blobs = bucket.list_blobs(prefix=f"fastqs/{bcl}")
    fastqs = [(blob.name, blob.size) for blob in fastq_blobs if blob.name.endswith(".fastq.gz")]
    mem_GBs = [int(sum(s for n,s in fastqs if "/"+i+"_S" in n) / 1e9) for i in series]
    return mem_GBs
    
if workflow == "cellranger-count":
    # Compute the FASTQ sizes
    mem_GBs = compute_fastq_sizes(df["RNAIndex"])
    mem_GBs = [1.5*mem for mem in mem_GBs]

elif workflow == "slide-tags":
    assert False
    
elif workflow in ["recon", "reconstruction"]:
    # Compte the FASTQ sizes
    assert_full_column(df["Index"])
    mem_GBs_fastq = compute_fastq_sizes(df["Index"])
    mem_GBs_fastq = [int(2*mem) for mem in mem_GBs_fastq]
    
    # Compute the intermediate matrix sizes
    recon_blobs = bucket.list_blobs(prefix=f"recon/{bcl}")
    recons = [(blob.name, blob.size) for blob in recon_blobs if blob.name.endswith("knn2.npz")]
    mem_GBs_mat = [max((s for n,s in recons if "/"+i+"/" in n), default=0) / 1e9 for i in df["Index"]]
    mem_GBs_mat = [int(25*mem) for mem in mem_GBs_mat]

    # Take the minimum available
    mem_GBs = [min(x,y) if min(x,y) > 0 else max(x,y) for x,y in zip(mem_GBs_fastq, mem_GBs_mat)]
    
else:
    sys.exit(f"Unknown workflow: {workflow}")

assert all(m > 0 for m in mem_GBs), f"Incomplete memory estimation: {mem_GBs} (missing input files)"
print(f"Memory (GB): {mem_GBs}")


# Compute the names, assert the jobs are not already running
job_names = ["_".join([workflow, idx, bcl]) for idx in df[idx_col]]
assert len(job_names) == len(set(job_names))

subs = fapi.list_submissions("testmybroad", "Macosko-Pipelines").json()
subs = [sub for sub in subs if sub["status"] not in ["Done","Aborted"]]
running_job_names = [sub["userComment"] for sub in subs if "userComment" in sub]
assert not set(job_names) & set(running_job_names), "Jobs already running!"


# Print + exit if testing
if dryrun:
    print(df)
    print("Dry run complete, no errors found - exiting...")
    sys.exit(0)

### Terra Submission ###########################################################

wnamespace = "testmybroad"
workspace = "Macosko-Pipelines"
cnamespace = "macosko-pipelines"

def run_cellranger_count(bcl, index, reference, mem_GB, disk_GB, params="", user_comment=""):
    body = fapi.get_workspace_config(wnamespace, workspace, cnamespace, "cellranger-count").json()
    body["inputs"]["cellranger_count.bcl"] = f'"{bcl}"'
    body["inputs"]["cellranger_count.index"] = f'"{index}"'
    body["inputs"]["cellranger_count.reference"] = f'"{reference}"'
    body["inputs"]["cellranger_count.mem_GB"] = f'"{mem_GB}"'
    body["inputs"]["cellranger_count.disk_GB"] = f'"{disk_GB}"'
    body["inputs"]["cellranger_count.params"] = f'"{params}"'
    body["inputs"]["cellranger_count.docker"] = f''
    res = fapi.update_workspace_config(wnamespace, workspace, cnamespace, "cellranger-count", body)
    assert res.status_code == 200, res.json()['message']
    
    submit("cellranger-count", user_comment)
    return True

def run_reconstruction(bcl, index, mem_GB, disk_GB, params="", user_comment=""):
    body = fapi.get_workspace_config(wnamespace, workspace, cnamespace, "reconstruction").json()  
    body["inputs"]["reconstruction.bcl"] = f'"{bcl}"'
    body["inputs"]["reconstruction.index"] = f'"{index}"'
    body["inputs"]["reconstruction.mem_GB"] = f'"{mem_GB}"'
    body["inputs"]["reconstruction.disk_GB"] = f'"{disk_GB}"'
    #body["inputs"]["reconstruction.params"] = f'"{params}"'
    #body["inputs"]["reconstruction.docker"] = f''
    res = fapi.update_workspace_config(wnamespace, workspace, cnamespace, "reconstruction", body)
    assert res.status_code == 200, res.json()['message']
    
    submit("reconstruction", user_comment)
    return True
    
def submit(config, user_comment=""):
    # Validate the configuration
    res = fapi.validate_config(wnamespace, workspace, cnamespace, config).json()
    assert res["extraInputs"] == [], f"ERROR: extra input: \n{res['extraInputs']}"
    assert res["invalidInputs"] == {}, f"ERROR: invalid input: \n{res['invalidInputs']}"
    assert res["invalidOutputs"] == {}, f"ERROR: invalid output: \n{res['invalidOutputs']}"
    assert res["missingInputs"] == [], f"ERROR: missing input: \n{res['missingInputs']}"
    
    # Submit the job
    res = fapi.create_submission(wnamespace, workspace, cnamespace, config, user_comment=user_comment)
    assert res.status_code == 201, res.status_code + ": " + res.json()['message']
    print(f"Submitted {config} {user_comment}")


if workflow == "cellranger-count":
    pass
elif workflow == "slide-tags":
    pass
elif workflow in ["recon", "reconstruction"]:
    for row, mem, job in zip(df.itertuples(index=False), mem_GBs, job_names):
        run_reconstruction(row.BCL, row.Index, mem, mem, params="", user_comment=job)
else:
    sys.exit(f"Unknown workflow: {workflow}")

### Terra Submission ###########################################################

# List all submissions
subs = fapi.list_submissions("testmybroad", "Macosko-Pipelines").json()
subs = [sub for sub in subs if sub["status"] not in ["Done","Aborted"]]
print(f"Currently running submissions: {len(subs)}")

# Abort all submissions
# ids = [sub["submissionId"] for sub in subs]
# [fapi.abort_submission(wnamespace, workspace, submission_id) for submission_id in ids]
