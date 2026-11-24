version 1.0

# Scatter/gather version of reconstruction.wdl. recon-count.jl is split into the
# recon-count_*.jl steps (see reconstruction/README.md) so no single VM has to hold every
# read at once; knn.py and recon.py then run as their own tasks with their own memory.

task setup {
    input {
        String bcl
        String index
        Int bc1
        Int bc2
        Int lanes
        String docker
        String bucket
        String branch
        Int? pr
        String? subfolder
        String? selection
    }
    command <<<
    set -euo pipefail

    # Resolve the ref once and hand the scripts to every task, so all shards run the same commit
    BRANCH="~{branch}"
    PR="~{default="" pr}"
    if [ -n "$PR" ]; then REF="refs/pull/$PR/head"; else REF="$BRANCH"; fi
    wget -q -O repo.tar.gz "https://codeload.github.com/MacoskoLab/Macosko-Pipelines/tar.gz/$REF"
    mkdir repo && tar xzf repo.tar.gz -C repo --strip-components=1

    # tools/ is only needed to cut a selection, so a plain run never fails over its absence
    SELECTION="~{default="" selection}"
    if [ -n "$SELECTION" ]; then
        cp repo/tools/puck-select.py repo/tools/subset-matrix.py repo/reconstruction/ 2>/dev/null || {
            echo "ERROR: tools/puck-select.py + tools/subset-matrix.py are not present at ref '$REF'"
            exit 1; }
    fi
    tar czf scripts.tar.gz -C repo/reconstruction .

    BUCKET="~{bucket}"

    SUBFOLDER="~{default="" subfolder}"
    if [ -z "$SUBFOLDER" ]; then
        if   [ -n "$PR" ];            then SUBFOLDER="pr-$PR"
        elif [ "$BRANCH" != "main" ]; then SUBFOLDER="$BRANCH"
        fi
    fi

    work_dir="gs://$BUCKET/recon/~{bcl}/~{index}-~{lanes}" ; work_dir=${work_dir%-12345678}
    if [ -n "$SUBFOLDER" ]; then work_dir="$work_dir/$SUBFOLDER"; fi
    echo "$work_dir" > work_dir.txt
    echo "Output directory: $work_dir"

    # A selection is a curated subset of the beads, cut by tools/puck-select.py and kept in its
    # own subfolder beside the full run it was cut from. Checked here so a missing sidecar fails
    # in minutes rather than after the counting stages
    sel_dir="$work_dir"
    if [ -n "$SELECTION" ]; then
        sel_dir="$work_dir/$SELECTION"
        gsutil -q stat "$sel_dir/selection.json" || {
            echo "ERROR: no $sel_dir/selection.json - run tools/puck-select.py first"; exit 1; }
        echo "Selection directory: $sel_dir"
    fi
    echo "$sel_dir" > sel_dir.txt

    # Only this run's own folder counts as a cache: unlike reconstruction.wdl there is no
    # fall-through to the root run, since the point of a scatter run is to count afresh
    cached=true
    for f in matrix.csv.gz sb1.txt.gz sb2.txt.gz metadata.csv; do
        gsutil -q stat "$work_dir/$f" || cached=false
    done
    if [ "$cached" = true ]; then
        meta=$(gcloud storage cat "$work_dir/metadata.csv")
        m1=$(grep -m1 '^R1_barcodes_manual,' <<< "$meta" | cut -d',' -f2 || true)
        m2=$(grep -m1 '^R2_barcodes_manual,' <<< "$meta" | cut -d',' -f2 || true)
        if ! { [ -n "$m1" ] && [ -n "$m2" ] && [ "$m1" -eq ~{bc1} ] && [ "$m2" -eq ~{bc2} ]; }; then
            cached=false
        fi
    fi
    echo "$cached" > cached.txt
    echo "Cached matrix: $cached"
    >>>
    output {
        File scripts = "scripts.tar.gz"
        String work_dir = read_string("work_dir.txt")
        String sel_dir = read_string("sel_dir.txt")
        Boolean cached = read_boolean("cached.txt")
    }
    runtime {
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 10 HDD"
        docker: docker
        preemptible: 2
    }
}

task find_fastqs {
    input {
        File scripts
        String bcl
        String index
        Int lanes
        String docker
        String bucket
    }
    command <<<
    set -euo pipefail
    mkdir src && tar xzf ~{scripts} -C src

    # learn-bead-types.jl only looks at the first 100000 reads of each FASTQ, so fetch just
    # the head of each file instead of downloading them whole. A byte range rather than
    # cat | head, since gcloud does not exit when head closes the pipe and the VM hangs.
    # 50 MB compressed holds ~7x the reads needed; the truncated gzip tail is discarded.
    BUCKET="~{bucket}"
    mkdir fastqs
    set +o pipefail
    for f in $(gcloud storage ls "gs://$BUCKET/fastqs/~{bcl}/~{index}_*_L00[~{lanes}]_*"); do
        out="fastqs/$(basename "$f")"
        gcloud storage cat -r 0-49999999 "$f" | (zcat 2>/dev/null || true) | head -n 400004 | gzip > "$out"
        n=$(zcat "$out" | wc -l)
        echo "$(basename "$f"): $n lines"
        # Fewer lines is fine for a small FASTQ, but never zero or a partial record
        [ "$n" -gt 0 ] && [ $((n % 4)) -eq 0 ] || { echo "ERROR: bad read of $f"; exit 1; }
    done
    set -o pipefail

    mkdir out
    julia src/recon-count_find-fastqs.jl fastqs out -r '_L00[~{lanes}]_'
    >>>
    output {
        File metadata = "out/metadata.csv"
        File fastqs1 = "out/fastqs1.txt"
        File fastqs2 = "out/fastqs2.txt"
        File filepaths_pdf = "out/filepaths.pdf"
    }
    runtime {
        cpu: 2
        memory: "8 GB"
        disks: "local-disk 20 HDD"
        docker: docker
        preemptible: 2
    }
}

task read_fastqs {
    input {
        File scripts
        File metadata
        File fastqs1
        File fastqs2
        String bcl
        String index
        Int lanes
        Float downsample_prop
        String R1_filter
        String R2_filter
        String docker
        String bucket
        Int mem_GB
        Int disk_GB
    }
    command <<<
    set -euo pipefail
    mkdir src && tar xzf ~{scripts} -C src

    BUCKET="~{bucket}"
    mkdir fastqs
    gcloud storage cp "gs://$BUCKET/fastqs/~{bcl}/~{index}_*_L00[~{lanes}]_*" fastqs

    # fastqs1/2.txt hold the absolute paths find_fastqs saw, so repoint them at this VM's copies
    mkdir in out
    cp ~{metadata} in/metadata.csv
    sed "s|.*/|$PWD/fastqs/|" ~{fastqs1} > in/fastqs1.txt
    sed "s|.*/|$PWD/fastqs/|" ~{fastqs2} > in/fastqs2.txt

    julia src/recon-count_read-fastqs.jl in out -p ~{downsample_prop} -f ~{R1_filter} -g ~{R2_filter}
    rm -rf fastqs
    >>>
    output {
        Array[File] reads_dfs = glob("out/reads_df.*.csv.gz")
        Array[File] rpsbs = glob("out/rpsb*.csv.gz")
        Array[File] metadatas = glob("out/metadata.*.csv")
    }
    runtime {
        cpu: 8
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 2
    }
}

task compute_whitelists {
    input {
        File scripts
        File metadata
        Array[File] rpsbs
        Int bc1
        Int bc2
        String docker
        Int mem_GB
        Int disk_GB
    }
    command <<<
    set -euo pipefail
    # Julia's logged "GB" is live memory after a GC, not the peak, so sample the VM's used
    # memory and report the real peak, even when the task is OOM-killed
    ( while sleep 10; do awk '/^MemTotal/{t=$2} /^MemAvailable/{a=$2} END{print int((t-a)/1024)}' /proc/meminfo >> mem.log; done ) & MON=$!
    trap 'kill $MON 2>/dev/null; echo "Peak memory used: $(sort -n mem.log 2>/dev/null | tail -1) MB"' EXIT
    mkdir src && tar xzf ~{scripts} -C src

    mkdir in out
    cp ~{metadata} in/metadata.csv
    while read -r f; do ln -s "$f" in/; done < ~{write_lines(rpsbs)}

    julia src/recon-count_compute-whitelists.jl in out -x ~{bc1} -y ~{bc2}
    >>>
    output {
        File metadata_out = "out/metadata.csv"
        File wl1 = "out/wl1.txt.gz"
        File wl2 = "out/wl2.txt.gz"
        File elbows_pdf = "out/elbows.pdf"
    }
    runtime {
        cpu: 2
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 2
    }
}

# Shared by match-whitelists, compute-chimeras and remove-chimeras: each links its inputs into
# one directory, runs the step for one prefix filter, and globs back whatever it wrote
task prefix_step {
    input {
        File scripts
        String step
        Array[File] inputs
        String R1_filter = ""
        String R2_filter = ""
        String docker
        Int mem_GB
        Int disk_GB
    }
    command <<<
    set -euo pipefail
    # Julia's logged "GB" is live memory after a GC, not the peak, so sample the VM's used
    # memory and report the real peak, even when the task is OOM-killed
    ( while sleep 10; do awk '/^MemTotal/{t=$2} /^MemAvailable/{a=$2} END{print int((t-a)/1024)}' /proc/meminfo >> mem.log; done ) & MON=$!
    trap 'kill $MON 2>/dev/null; echo "Peak memory used: $(sort -n mem.log 2>/dev/null | tail -1) MB"' EXIT
    mkdir src && tar xzf ~{scripts} -C src

    mkdir in out
    while read -r f; do ln -s "$f" in/; done < ~{write_lines(inputs)}

    FILTERS=""
    if [ -n "~{R1_filter}" ]; then FILTERS="$FILTERS -f ~{R1_filter}"; fi
    if [ -n "~{R2_filter}" ]; then FILTERS="$FILTERS -g ~{R2_filter}"; fi
    julia src/recon-count_~{step}.jl in out $FILTERS
    >>>
    output {
        Array[File] outputs = glob("out/*")
    }
    runtime {
        cpu: 2
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 2
    }
}

task count_umis {
    input {
        File scripts
        Array[File] inputs
        Int bc1
        Int bc2
        Float downsample_prop
        String work_dir
        String docker
        Int mem_GB
        Int disk_GB
    }
    command <<<
    set -euo pipefail
    # Julia's logged "GB" is live memory after a GC, not the peak, so sample the VM's used
    # memory and report the real peak, even when the task is OOM-killed
    ( while sleep 10; do awk '/^MemTotal/{t=$2} /^MemAvailable/{a=$2} END{print int((t-a)/1024)}' /proc/meminfo >> mem.log; done ) & MON=$!
    trap 'kill $MON 2>/dev/null; echo "Peak memory used: $(sort -n mem.log 2>/dev/null | tail -1) MB"' EXIT
    mkdir src && tar xzf ~{scripts} -C src

    mkdir in out
    while read -r f; do ln -s "$f" in/; done < ~{write_lines(inputs)}

    julia src/recon-count_count-umis.jl in out -x ~{bc1} -y ~{bc2} -p ~{downsample_prop}

    # Uploaded here rather than after recon.py, so a later OOM does not discard the counting
    echo "----- Uploading recon-count output -----"
    gcloud storage cp out/* "~{work_dir}/"
    >>>
    output {
        String done = "~{work_dir}"
    }
    runtime {
        cpu: 2
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 2
    }
}

task select_beads {
    input {
        File scripts
        String base_dir
        String sel_dir
        String docker
        Int mem_GB
        Int disk_GB
    }
    command <<<
    set -euo pipefail
    ( while sleep 10; do awk '/^MemTotal/{t=$2} /^MemAvailable/{a=$2} END{print int((t-a)/1024)}' /proc/meminfo >> mem.log; done ) & MON=$!
    trap 'kill $MON 2>/dev/null; echo "Peak memory used: $(sort -n mem.log 2>/dev/null | tail -1) MB"' EXIT
    mkdir src && tar xzf ~{scripts} -C src
    PY="/root/.local/bin/micromamba run python"

    mkdir base cache
    gcloud storage cp "~{base_dir}/matrix.csv.gz" "~{base_dir}/sb1.txt.gz" \
                      "~{base_dir}/sb2.txt.gz"    "~{base_dir}/metadata.csv" \
                      "~{sel_dir}/selection.json" base

    # Replay the cut recorded by tools/puck-select.py against the puck it was drawn on, which
    # recon.py wrote inside its own UMAP..._ne2000 folder of the full run
    SOURCE_PUCK=$($PY -c "import json; print(json.load(open('base/selection.json'))['source_puck'])")
    gcloud storage cp "~{base_dir}/$SOURCE_PUCK/Puck.csv" base/Puck.csv
    $PY src/puck-select.py base/Puck.csv -o base/Puck-selected.csv --non_interactive \
        --from_json base/selection.json
    $PY src/subset-matrix.py -i base -o cache -p base/Puck-selected.csv -b 2 -l 1
    cp base/Puck-selected.csv base/selection.json cache/
    gcloud storage cp cache/* "~{sel_dir}/"
    >>>
    output {
        String done = "~{sel_dir}"
    }
    runtime {
        cpu: 2
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 2
    }
}

task knn {
    input {
        File scripts
        String work_dir
        String docker
        Int mem_GB
        Int cpu
        Int disk_GB
    }
    command <<<
    set -euo pipefail
    # Julia's logged "GB" is live memory after a GC, not the peak, so sample the VM's used
    # memory and report the real peak, even when the task is OOM-killed
    ( while sleep 10; do awk '/^MemTotal/{t=$2} /^MemAvailable/{a=$2} END{print int((t-a)/1024)}' /proc/meminfo >> mem.log; done ) & MON=$!
    trap 'kill $MON 2>/dev/null; echo "Peak memory used: $(sort -n mem.log 2>/dev/null | tail -1) MB"' EXIT
    mkdir src && tar xzf ~{scripts} -C src
    PY="/root/.local/bin/micromamba run python"

    mkdir cache
    gcloud storage cp "~{work_dir}/matrix.csv.gz" "~{work_dir}/sb1.txt.gz" \
                      "~{work_dir}/sb2.txt.gz"    "~{work_dir}/metadata.csv" cache

    $PY src/knn.py -i cache -o cache -b 2 -k 2
    gcloud storage cp cache/knn2.npz "~{work_dir}/"
    >>>
    output {
        String done = "~{work_dir}"
    }
    runtime {
        cpu: cpu
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 0
    }
}

task recon {
    input {
        File scripts
        String work_dir
        String params
        Boolean selected
        String docker
        Int mem_GB
        Int cpu
        Int disk_GB
    }
    command <<<
    set -euo pipefail
    # Julia's logged "GB" is live memory after a GC, not the peak, so sample the VM's used
    # memory and report the real peak, even when the task is OOM-killed
    ( while sleep 10; do awk '/^MemTotal/{t=$2} /^MemAvailable/{a=$2} END{print int((t-a)/1024)}' /proc/meminfo >> mem.log; done ) & MON=$!
    trap 'kill $MON 2>/dev/null; echo "Peak memory used: $(sort -n mem.log 2>/dev/null | tail -1) MB"' EXIT
    mkdir src && tar xzf ~{scripts} -C src
    PY="/root/.local/bin/micromamba run python"

    mkdir cache
    gcloud storage cp "~{work_dir}/matrix.csv.gz" "~{work_dir}/sb1.txt.gz" \
                      "~{work_dir}/sb2.txt.gz"    "~{work_dir}/metadata.csv" \
                      "~{work_dir}/knn2.npz" cache

    # estimate_diameter() bins on total bead count, so a large cut silently rescales the puck.
    # subset-matrix.py records the pre-cut bin; carry it over unless params already sets one
    DIAM=""
    if [ "~{selected}" = true ] && ! grep -q -- '-D' <<< "~{params}"; then
        prev=$(grep -m1 '^subset_diameter_prev,' cache/metadata.csv | cut -d',' -f2 || true)
        new=$( grep -m1 '^subset_diameter_new,'  cache/metadata.csv | cut -d',' -f2 || true)
        if [ -n "$prev" ] && [ "$prev" != "None" ] && [ "$prev" != "$new" ]; then
            DIAM="-D $prev"
            echo "NOTE: subsetting changed the estimated diameter ${prev} -> ${new}um, preserving ${prev}um"
        fi
    fi

    $PY src/recon.py -i cache -o output -b 2 $DIAM ~{params}

    echo "----- Uploading results -----"
    gcloud storage cp -r output/* "~{work_dir}/"
    >>>
    runtime {
        cpu: cpu
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 0
    }
}

workflow reconstruction {
    input {
        String bcl
        String index
        # Not used: submit.py always writes mem_GB (sized for the serial run, so FASTQ-dominated
        # and far too big here), and the per-stage inputs below set every task's memory instead
        Int? mem_GB
        Int disk_GB
        Int bc1 = 0
        Int bc2 = 0
        Int lanes = 12345678
        String params = ""
        String docker = "us-central1-docker.pkg.dev/velina-208320/terra/pipeline-image:latest"
        # Workspace bucket holding fastqs/ and recon/; the BICAN config sets its own
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        Float downsample_prop = 1.0
        String branch = "main"
        Int? pr
        String? subfolder
        # Name of a tools/puck-select.py cut under <run dir>/<selection>/selection.json. The full
        # run's matrix is reused if present (counted first if not), then subset, then knn + recon
        String? selection

        # Per-stage memory: about 1.3-1.6x the "Peak memory used" each task printed on D707D504
        # (514 GB of FASTQs, 100% of reads; peaks 27/14/53/29/60/63/40 GB). Julia's own logged
        # "GB" is live memory after a GC and far below the peak, so don't size from it. Larger
        # pucks may need these raised in the config. read_fastqs has no peak measurement yet and
        # keeps the 64 GB it has always run with
        Int mem_read_GB = 64
        Int mem_whitelist_GB = 40
        Int mem_match_GB = 24
        Int mem_chimera_GB = 72
        Int mem_remove_GB = 40
        Int mem_count_GB = 80
        Int mem_knn_GB = 96
        Int mem_recon_GB = 64
        Int mem_select_GB = 32

        # knn/recon are multi-threaded, so their cores are set explicitly rather than left to
        # follow memory. recon's UMAP embedding scaled well (63 -> 29 min from ~20 to 48 vCPUs);
        # knn barely did (2h01 -> 1h35 with all 48 busy), so it gets fewer
        Int cpu_knn = 32
        Int cpu_recon = 48

        # disk_GB (sized off the FASTQs) only goes to read_fastqs, which downloads them all;
        # every later task works on a slice of the reads and gets disk_work_GB
        Int disk_work_GB = 200
    }

    Array[String] bases = ["A", "C", "G", "T"]
    Array[String] bases_n = ["A", "C", "G", "T", "N"]

    call setup {
        input:
            bcl = bcl, index = index, bc1 = bc1, bc2 = bc2, lanes = lanes,
            docker = docker, bucket = bucket, branch = branch, pr = pr, subfolder = subfolder, selection = selection
    }

    if (!setup.cached) {
        call find_fastqs {
            input:
                scripts = setup.scripts, bcl = bcl, index = index, lanes = lanes, docker = docker, bucket = bucket
        }

        # Every prefix pair, N included, so N reads are still tallied in the metadata
        scatter (p in cross(bases_n, bases_n)) {
            call read_fastqs {
                input:
                    scripts = setup.scripts, metadata = find_fastqs.metadata,
                    fastqs1 = find_fastqs.fastqs1, fastqs2 = find_fastqs.fastqs2,
                    bcl = bcl, index = index, lanes = lanes, downsample_prop = downsample_prop,
                    R1_filter = p.left, R2_filter = p.right,
                    docker = docker, bucket = bucket, mem_GB = mem_read_GB, disk_GB = disk_GB
            }
        }

        call compute_whitelists {
            input:
                scripts = setup.scripts, metadata = find_fastqs.metadata,
                rpsbs = flatten(read_fastqs.rpsbs), bc1 = bc1, bc2 = bc2,
                docker = docker, mem_GB = mem_whitelist_GB, disk_GB = 50
        }

        # Shards are addressed by index rather than handed every file, so each VM only
        # localizes its own slice: read_fastqs shard 5*a+b holds prefix pair (a, b) over
        # ACGTN, and match/remove shard 4*a+b holds the same pair over ACGT.

        # N reads can never match a whitelist, so only the ACGT pairs go further
        # Terra's WDL parser only takes a plain identifier inside [...], hence the named indices
        scatter (i in range(16)) {
            Int a = i / 4
            Int b = i % 4
            Int r = 5 * a + b
            call prefix_step as match_whitelists {
                input:
                    scripts = setup.scripts, step = "match-whitelists",
                    inputs = flatten([read_fastqs.reads_dfs[r], read_fastqs.metadatas[r],
                                      [compute_whitelists.wl1, compute_whitelists.wl2]]),
                    R1_filter = bases[a], R2_filter = bases[b],
                    docker = docker, mem_GB = mem_match_GB, disk_GB = disk_work_GB
            }
        }

        # Chimeras are called per sb1 prefix (across every sb2) and per sb2 prefix
        scatter (c in range(4)) {
            Int r0 = 4 * c
            Int r1 = 4 * c + 1
            Int r2 = 4 * c + 2
            Int r3 = 4 * c + 3
            Int c1 = 4 + c
            Int c2 = 8 + c
            Int c3 = 12 + c
            call prefix_step as compute_chimeras1 {
                input:
                    scripts = setup.scripts, step = "compute-chimeras",
                    inputs = flatten([match_whitelists.outputs[r0], match_whitelists.outputs[r1],
                                      match_whitelists.outputs[r2], match_whitelists.outputs[r3]]),
                    R1_filter = bases[c], docker = docker, mem_GB = mem_chimera_GB, disk_GB = disk_work_GB
            }
            call prefix_step as compute_chimeras2 {
                input:
                    scripts = setup.scripts, step = "compute-chimeras",
                    inputs = flatten([match_whitelists.outputs[c],  match_whitelists.outputs[c1],
                                      match_whitelists.outputs[c2], match_whitelists.outputs[c3]]),
                    R2_filter = bases[c], docker = docker, mem_GB = mem_chimera_GB, disk_GB = disk_work_GB
            }
        }

        # The chimeric files only list chimeric reads, so every shard gets all of them, as in
        # test_recon_count_scatter.sh
        scatter (i in range(16)) {
            Int ra = i / 4
            Int rb = i % 4
            call prefix_step as remove_chimeras {
                input:
                    scripts = setup.scripts, step = "remove-chimeras",
                    inputs = flatten([match_whitelists.outputs[i], flatten(compute_chimeras1.outputs),
                                      flatten(compute_chimeras2.outputs)]),
                    R1_filter = bases[ra], R2_filter = bases[rb],
                    docker = docker, mem_GB = mem_remove_GB, disk_GB = disk_work_GB
            }
        }

        # Mirrors the inputs test_recon_count_scatter.sh hands count-umis: the N-shard
        # metadata from read-fastqs (the ACGT ones are superseded by remove-chimeras')
        scatter (i in range(25)) {
            if (i / 5 == 4 || i % 5 == 4) {
                Array[File] n_metadata = read_fastqs.metadatas[i]
            }
        }

        call count_umis {
            input:
                scripts = setup.scripts,
                inputs = flatten([[find_fastqs.filepaths_pdf, compute_whitelists.elbows_pdf,
                                   compute_whitelists.metadata_out],
                                  flatten(select_all(n_metadata)), flatten(remove_chimeras.outputs)]),
                bc1 = bc1, bc2 = bc2, downsample_prop = downsample_prop, work_dir = setup.work_dir,
                docker = docker, mem_GB = mem_count_GB, disk_GB = disk_work_GB
        }
    }

    # The base matrix comes from count_umis on a fresh run, or is already in work_dir
    String base_dir = select_first([count_umis.done, setup.work_dir])

    if (defined(selection)) {
        call select_beads {
            input:
                scripts = setup.scripts, base_dir = base_dir, sel_dir = setup.sel_dir,
                docker = docker, mem_GB = mem_select_GB, disk_GB = disk_work_GB
        }
    }

    call knn {
        input:
            scripts = setup.scripts, work_dir = select_first([select_beads.done, base_dir]),
            docker = docker, mem_GB = mem_knn_GB, cpu = cpu_knn, disk_GB = disk_work_GB
    }

    call recon {
        input:
            scripts = setup.scripts, work_dir = knn.done, params = params, selected = defined(selection),
            docker = docker, mem_GB = mem_recon_GB, cpu = cpu_recon, disk_GB = disk_work_GB
    }
}
