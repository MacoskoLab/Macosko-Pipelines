version 1.0

task recon {
    input {
        String bcl
        String index
        Int mem_GB
        Int disk_GB
        Int bc1
        Int bc2
        Int lanes
        String params
        String docker
        Float downsample_prop
        String branch
        Int? pr
        String? subfolder
        String? selection
    }
    command <<<
    set -euo pipefail

    # Resolve the git ref to pull scripts from: a PR number (refs/pull/<n>/head) wins,
    # otherwise ~{branch} is used as a raw ref (branch/tag name or commit SHA).
    BRANCH="~{branch}"
    PR="~{default="" pr}"
    if [ -n "$PR" ]; then REF="refs/pull/$PR/head"; else REF="$BRANCH"; fi

    # A selection is a curated subset of the beads, living in its own subfolder alongside
    # the full reconstruction it was cut from. Resolved here because it decides whether the
    # tools/ scripts below are needed at all.
    SELECTION="~{default="" selection}"

    # Download the entire reconstruction/ folder for this ref into the working dir
    # (recon-count.jl includes other .jl files, so fetching individual files is not enough)
    wget -q -O repo.tar.gz "https://codeload.github.com/MacoskoLab/Macosko-Pipelines/tar.gz/$REF"
    mkdir repo && tar xzf repo.tar.gz -C repo --strip-components=1
    cp repo/reconstruction/* .

    # tools/ is only needed to build a selection. Guarded so that a plain run never fails
    # over scripts it will not execute, and so a selection run says why if they are absent.
    if [ -n "$SELECTION" ]; then
        cp repo/tools/puck-select.py repo/tools/subset-matrix.py . 2>/dev/null || {
            echo "ERROR: tools/puck-select.py + tools/subset-matrix.py are not present at ref '$REF'"
            echo "       commit and push them, or point branch/pr at a ref that has them"
            exit 1; }
    fi

    rm -rf repo repo.tar.gz

    PY="/root/.local/bin/micromamba run python"

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/~{bcl}"

    # Determine the output subfolder: explicit override wins, else a PR names it pr-<n>,
    # else a non-main branch names it, else none.
    SUBFOLDER="~{default="" subfolder}"
    if [ -z "$SUBFOLDER" ]; then
        if   [ -n "$PR" ];            then SUBFOLDER="pr-$PR"
        elif [ "$BRANCH" != "main" ]; then SUBFOLDER="$BRANCH"
        fi
    fi

    base_dir="gs://$BUCKET/recon/~{bcl}/~{index}-~{lanes}" ; base_dir=${base_dir%-12345678}
    if [ -n "$SUBFOLDER" ]; then
        base_dir="$base_dir/$SUBFOLDER"
    fi

    # The selection subfolder is appended after the PR/branch one, so a PR test and a
    # selection compose rather than clobbering each other.
    if [ -n "$SELECTION" ]; then work_dir="$base_dir/$SELECTION"; else work_dir="$base_dir"; fi

    echo "==================== START RECONSTRUCTION ===================="
    echo "Output directory: $work_dir"

    # True when $1 holds a complete recon-count.jl output set matching the requested barcode
    # counts. The greps are guarded because a missing key would otherwise trip set -e.
    count_cached() {
        local dir="$1" meta m1 m2
        for f in matrix.csv.gz sb1.txt.gz sb2.txt.gz metadata.csv; do
            gsutil -q stat "$dir/$f" || return 1
        done
        meta=$(gcloud storage cat "$dir/metadata.csv")
        m1=$(grep -m1 '^R1_barcodes_manual,' <<< "$meta" | cut -d',' -f2 || true)
        m2=$(grep -m1 '^R2_barcodes_manual,' <<< "$meta" | cut -d',' -f2 || true)
        [ -n "$m1" ] && [ -n "$m2" ] && [ "$m1" -eq ~{bc1} ] && [ "$m2" -eq ~{bc2} ]
    }

    ##### Stage 1: the diffusion matrix #####

    mkdir cache
    if count_cached "$work_dir"; then
        echo "----- Using cached matrix from $work_dir -----"
        gcloud storage cp "$work_dir/matrix.csv.gz" "$work_dir/sb1.txt.gz" \
                          "$work_dir/sb2.txt.gz"    "$work_dir/metadata.csv" cache
        FRESH_MATRIX=false
    else
        FRESH_MATRIX=true
        mkdir base

        if [ "$work_dir" != "$base_dir" ] && count_cached "$base_dir"; then
            echo "----- Using cached matrix from $base_dir -----"
            gcloud storage cp "$base_dir/matrix.csv.gz" "$base_dir/sb1.txt.gz" \
                              "$base_dir/sb2.txt.gz"    "$base_dir/metadata.csv" base
        else
            echo "----- Running recon-count.jl -----"
            mkdir fastqs
            gcloud storage cp "$fastq_dir/~{index}_*_L00[~{lanes}]_*" fastqs
            ls -1 fastqs
            julia --threads 8 recon-count.jl fastqs base -x ~{bc1} -y ~{bc2} -r '_L00[~{lanes}]_' -p ~{downsample_prop}
            rm -rf fastqs
            # Uploaded before knn.py gets a chance to fail, so a later OOM does not
            # discard hours of Julia work
            echo "----- Uploading recon-count.jl output -----"
            gcloud storage cp base/* "$base_dir/"
        fi

        if [ -n "$SELECTION" ]; then
            echo "----- Running puck-select.py + subset-matrix.py -----"
            gsutil -q stat "$work_dir/selection.json" || {
                echo "ERROR: no $work_dir/selection.json - run tools/puck-select.py first"; exit 1; }
            gcloud storage cp "$work_dir/selection.json" base/

            # Replay the circle recorded by tools/puck-select.py
            read -r SOURCE_PUCK CX CY RAD INV < <($PY -c "
import json
d = json.load(open('base/selection.json'))
print(d['source_puck'], d['center'][0], d['center'][1], d['radius'],
      '--invert' if d.get('invert') else '')")

            gcloud storage cp "$base_dir/$SOURCE_PUCK/Puck.csv" base/Puck.csv
            $PY puck-select.py base/Puck.csv -o base/Puck-selected.csv --non_interactive \
                --center "$CX" "$CY" --radius "$RAD" $INV
            $PY subset-matrix.py -i base -o cache -p base/Puck-selected.csv -b 2
            cp base/Puck-selected.csv base/selection.json cache/
            gcloud storage cp cache/* "$work_dir/"
        else
            mv base/* cache/  # work_dir == base_dir, so this is already uploaded
        fi

        rm -rf base  # free the disk the base matrix held before knn.py allocates
    fi

    ##### Stage 2: the KNN graph #####

    # knn2.npz is only reusable if the matrix underneath it was reused too
    if [ "$FRESH_MATRIX" = false ] && gsutil -q stat "$work_dir/knn2.npz"; then
        echo "----- Using cached knn.py output -----"
        gcloud storage cp "$work_dir/knn2.npz" cache
    else
        echo "----- Running knn.py -----"
        $PY knn.py -i cache -o cache -b 2 -k 2
        gcloud storage cp cache/knn2.npz "$work_dir/"
    fi
    ls -1 cache

    ##### Stage 3: the embedding #####

    # estimate_diameter() bins on total bead count, so a large cut silently rescales the puck.
    # subset-matrix.py records the pre-cut bin; carry it over unless params already sets one.
    DIAM=""
    if [ -n "$SELECTION" ] && ! grep -q -- '-D' <<< "~{params}"; then
        prev=$(grep -m1 '^subset_diameter_prev,' cache/metadata.csv | cut -d',' -f2 || true)
        new=$( grep -m1 '^subset_diameter_new,'  cache/metadata.csv | cut -d',' -f2 || true)
        if [ -n "$prev" ] && [ "$prev" != "None" ] && [ "$prev" != "$new" ]; then
            DIAM="-D $prev"
            echo "NOTE: subsetting changed the estimated diameter ${prev} -> ${new}um, preserving ${prev}um"
        fi
    fi

    echo "----- Running recon.py -----"
    $PY recon.py -i cache -o output -b 2 $DIAM ~{params}

    echo "----- Uploading results -----"
    gcloud storage cp -r output/* "$work_dir/"

    echo "==================== END RECONSTRUCTION ===================="

    >>>
    runtime {
        cpu: 8
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
        Int mem_GB
        Int disk_GB
        Int bc1 = 0
        Int bc2 = 0
        Int lanes = 12345678
        String params = ""
        String docker = "us-central1-docker.pkg.dev/velina-208320/terra/pipeline-image:latest"
        Float downsample_prop = 1.0
        String branch = "main"
        Int? pr
        String? subfolder
        String? selection
    }
    call recon {
        input:
            bcl = bcl,
            index = index,
            mem_GB = mem_GB,
            disk_GB = disk_GB,
            bc1 = bc1,
            bc2 = bc2,
            lanes = lanes,
            params = params,
            docker = docker,
            downsample_prop = downsample_prop,
            branch = branch,
            pr = pr,
            subfolder = subfolder,
            selection = selection
    }
}
