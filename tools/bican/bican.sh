#!/usr/bin/env bash
# BICAN Data Management Utility
#
# Defines a `bican` shell function for uploading BICAN data (gene expression,
# puck files, and FASTQ files) to Google Cloud Storage to run jobs on Terra.
# This file is self-contained: it can be sourced for interactive one-row-at-a-
# time use, or executed directly to batch-replay many rows from a file.
#
# INTERACTIVE USE (one row at a time):
#   1. Set required environment variables (see below)
#   2. source bican.sh
#   3. bican [--dryrun] <command> [<command>...]
#      where <command> is: gex, puck, fastq, or all
#      --dryrun: preview commands without executing them (can go anywhere)
#
# REQUIRED ENVIRONMENT VARIABLES:
#   bcl           - BCL flowcell identifier
#   rnaidx        - RNA index
#   puck          - Puck recon directory name
#   spidx         - Spatial index
#   spfastqprefix - Spatial FASTQ prefix
#   gex           - Gene expression h5ad file
#   gcp_bucket    - Google Cloud Storage bucket
#
# OPTIONAL ENVIRONMENT VARIABLES:
#   fastqpost     - Optional suffix for fastq directory (default: empty/blank)
#                   Example: "_s62_unsorted" to use fastq_s62_unsorted directory
#   spfastqpath   - Optional custom fastq directory (default: empty/blank)
#                   If empty, uses: /broad/bican_um1_mccarroll_storage/slide_tags/data/flowcells/${bcl}/fastq${fastqpost}
#                   If set, uses: ${spfastqpath}${fastqpost}  (fastqpost appended as suffix; include the trailing 'fastq' dir in the value)
#   gexpath       - Optional custom directory path for gene expression file (default: empty/blank)
#                   If empty, uses: /broad/bican_um1_mccarroll/slide_tags/data/reprocessing/h5ad_full_annotations/${gex}
#                   If set, uses: ${gexpath}/${gex}
#
# USAGE EXAMPLES (interactive):
#   bican gex
#   bican puck
#   bican fastq
#   bican all
#   bican gex --dryrun
#   bican --dryrun all
#
# BATCH-REPLAY USE (many rows at once):
#   bican.py writes <BCL>.txt as plain data: one block of `KEY=VALUE` lines
#   per matching spreadsheet row, blocks separated by a blank line. Run this
#   script directly (not sourced) to replay every row through `bican`:
#
#   USAGE:
#     ./bican.sh <file.txt> [--rnaindex <RNAIndex>] <command> [<command>...] [--dryrun]
#
#     --rnaindex <RNAIndex>  Only run the row/chunk whose rnaidx matches
#                            (default: run every row in the file)
#
#   EXAMPLE WORKFLOW:
#     python bican.py 232HW5LT4 all                  # writes 232HW5LT4.txt
#     ./bican.sh 232HW5LT4.txt all --dryrun           # preview every matching row
#     ./bican.sh 232HW5LT4.txt all                    # actually run every matching row
#     ./bican.sh 232HW5LT4.txt --rnaindex SI-TT-B2 all --dryrun  # preview just one rxn

# Accept as argument a string and search in location to match last two terms and output any folder that matches as string
puck_locator(){
	local puck_input="$1"

	# Validate input argument
	if [ -z "$puck_input" ]; then
		echo "Error: No puck identifier provided" >&2
		return 1
	fi

	# Extract last two underscore-separated terms
	# e.g., from "Puck_20250723_Q16_JM" get "Q16_JM"
	local last_two=$(echo "$puck_input" | awk -F'_' '{print $(NF-1)"_"$NF}')

	# Search directory
	local recon_dir="/broad/bican_um1_mccarroll/slide_tags/data/recon"

	# Check if recon directory exists
	if [ ! -d "$recon_dir" ]; then
		echo "Error: Directory not found: $recon_dir" >&2
		return 1
	fi

	# Find matching directories and store in array (relative paths)
	local matches=()
	while IFS= read -r dir; do
		# Get relative path from recon_dir
		local rel_path="${dir#$recon_dir/}"
		matches+=("$rel_path")
	done < <(find "$recon_dir" -maxdepth 2 -type d -name "*_${last_two}")

	# Check number of matches
	local match_count=${#matches[@]}

	if [ $match_count -eq 0 ]; then
		echo "Error: No matching directories found for pattern '*_${last_two}'" >&2
		return 1
	elif [ $match_count -gt 1 ]; then
		echo "Multiple directories found with last two terms ('*_${last_two}'):" >&2
		printf '  %s\n' "${matches[@]}" >&2
		echo "Trying with last three terms..." >&2

		# Extract last three underscore-separated terms
		# e.g., from "Puck_20250723_Q16_JM" get "20250723_Q16_JM"
		local last_three=$(echo "$puck_input" | awk -F'_' '{print $(NF-2)"_"$(NF-1)"_"$NF}')

		# Find matching directories with last three terms (relative paths)
		local matches_three=()
		while IFS= read -r dir; do
			# Get relative path from recon_dir
			local rel_path="${dir#$recon_dir/}"
			matches_three+=("$rel_path")
		done < <(find "$recon_dir" -maxdepth 2 -type d -name "*_${last_three}")

		local match_count_three=${#matches_three[@]}

		if [ $match_count_three -eq 0 ]; then
			echo "Error: No matching directories found for pattern '*_${last_three}'" >&2
			return 1
		elif [ $match_count_three -gt 1 ]; then
			echo "Error: Multiple directories still found matching pattern '*_${last_three}':" >&2
			printf '  %s\n' "${matches_three[@]}" >&2
			return 1
		fi

		# Exactly one match with three terms - output it
		echo "${matches_three[0]}"
		return 0
	fi

	# Exactly one match - output it
	echo "${matches[0]}"
}

bican(){
	# Parse arguments for --dryrun flag and commands
	# --dryrun anywhere in the command applies to ALL commands
	local cmd_args=()
	local global_dryrun=false

	# Parse arguments: check if --dryrun appears anywhere
	for arg in "$@"; do
		if [ "$arg" = "--dryrun" ]; then
			global_dryrun=true
		elif [ "$arg" != "--dryrun" ]; then
			cmd_args+=("$arg")
		fi
	done

	# If "all" is specified, expand it to all three commands
	local final_cmd_args=()
	for cmd in "${cmd_args[@]}"; do
		if [ "$cmd" = "all" ]; then
			final_cmd_args+=("fastq" "gex" "puck")
		else
			final_cmd_args+=("$cmd")
		fi
	done

	# List of required environment variables
	local required_vars=("gcp_bucket" "bcl" "rnaidx" "puck" "spidx" "spfastqprefix" "gex")

	echo "Checking required environment variables..."

	# Check and print environment variables. Exit function if any are NULL.
	for var_name in "${required_vars[@]}"; do
		local var_value="${!var_name}" # Get the value of the variable named by var_name

		# Print the variable and its value, or "NULL"
		printf "%-15s: %s\n" "$var_name" "${var_value:-NULL}"

		# Check if the variable is null or unset
		if [ -z "$var_value" ]; then
			echo "Error: Environment variable '$var_name' is NULL or unset. Exiting function."
			return 1 # Use 'return' to exit the function with a non-zero status
		fi
	done

	# Print optional variables
	echo ""
	echo "Optional environment variables:"
	printf "%-15s: %s\n" "fastqpost" "${fastqpost:-(empty)}"
	printf "%-15s: %s\n" "spfastqpath" "${spfastqpath:-(empty)}"
	printf "%-15s: %s\n" "gexpath" "${gexpath:-(empty)}"

	# Resolve puck directory if puck command is requested OR if no commands (validation mode)
	local puck_dir=""
	local needs_puck=false

	# Check if puck command is in the list
	for cmd in "${final_cmd_args[@]}"; do
		if [ "$cmd" = "puck" ]; then
			needs_puck=true
			break
		fi
	done

	# Also resolve puck if no commands provided (validation/info mode)
	if [ ${#final_cmd_args[@]} -eq 0 ]; then
		needs_puck=true
	fi

	if [ "$needs_puck" = true ]; then
		echo ""
		echo "Resolving puck directory(ies)..."
		printf "%-15s: %s\n" "puck (input)" "$puck"
		local puck_dirs=()
		IFS=',' read -ra puck_list <<< "$puck"
		for p in "${puck_list[@]}"; do
			p=$(echo "$p" | xargs)  # trim whitespace
			local resolved
			resolved=$(puck_locator "$p")
			if [ $? -ne 0 ] || [ -z "$resolved" ]; then
				echo "Error: Failed to resolve puck directory for '$p'. Exiting function."
				return 1
			fi
			puck_dirs+=("$resolved")
		done
		local puck_dirs_joined=$(IFS=', '; echo "${puck_dirs[*]}")
		printf "%-15s: %s\n" "puck (resolved)" "\"$puck_dirs_joined\""
	fi

	if [ -z "$spfastqpath" ]; then
		fqdir="/broad/bican_um1_mccarroll_storage/slide_tags/data/flowcells/$bcl/fastq${fastqpost}"
	else
		fqdir="${spfastqpath}${fastqpost}"
	fi
	desc=slide-tags-bican-${sample}-${rnaidx}-${bcl}
	fq_opdir="${gcp_bucket}/fastqs/${bcl}"

	printf "\n"
	printf "%-15s: %s\n" "fqdir" "${fqdir}"
	printf "%-15s: %s\n" "fq pattern" "*${spfastqprefix}*"
	printf "%-15s: %s\n" "desc" "${desc}"
	printf "%-15s: %s\n" "fq_opdir" "${fq_opdir}"
	if [ "$needs_puck" = true ]; then
		local puck_basenames=()
		for pd in "${puck_dirs[@]}"; do
			puck_basenames+=("$(basename "$pd")")
		done
		local puck_basenames_joined=$(IFS=', '; echo "${puck_basenames[*]}")
		printf "%-15s: %s\n" "puck_opdir" "\"$puck_basenames_joined\""
	fi
	printf "\n"

	# Display dry run mode if enabled
	if [ "$global_dryrun" = true ]; then
		echo "**DRY RUN MODE** - Commands will be displayed but not executed"
	fi

	# If no commands provided, just print info and exit
	if [ ${#final_cmd_args[@]} -eq 0 ]; then
		echo "No commands specified. Environment validated and paths displayed above."
		return 0
	fi

	# Define commands
	# Determine gex source path based on gexpath variable
	local gex_source_path
	if [ -z "$gexpath" ]; then
		gex_source_path="/broad/bican_um1_mccarroll/slide_tags/data/reprocessing/h5ad_full_annotations/${gex}"
	else
		gex_source_path="${gexpath}/${gex}"
	fi
	local cmdgex="gcloud storage cp ${gex_source_path} ${gcp_bucket}/gene-expression/${bcl}/${rnaidx}/outs/"
	# cmdpuck is built dynamically per puck_dir in the puck command block
	local cmdfastq='for file in ${fqdir}/*${spfastqprefix}*;do gcloud storage cp $file $fq_opdir/${spidx}_$(basename $file); done'
	local cmdfastq_dryrun='for file in ${fqdir}/*${spfastqprefix}*;do echo gcloud storage cp $file $fq_opdir/${spidx}_$(basename $file); done'

	echo "Environment variables checked successfully. Proceeding with commands."
	printf "%-15s: %s\n" "Commands" "${final_cmd_args[*]}"
	printf "\n"

	# Process each command
	for cmd_arg in "${final_cmd_args[@]}"; do
		if [ "$cmd_arg" = "gex" ]; then
			if [ "$global_dryrun" = true ]; then
				printf "%-15s: %s\n" "Running (DRY)" "gex"
			else
				printf "%-15s: %s\n" "Running" "gex"
			fi
			echo ${cmdgex}
			if [ "$global_dryrun" = true ]; then
				# Check if gex file exists
				if [ -f "$gex_source_path" ]; then
					echo "✓ File exists: $gex_source_path"
				else
					echo "✗ File NOT found: $gex_source_path"
				fi
			else
				eval ${cmdgex}
				local exit_code=$?
				if [ $exit_code -eq 0 ]; then
					echo "✓ SUCCESS: Gene expression file uploaded"
				else
					echo "✗ FAILED: Gene expression upload failed with exit code $exit_code"
				fi
			fi
			printf "\n"
		elif [ "$cmd_arg" = "puck" ]; then
			if [ "$global_dryrun" = true ]; then
				printf "%-15s: %s\n" "Running (DRY)" "puck"
			else
				printf "%-15s: %s\n" "Running" "puck"
			fi
			for pd in "${puck_dirs[@]}"; do
				local cmdpuck="gcloud storage cp -r /broad/bican_um1_mccarroll/slide_tags/data/recon/${pd} ${gcp_bucket}/recon"
				echo ${cmdpuck}
				if [ "$global_dryrun" = true ]; then
					local puck_path="/broad/bican_um1_mccarroll/slide_tags/data/recon/${pd}"
					if [ -d "$puck_path" ]; then
						echo "✓ Directory exists: $puck_path"
					else
						echo "✗ Directory NOT found: $puck_path"
					fi
				else
					eval ${cmdpuck}
					local exit_code=$?
					if [ $exit_code -eq 0 ]; then
						echo "✓ SUCCESS: Puck directory uploaded (${pd})"
					else
						echo "✗ FAILED: Puck upload failed with exit code $exit_code (${pd})"
					fi
				fi
			done
			printf "\n"
		elif [ "$cmd_arg" = "fastq" ]; then
			if [ "$global_dryrun" = true ]; then
				printf "%-15s: %s\n" "Running (DRY)" "fastq"
			else
				printf "%-15s: %s\n" "Running" "fastq"
			fi
			eval ${cmdfastq_dryrun}
			if [ "$global_dryrun" = true ]; then
				# Check if fastq files exist and count them
				local file_count=0
				for file in ${fqdir}/*${spfastqprefix}*; do
					if [ -f "$file" ]; then
						file_count=$((file_count + 1))
					fi
				done
				if [ $file_count -gt 0 ]; then
					echo "✓ Found $file_count FASTQ file(s) matching: ${fqdir}/*${spfastqprefix}*"
				else
					echo "✗ No FASTQ files found matching: ${fqdir}/*${spfastqprefix}*"
				fi
			else
				eval ${cmdfastq}
				local exit_code=$?
				if [ $exit_code -eq 0 ]; then
					# Count files that were uploaded
					local file_count=0
					for file in ${fqdir}/*${spfastqprefix}*; do
						if [ -f "$file" ]; then
							file_count=$((file_count + 1))
						fi
					done
					echo "✓ SUCCESS: $file_count FASTQ file(s) uploaded"
				else
					echo "✗ FAILED: FASTQ upload failed with exit code $exit_code"
				fi
			fi
			printf "\n"
		else
			printf "%-15s: %s\n" "Unknown command" "$cmd_arg"
			echo "Warning: Command '$cmd_arg' does not match 'gex', 'puck', or 'fastq'. Skipping."
			printf "\n"
		fi
	done

	# Final summary
	if [ "$global_dryrun" = false ] && [ ${#final_cmd_args[@]} -gt 0 ]; then
		echo "=========================================="
		echo "All commands completed"
		echo "=========================================="
		echo "Summary:"
		printf "  %-15s: %s\n" "desc" "$desc"
		printf "  %-15s: %s\n" "bcl" "$bcl"
		printf "  %-15s: %s\n" "sample" "${sample:-N/A}"
		printf "  %-15s: %s\n" "rnaidx" "$rnaidx"
		printf "  %-15s: %s\n" "spidx" "$spidx"
		if [ "$needs_puck" = true ]; then
			local puck_dirs_summary=$(IFS=','; echo "${puck_dirs[*]}")
			printf "  %-15s: %s\n" "puck" "$puck_dirs_summary"
		fi
		printf "  %-15s: %s\n" "fastqpost" "${fastqpost:-(empty)}"
		printf "  %-15s: %s\n" "spfastqpath" "${spfastqpath:-(empty)}"
		printf "  %-15s: %s\n" "Commands run" "${final_cmd_args[*]}"
		printf "  %-15s: %s\n" "Bucket" "$gcp_bucket"
		printf "  %-15s: %s\n" "GEX location" "${gcp_bucket}/gene-expression/${bcl}/${rnaidx}/outs/"
		if [ "$needs_puck" = true ]; then
			local puck_locs=()
			for pd in "${puck_dirs[@]}"; do
				puck_locs+=("https://console.cloud.google.com/storage/browser/${gcp_bucket#gs://}/recon/$(basename "$pd")")
			done
			local puck_locs_joined=$(IFS=','; echo "${puck_locs[*]}")
			printf "  %-15s: %s\n" "Puck location" "$puck_locs_joined"
		fi
		printf "  %-15s: %s\n" "FASTQ location" "${gcp_bucket}/fastqs/${bcl}"
		echo "=========================================="
	fi
}

# Batch-replay mode: only runs when this file is executed directly (./bican.sh ...),
# not when it's `source`d for interactive use — sourcing must only define the
# functions above, since `exit` here would otherwise kill the sourcing shell.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	# No -e: `bican`'s error handling (checking $? after eval/puck_locator calls)
	# assumes a non-strict shell, and one bad row shouldn't abort the whole batch.
	# No -u: bash <4.4 (e.g. RHEL's 4.2, used on the Broad cluster) treats
	# "${empty_array[@]}" as an unbound-variable error under nounset, and this
	# script expands several arrays (cmd_args, puck_dirs, matches, ...) that can
	# legitimately be empty.
	set -o pipefail

	if [ $# -lt 2 ]; then
		echo "Usage: $0 <file.txt> [--rnaindex <RNAIndex>] <command> [<command>...] [--dryrun]" >&2
		exit 1
	fi

	file="$1"; shift

	if [ ! -f "$file" ]; then
		echo "Error: file not found: $file" >&2
		exit 1
	fi

	rnaindex_filter=""
	cmd_args=()
	while [ $# -gt 0 ]; do
		case "$1" in
			--rnaindex)
				rnaindex_filter="$2"
				shift 2
				;;
			--rnaindex=*)
				rnaindex_filter="${1#--rnaindex=}"
				shift
				;;
			*)
				cmd_args+=("$1")
				shift
				;;
		esac
	done

	if [ ${#cmd_args[@]} -lt 1 ]; then
		echo "Usage: $0 <file.txt> [--rnaindex <RNAIndex>] <command> [<command>...] [--dryrun]" >&2
		exit 1
	fi

	row=0
	block=""

	run_block(){
		[ -z "$block" ] && return 0
		while IFS= read -r line; do
			[ -z "$line" ] && continue
			eval "export $line"
		done <<< "$block"
		if [ -n "$rnaindex_filter" ] && [ "$rnaidx" != "$rnaindex_filter" ]; then
			return 0
		fi
		row=$((row + 1))
		echo "=========================================="
		echo "Row $row (rnaidx=$rnaidx)"
		echo "=========================================="
		bican "$@"
	}

	while IFS= read -r line || [ -n "$line" ]; do
		if [ -z "$line" ]; then
			run_block "${cmd_args[@]}"
			block=""
		else
			block+="$line"$'\n'
		fi
	done < "$file"
	run_block "${cmd_args[@]}"

	if [ -n "$rnaindex_filter" ] && [ "$row" -eq 0 ]; then
		echo "Error: no row with rnaidx '$rnaindex_filter' found in $file" >&2
		exit 1
	fi

	echo "=========================================="
	echo "Processed $row row(s) from $file"
	echo "=========================================="
fi
