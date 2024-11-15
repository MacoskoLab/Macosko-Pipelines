bamscrape
---------

For each line in the BAM that has (CB, UB, RNAME), this script records:
* FLAG, RNAME, POS, MAPQ
* BAM tags: CB, UB, xf, RE, ts, pa
* Alternate base + quality at discrepancies between the matching interval and reference
* Insertions (pos, str)
* Deletions (pos, len)
* Reference skips (pos, len)
* Matching intervals (start, end)
* Soft-clipped sequences (pos, str)
* Metadata: #reads, #reads without CB, #reads without UB, #reads without RNAME

Does not record:
* Anything about unmapped reads
* BAM tags other than CB, UB, xf, RE, ts, pa
* QNAME, RNEXT, PNEXT, TLEN
* Base quality at locations where the reference was matched

summary
-------

For each (CB, UB) tuple, find the dominant (RNAME, strand) and compute:
* Number of reads and average MAPQ
* Number of HQ and LQ reads for each SNV, as well as the total number of covers
* Union of matching intervals
* Metadata: #reads chimeric, #umis chimeric

Does not compute:
* Anything related to insertions, deletions, reference skips, or soft-clipped sequences
* Anything related to non-CB/UB tags (xf, RE, ts, pa)
* Anything about the FLAG (except the strand)

umis.h5
-------

data
* umi - umi number (used as database key)
* cb_i - whitelists.h5 cb index (0-indexed)
* ub_i - whitelists.h5 ub index (0-indexed)
* rname_i - whitelists.h5 rname index (0-indexed)
* strand - bit 0x10 ^ 0x80 of FLAG (0 normal, 1 reverse complemented)
* reads - number of BAM records collapsed into the umi
* mapq_avg - average MAPQ for all reads for the umi, rounded to nearest integer

snv
* umi - umi number (used as database key)
* pos - position of SNV
* alt - ALT base
* hq - high-quality observations (#reads where snv base quality > QUAL_CUTOFF)
* lq - low-quality observations (#reads where snv base quality ≤ QUAL_CUTOFF)
* total - total number of reads covering the snv position

match
* umi - umi number (used as database key)
* start - 0-indexed position of matching interval start (inclusive)
* end - 0-indexed position of range end (exclusive)

metadata
* reads_chimeric - number of reads not belonging to a dominant (RNAME, strand) for a (CB, UB)
* umis_chimeric - number of (CB, UB) tuples discarded for not having a dominant (RNAME, strand)
