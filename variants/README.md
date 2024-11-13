bamscrape
---------

For each read in the BAM, this script records:
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

For each (CB, UMI) tuple, find the dominant (RNAME, strand) and compute:
* Number of reads and average MAPQ
* Number of HQ and LQ reads for each SNV, as well as the total number of covers
* Union of matching intervals
* Metadata: #reads chimeric, #umis chimeric

Does not compute:
* Anything related to insertions, deletions, reference skips, or soft-clipped sequences
* Anything related to non-CB/UB tags (xf, RE, ts, pa)
* Anything about the FLAG (except the strand)
