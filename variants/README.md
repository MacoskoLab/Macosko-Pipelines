bamscrape
---------

Records:
* FLAG, RNAME, POS, MAPQ
* BAM tags: CB, UB, xf, RE, ts, pa
* Discrepancies between the matching interval and the reference ()
* Matching intervals (start, end)
* Insertions (pos, str)
* Deletions (pos, len)
* Reference skips (pos, len)
* Soft-clipped sequences (pos, str)
* Metadata: total number of reads, #reads without CB, #reads without UB, #reads without RNAME
Does not record:
* Unmapped reads
* BAM tags other than CB, UB, xf, RE, ts, pa
* QNAME, RNEXT, PNEXT, TLEN
* Base quality at locations where the reference was matched

summary
-------

For each (CB, UMI) tuple, find the dominant (RNAME, strand) and compute:
* Number of reads and average MAPQ
* Number of HQ and LQ reads for each SNV, as well as the total overlaps
* Union of matching intervals
Does not compute:
* Anything related to insertions, deletions, reference skips, and soft-clipped sequences
