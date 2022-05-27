# bam2bedpe

## Summary
- scripts for converting BAM into BEDPE format
- BAM file size prohibitive large
- BEDPE retains essential BAM information while having 6x size reduction
- bedtools bam2bedpe conversion cannot handle secondary alignments
- bam2bedpe script here preserves all alignments

## Workflow overview

![bam2bedpe](https://user-images.githubusercontent.com/98410560/170511327-ca6435bf-9d32-4676-890e-caf9f067f1ad.png)

## BEDPE7+12 file format
- inspired by bedtools definition (https://bedtools.readthedocs.io/en/latest/content/general-usage.html)
- originally designed to describe disjoint genome features (e.g. structural variants), also used for paired end reads which cannot be handled by traditional BED files
- first read in pair is always READ1, i.e. column order is READ1 READ2 (as oppose to coordinate order)
- follows 0-based coordinate system ( first base in chromosome is numbered 0 )
- chromosome name can be any character ( e.g. chr6_KQ090016v1_fix, chr10_GL383545v1_alt, chr11_KI270721v1_random, chrUn_GL000195v1, scaffold123, contig111 )
- TLEN is identical to ones found in BAM (calculated where length only include M CIGAR)

### BEDPE7+12 fields

| Col  | BED Field        | Type   | Regex or range           | Brief Description                | BAM Equivalent            |
| ---- | ---------------- | ------ | ------------------------ | -------------------------------- | ------------------------- |
| 1    | read1_chrom      | String | `[[:alnum:]_]{1,255}^4`  | read1 chromosome name            | RNAME                     |
| 2    | read1_chromStart | Int    | `.\|([0,264−1])`         | read1 start position             | POS                       |
| 3    | read1_chromEnd   | Int    | `.\|([0,264−1])`         | read1 end position               | None, calculated by pysam |
| 4    | read2_chrom      | String | `[[:alnum:]_]{1,255}^4`  | read2 chromosome name            | RNAME                     |
| 5    | read2_chromStart | Int    | `.\|([0,264−1])`         | read2 start position             | POS                       |
| 6    | read2_chromEnd   | Int    | `.\|([0,264−1])`         | read2 end position               | None, calculated by pysam |
| 7    | fragment_id      | String | `[\x20-\x7e]{1,255}`     | fragment id                      | QNAME                     |
| 8    | read1_MAPQ       | Int    | `[0, 2^31 - 1]`          | read1 mapping quality            | MAPQ                      |
| 9    | read2_MAPQ       | Int    | `[0, 2^31 - 1]`          | read2 mapping quality            | MAPQ                      |
| 10   | read1_STRAND     | String | `[-+.]`                  | read1 strand                     | None, calculated by pysam |
| 11   | read2_STRAND     | String | `[-+.]`                  | read2 strand                     | None, calculated by pysam |
| 12   | read1_CIGAR      | String | `*\|([0-9]+[MIDNSHPX=])+`| read1 CIGAR string               | CIGAR                     |
| 13   | read2_CIGAR      | String | `*\|([0-9]+[MIDNSHPX=])+`| read2 CIGAR string               | CIGAR                     |
| 14   | read1_FLAG       | Int    | `[0, 2^16 -1]`           | read1 bitwise FLAG               | FLAG                      |
| 15   | read2_FLAG       | Int    | `[0, 2^16 -1]`           | read2 bitwise FLAG               | FLAG                      |
| 16   | read1_TLEN       | Int    | `[-2^31 + 1, 2^31 -1]`   | read1 template length            | TLEN                      |
| 17   | read2_TLEN       | Int    | `[-2^31 + 1, 2^31 -1]`   | read2 template length            | TLEN                      |
| 18   | read1_NM_TAG     | Int    | `NA\|([0-9]*)`           | read1 edit distance to reference | NM                        |
| 19   | read2_NM_TAG     | Int    | `NA\|([0-9]*)`           | read2 edit distance to reference | NM                        |

Please see [wiki](https://github.com/mhanbioinfo/bam2bedpe/wiki) for tutorial on settings up and running module.

Questions please contact ming.han@uhn.ca.
