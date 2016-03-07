#!/bin/bash

set -e

samFile=$(tempfile).sam
bamFile=${samFile%.sam}.bam
sortBamFile=${samFile%.sam}_sort.bam
#qname flag rname pos mapq cigar mcrap mpos isize seq qual tags
echo "@HD	VN:1.0	SO:coordinate	
@SQ	SN:chr1	LN:1000
read01	99	chr1	10	40	10M	=	18	0	ACACACACAC	>>>>>>>>>>
read02	1	chr1	10	40	10M	=	1	0	ACACACACAC	>>>>>>>>>>
read03	1	chr1	10	40	10M	=	1	0	ACACACACAC	>>>>>>>>>>
read01	147	chr1	19	40	10M	=	1	0	GGGGGTTTTT	>>>>>>>>>>
">$samFile
samtools/samtools view -bS $samFile >$bamFile 2>/dev/null
samtools sort -f $bamFile ${sortBamFile}
samtools index $sortBamFile


./bedCount 2>/dev/null && { echo "Missing files did not fail"; exit 1; }
./bedCount -h |grep Usage >/dev/null|| { echo "-h did not generate usage"; exit 1; }
./bedCount -Z 2>/dev/null&& { echo "Weird arg did not fail"; exit 1; }
./bedCount -$'\05' 2>/dev/null&& { echo "Weird arg did not fail"; exit 1; }

./bam2depth 2>/dev/null && { echo "Missing files did not fail"; exit 1; }
./bam2depth -h |grep Usage >/dev/null|| { echo "-h did not generate usage"; exit 1; }
./bam2depth -Z 2>/dev/null&& { echo "Weird arg did not fail"; exit 1; }
./bam2depth -$'\05' 2>/dev/null&& { echo "Weird arg did not fail"; exit 1; }
./bam2depth $bamFile |head -1|grep "chr1	10	3" >/dev/null|| { echo "Unexpected output for bamFile"; exit 1; }
./bam2depth $sortBamFile -r chr1:19-20|head -1|grep "chr1	19	4" >/dev/null|| { echo "Unexpected output for region"; exit 1; }
exit 0
