#!/bin/bash

set -e

samFile=$(tempfile).sam
bamFile=${samFile%.sam}.bam
sortBamFile=${samFile%.sam}_sort.bam
#qname flag rname pos mapq cigar mcrap mpos isize seq qual tags
echo "@HD	VN:1.0	SO:coordinate	
@SQ	SN:chr1	LN:1000
read01	99	chr1	10	40	10M	=	18	0	ACACACACAC	>>>>>>>>>>
read02	1	chr1	11	35	11M	=	1	0	ACACACACACA	>>>>>>>>>>>
read03	1	chr1	10	40	10M	=	1	0	ACACACACAC	>>>>>>>>>>
read01	147	chr1	19	40	10M	=	1	0	GGGGGTTTTT	>>>>>>>>>>
">$samFile
samtools/samtools view -bS $samFile >$bamFile 2>/dev/null
samtools/samtools sort -f $bamFile $sortBamFile
samtools/samtools index $sortBamFile

samFile2=$(tempfile).sam
bamFile2=${samFile2%.sam}.bam
sortBamFile2=${samFile2%.sam}_sort.bam
#qname flag rname pos mapq cigar mcrap mpos isize seq qual tags
echo "@HD	VN:1.0	SO:coordinate	
@SQ	SN:chr1	LN:1000
read01	99	chr1	10	40	10M	=	18	0	ACACACACAC	>>>>>>>>>>
read02	1	chr1	10	35	10M	=	1	0	ACACACACAC	>>>>>>>>>>
read03	1	chr1	50	40	10M	=	1	0	ACACACACAC	>>>>>>>>>>
read01	147	chr1	19	40	10M	=	1	0	GGGGGTTTTT	>>>>>>>>>>
">$samFile2
samtools/samtools view -bS $samFile2 >$bamFile2 2>/dev/null
samtools/samtools sort -f $bamFile2 $sortBamFile2
samtools/samtools index $sortBamFile2


bedFile=$(tempfile).bed
echo "chr1 19 20
chr1 24 25">$bedFile


echo "Testing bedCount"
./bedCount 2>/dev/null && { echo "Missing files did not fail"; exit 1; }
./bedCount -h |grep Usage >/dev/null|| { echo "-h did not generate usage"; exit 1; }
./bedCount -Z 2>/dev/null&& { echo "Weird arg did not fail"; exit 1; }
./bedCount -$'\05' 2>/dev/null&& { echo "Weird arg did not fail"; exit 1; }

echo "Testing bam2depth"
./bam2depth 2>/dev/null && { echo "Missing files did not fail"; exit 1; }
./bam2depth -h |grep Usage >/dev/null|| { echo "-h did not generate usage"; exit 1; }
./bam2depth -Z 2>/dev/null&& { echo "Weird arg did not fail"; exit 1; }
./bam2depth -$'\05' 2>/dev/null&& { echo "Weird arg did not fail"; exit 1; }
./bam2depth $sortBamFile |head -1|grep "chr1	10	2" >/dev/null|| { echo "Unexpected output for bamFile"; exit 1; }
./bam2depth $sortBamFile -r chr1:19-20|head -1|grep "chr1	19	4" >/dev/null|| { echo "Unexpected output for region"; exit 1; }
./bam2depth $sortBamFile -b $bedFile|head -1|grep "chr1	20	2" >/dev/null|| { echo "Unexpected start output for bedFile"; exit 1; }
./bam2depth $sortBamFile -b $bedFile|tail -1|grep "chr1	25	1" >/dev/null|| { echo "Unexpected end output for bedFile"; exit 1; }
./bam2depth $sortBamFile -r chr1:10-11 -Q 35|tail -1|grep "chr1	11	3" >/dev/null|| { echo ">35 map Q count messed up"; exit 1; }
./bam2depth $sortBamFile -r chr1:10-11 -Q 36|tail -1|grep "chr1	11	2" >/dev/null|| { echo ">36 map Q count messed up"; exit 1; }
./bam2depth $sortBamFile -r chr1:11-19 -d100|head -1|grep "chr1	11	3" >/dev/null|| { echo "Max depth count messed up"; exit 1; }
./bam2depth $sortBamFile $sortBamFile2 -r chr1:10-11 |tail -1|grep "chr1	11	3	2" >/dev/null|| { echo "Multiple depth count messed up"; exit 1; }
echo "All clear"
exit 0
