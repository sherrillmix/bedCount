# bedCount: Simple tools to efficiently parse bam files
[![Build Status](https://travis-ci.org/sherrillmix/bedCount.svg?branch=master)](https://travis-ci.org/sherrillmix/bedCount)
[![codecov.io](https://codecov.io/github/sherrillmix/bedCount/coverage.svg?branch=master)](https://codecov.io/github/sherrillmix/bedCount?branch=master)

## Introduction
Some preliminary tools using the [samtools](http://www.htslib.org/) library to analyze bam files.

## Installation
To install, [download the repository](https://github.com/sherrillmix/bedCount/archive/master.zip) (or git clone) and run make in the resulting directory:

```
wget https://github.com/sherrillmix/bedCount/archive/master.zip
unzip master.zip
cd bedCount-master
make
```

## Usage
### bedCount
```
Usage: ./bedCount [-r reg] [-q baseQthres] [-Q mapQthres] [-b in.bed] [-s] [-G] [-v] [-h] <in1.bam> [...]
Output: streams to standard out a tab separated file with columns; region, trimmedRegion, count1, count2 ...
  where countX is is the count for that region in the Xth file argument
Arguments:
 first and additional arguments: bam files to be parsed
 -Q: only count reads with a map quality greater than or equal this (default:0) 
 -B: don't count reads only falling within this number of bases of the borders of a region (default: 15)
 -t: number of threads to use (default: 1)
 -s: Data is single end reads. Do not filter out reads not in a well mapped pair. Assume strandedness behaves as first read (default: paired reads)
 -G report the total unique reads combined over all the regions
 -v: increase verbosity
 -h: display this message and exit
```

### bam2depth
```
Usage: ./bam2depth [-r reg] [-q baseQthres] [-Q mapQthres] [-b in.bed] <in1.bam> [...]
Output: streams to standard out a tab separated file with columns; chromosome, position, count1, count2 ...
  where countX is the count for that position in the Xth file argument
Arguments:
 first and additional optional arguments: bam files to be parsed
 -r: region to get coverage for in samtools format e.g. chr1:1000-1029 (default: all positions in the reference)
 -b: bed file specifying multiple regions
 -q: only count positions with a quality greater than or equal this (default:0)
 -Q: only count reads with a map quality greater than or equal this (default:0) 
 -d: approximate maximum depth counted for a base. In samtools version this is 8000. (default: INT_MAX)
 -h: (optional) display this message and exit
```


