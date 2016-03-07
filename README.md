# bamTools: Simple tools to efficiently parse bam files

## Introduction
Tools to analyze bam files.

## Installation
To install, [download the repository](https://github.com/sherrillmix/suffixr/archive/master.zip) (or git clone) and run make in the resulting directory:

```
wget https://github.com/sherrillmix/bedCount/archive/master.zip
unzip master.zip
cd bedCount-master
make
```

## Usage
### bedCount
```
Usage: ./bam2depth [-r reg] [-q baseQthres] [-Q mapQthres] [-b in.bed] <in1.bam> [...]
 first and additional arguments: bam files to be parsed
 -r: region to get coverage for in samtools format e.g. chr1:1000-1029 (default: all positions in the reference)
  -q: only count positions with a quality greater than or equal this (default:0)
 -Q: only count reads with a map quality greater than or equal this (default:0) 
 -d: approximate maximum depth counted for a base. In samtools version this is 8000. (default: INT_MAX)
 -h: (optional) display this message and exit
```

### bam2depth
```
Usage: ./bedCount [-r reg] [-q baseQthres] [-Q mapQthres] [-b in.bed] <in1.bam> [...]
  first and additional arguments: bam files to be parsed
 -Q: only count reads with a map quality greater than or equal this (default:0) 
 -B: don't count reads only falling within this number of bases of the borders of a region
 -t: number of threads to use
 -s: only report good pairs -G report the total unique reads combined over all the regions
 -v: increase verbosity
 -h: (optional) display this message and exit
```


