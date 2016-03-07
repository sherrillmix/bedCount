SAMDIR=/home/shescott/installs/samtools
all : bam2depth bedCount


bam2depth : bam2depth.c
	gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -L$(SAMDIR) -lbam -lz -I$(SAMDIR) 

bedCount : bedCount.c
	gcc -g -O2 -Wall -o bedCount -D_MAIN_BAM2DEPTH bedCount.c -L$(SAMDIR) -lbam -lz -lpthread -I$(SAMDIR) 

