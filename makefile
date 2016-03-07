SAMDIR=samtools
all : bam2depth bedCount

samtools/libbam.a: samtools-0.1.19.tar.bz2
	tar xvfj samtools-0.1.19.tar.bz2
	mv samtools-0.1.19 samtools
	cd samtools &&  make

bam2depth : bam2depth.c samtools/libbam.a
	gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -L$(SAMDIR) -lbam -lz -I$(SAMDIR) -lpthread

bedCount : bedCount.c samtools/libbam.a
	gcc -g -O2 -Wall -o bedCount -D_MAIN_BAM2DEPTH bedCount.c -L$(SAMDIR) -lbam -lz -lpthread -I$(SAMDIR) 

