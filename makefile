SAMDIR=samtools
all : bam2depth bedCount
.PHONY: all test
FLAGS=-g -O2 -Wall -L$(SAMDIR) -lbam -lz -I$(SAMDIR) -lpthread

samtools/libbam.a: samtools-0.1.19.tar.bz2
	tar xvfj samtools-0.1.19.tar.bz2
	mv samtools-0.1.19 samtools
	cd samtools &&  make

TEMP := $(shell mktemp)
TEMP2 := $(shell mktemp)
README.md: bam2depth bedCount README.template makefile
	./bedCount -h>$(TEMP)
	./bam2depth -h>$(TEMP2)
	sed -e "/##BEDCOUNT_USAGEHERE##/{r $(TEMP)" -e "d;}" -e "/##BAM2DEPTH_USAGEHERE##/{r $(TEMP2)" -e "d;}" <README.template >README.md

bam2depth : bam2depth.c samtools/libbam.a
	$(CC) -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c $(FLAGS)

bedCount : bedCount.c samtools/libbam.a
	$(CC) -o bedCount -D_MAIN_BAM2DEPTH bedCount.c $(FLAGS)

#testBed.c 
test: bam2depth.c bedCount.c tests.bash makefile
	$(CC) -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c $(FLAGS) -coverage
	$(CC) -o bedCount -D_MAIN_BAM2DEPTH bedCount.c $(FLAGS) -coverage
	#$(CC) -Wall -o testBed testBed.c -lz -lpthread -coverage 
	#./testBed
	bash tests.bash
	#gcov main.c tree.h tree.c
	rm bam2depth bedCount
	make all
