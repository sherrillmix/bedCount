SAMDIR=/home/shescott/installs/samtools
all : bamSplice bamSplice2 bam2depth bedCount

bamSplice : bamSplice.c logFactorial.o
	gcc -g -O2 -Wall -o bamSplice -D_MAIN_BAM2DEPTH bamSplice.c logFactorial.o -L$(SAMDIR) -lbam -lz -lm -lpthread -I$(SAMDIR) 

logFactorial.o : logFactorial.c
	gcc -c logFactorial.c

bamSplice2 : bamSplice2.c logFactorial.o
	gcc -g -O2 -Wall -o bamSplice2 -D_MAIN_BAM2DEPTH bamSplice2.c logFactorial.o -L$(SAMDIR) -lbam -lz -lm -lpthread -I$(SAMDIR) 

bam2depth : bam2depth.c
	gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -L$(SAMDIR) -lbam -lz -I$(SAMDIR) 

bedCount : bedCount.c
	gcc -g -O2 -Wall -o bedCount -D_MAIN_BAM2DEPTH bedCount.c -L$(SAMDIR) -lbam -lz -lpthread -I$(SAMDIR) 


oldBedCount : oldBedCount.c
	gcc -g -O2 -Wall -o oldBedCount -D_MAIN_BAM2DEPTH oldBedCount.c -L$(SAMDIR) -lbam -lz -lpthread -I$(SAMDIR) 
