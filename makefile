SAMDIR=/home/shescott/installs/samtools
bamSplice : bamSplice.c logFactorial.c
	gcc -c logFactorial.c
	gcc -g -O2 -Wall -o bamSplice -D_MAIN_BAM2DEPTH bamSplice.c logFactorial.o -L$(SAMDIR) -lbam -lz -lm -lpthread -I$(SAMDIR) 

