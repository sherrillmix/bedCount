SAMDIR=/home/shescott/installs/samtools
all : bamSplice bamSplice2 

bamSplice : bamSplice.c logFactorial.o
	gcc -g -O2 -Wall -o bamSplice -D_MAIN_BAM2DEPTH bamSplice.c logFactorial.o -L$(SAMDIR) -lbam -lz -lm -lpthread -I$(SAMDIR) 

logFactorial.o : logFactorial.c
	gcc -c logFactorial.c

bamSplice2 : bamSplice2.c logFactorial.o
	gcc -g -O2 -Wall -o bamSplice2 -D_MAIN_BAM2DEPTH bamSplice2.c logFactorial.o -L$(SAMDIR) -lbam -lz -lm -lpthread -I$(SAMDIR) 
