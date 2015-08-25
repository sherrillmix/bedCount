/* This program demonstrates how to generate pileup from multiple BAMs * simutaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -L. -lbam -lz
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include "bam.h"

#define WORSTPVAL -DBL_MAX
#define NOCHROM 2
#define NOREADS 3
#define TOOSHORT 4

#define BED_READ_LENGTH 3000
#define DEBUG 0
//assuming chr name <200 characters
#define MAX_CHR_LENGTH 200


typedef struct {
	int32_t nNames;
	int32_t nBuffers;
	char** names;
} nameBuffer;

#define NEW_NAME_BUFFER(X) nameBuffer X = {.nNames=0, .nBuffers=0}

void addBuffersToNameBuffer(nameBuffer* buffer, int32_t n){
	int32_t newN;
	newN=buffer->nBuffers+n;
	char** newBuffer;
	int ii;
	newBuffer=(char **)malloc(newN*sizeof(char *));
	for(ii=0;ii<buffer->nNames;ii++)newBuffer[ii]=buffer->names[ii];
	if(buffer->nNames>0)free(buffer->names);
	buffer->names=newBuffer;
	buffer->nBuffers=newN;
}

void clearBuffer(nameBuffer* buffer){
	int ii;
	for(ii=0;ii<buffer->nNames;ii++)free(buffer->names[ii]);
	buffer->nNames=0;
}

void addNameToBuffer(nameBuffer* buffer, const char * name){
	if(buffer->nNames+1>buffer->nBuffers){
		addBuffersToNameBuffer(buffer,10000);
	}
	//don't need nNames+1 because 0 based array
	buffer->names[buffer->nNames]=(char *)malloc((strlen(name))*sizeof(char));
	strcpy(buffer->names[buffer->nNames],name);
	buffer->nNames++;
}

void copyBufferToBuffer(nameBuffer* buffer1,nameBuffer* buffer2){
	int ii;
	for(ii=0;ii<buffer1->nNames;ii++){
		addNameToBuffer(buffer2,buffer1->names[ii]);
	}
}

void destroyNameBuffer(nameBuffer* buffer){
	if(buffer->nBuffers>0){
		clearBuffer(buffer);
		free(buffer->names);
		buffer->nBuffers=0;
	}
}


int cstring_cmp(const void *a, const void *b){
	const char **ia = (const char **)a;
	const char **ib = (const char **)b;
	return strcmp(*ia, *ib);
}

int32_t countUniqueNamesInBuffer(nameBuffer* buffer){
	int ii;
	//don't need to do anything if 0 or 1 names
	if(buffer->nNames<2)return(buffer->nNames);
	//sort the names in buffer
	qsort(buffer->names,buffer->nNames,sizeof(char*),cstring_cmp);
	//count uniques in sorted buffer
	int32_t uniqCount=1;
	//for(ii=0;ii<buffer->nNames;ii++)printf("%d %s\n",ii,buffer->names[ii]);
	for(ii=1;ii<buffer->nNames;ii++){
		//printf("%s==%s\n",buffer->names[ii-1],buffer->names[ii]);
		if(strcmp(buffer->names[ii-1],buffer->names[ii])!=0)uniqCount++;
	}
	return(uniqCount);
}


typedef struct {     // auxiliary data structure
	bamFile fp;      // the file handler
	bam_iter_t iter; // NULL if a region not specified
	int min_mapQ;    // mapQ filter
} aux_t;

typedef struct {
	char chr[MAX_CHR_LENGTH];
	int start, tStart;
	int end;
	int32_t tid;
} region;

void setRegionStart(region *reg,int start){
	reg->start=start;
	reg->tStart=start-1;
}
int sprintRegion(char *str,const region *reg){
	sprintf(str,"%s:%d-%d",reg->chr,reg->start,reg->end);
	return(0);
}
char* printRegion(const region *reg){//careful with this one when parallel
	static char text[BED_READ_LENGTH];	
	sprintRegion(text,reg);
	return(text);
}


struct regionArray{
	region* regions;
	int num;
};
int destroyRegionArray(struct regionArray x){
	free(x.regions);
	return(0);
}

typedef struct {
	int *startStop[2];
	int num;
} startStopArray;
void destroyStartStopArray(startStopArray *startStops){
	int ii;
	for(ii=0;ii<2;ii++)free(startStops->startStop[ii]);
}

struct regionList{
	region reg;
	struct regionList* next;
};

typedef struct {
	char reg[300];//assuming chrx:98987898-78998789 less than 300 characters
	int *counts;
} regAndCounts;
typedef struct {
	char parent[300];
	int exonNum;
	int countNum;
	regAndCounts *exonCounts;
} exonCountArray;
int fprintExonCountArray(FILE *fp,exonCountArray *exons,int numGenes){
	int ii,jj,kk;
	for(ii=0;ii<numGenes;ii++){
		for(jj=0;jj<exons[ii].exonNum;jj++){
			fprintf(fp, "%s\t%s",exons[ii].parent,exons[ii].exonCounts[jj].reg);
			for(kk=0;kk<exons[ii].countNum;kk++){
				fprintf(fp, "\t%d",exons[ii].exonCounts[jj].counts[kk]);
			}
			fprintf(fp,"\n");
		}
	}
	return(0);
}



int getBedLine(FILE *fp,region* out){
	char buffer[BED_READ_LENGTH]; //set this outside so we don't have to create every time? 
	char tmp[BED_READ_LENGTH]; 
	if(fgets(buffer,BED_READ_LENGTH,fp) == NULL){
		if(feof(fp))return(1);
		fprintf(stderr,"Problem reading file %s\n",strerror(ferror(fp)));
		return(2);
	}

	//counters
	int ii;
	int tmpPos=0;
	int field=0;
	for(ii=0;ii<strlen(buffer);ii++){
		if(buffer[ii]=='\t'||buffer[ii]=='\n'){
			tmp[tmpPos]='\0';
			switch (field) {
				case 0: strcpy(out->chr,tmp); break;
				case 1: setRegionStart(out, atoi(tmp)+1); break; //+1 to deal with ucsc 0-based starts
				case 2: out->end = atoi(tmp); break;
			}

			field++;
			tmpPos=0;
			if(field>2)break;
		}else{
			tmp[tmpPos]=buffer[ii];
			tmpPos++;
		}
	}
	//make sure we get the whole line (discarding anything in other fields)
	while(buffer[strlen(buffer)-1]!='\n' && fgets(buffer,BED_READ_LENGTH,fp) != NULL){}
	return(0);
}

int readBed(FILE *bedFile, struct regionArray* out){
	region location;
	struct regionList *head=NULL, *current=NULL, *new;
	int count=0;
	int ii;
	while(getBedLine(bedFile,&location)==0){
		new=(struct regionList*)malloc(sizeof(struct regionList));
		new->reg=location;
		new->next=NULL;
		if(count==0) head=new;
		else current->next=new;
		current=new;
		count++;
	}
	if(head==NULL||count==0)return(1);
	out->regions=(region*)malloc(count*sizeof(region));
	out->num=count;
	current=head;
	for(ii=0;ii<count;ii++){
		out->regions[ii]=current->reg;
		new=current->next;
		free(current);
		current=new;
	}
	return(0);
}

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b){ // read level filters better go here to avoid pileup
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
	if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
	return ret;
}

int assignTidToLocations(struct regionArray *locations, bam_header_t *h){
	int anyAssigned=0;
	int dummy1, dummy2;
	int ii;
	for(ii=0;ii<locations->num;ii++){
		if(bam_parse_region(h, printRegion(&locations->regions[ii]), &locations->regions[ii].tid, &dummy1, &dummy2)!=0){
			if(DEBUG)fprintf(stderr,"Couldn't figure out %s. Called %d\n",printRegion(&locations->regions[ii]),locations->regions[ii].tid);
			//return(-1);
		}else{
			anyAssigned=1;
		}
	}
	return(!anyAssigned);	
}



//data structure to pass into bam_fetch
typedef struct {
	int *counter;
	struct regionArray regionArray;
	nameBuffer names;
	int onlyPaired;
} fetchData;
#define NEW_FETCH_DATA(X) fetchData X;X.names.nNames=0;X.names.nBuffers=0;X.onlyPaired=1;


// callback for bam_fetch()  
int fetchFunc(const bam1_t *b, void *data){  
	//convert data into regionArray
	fetchData *regionArrayAndCounter=(fetchData*)data;
	struct regionArray *region=&regionArrayAndCounter->regionArray;
	int onlyPaired=regionArrayAndCounter->onlyPaired;
	int ii,jj;
	uint32_t operation,length;
	uint32_t *cigarBuffer=bam1_cigar(b);
	int genomePos=(int)b->core.pos;
	int operationEnd;
	int isInsideTarget=0;
	//flag 256 => not primary alignment. samtools depth appears to ignore these by default so I guess I'll follow
	//flag 1 => read paired. Only keep pairs (could add option)
	//flag 2 => read mapped in proper pair. Only keep proper pairs
	if(b->core.flag & BAM_FSECONDARY || ( onlyPaired && (!(b->core.flag & BAM_FPAIRED) || !(b->core.flag & BAM_FPROPER_PAIR))))return(0);

	for(ii=0;ii<b->core.n_cigar;ii++){
		operation=((cigarBuffer[ii])&0xf);
		length=((cigarBuffer[ii])>>4);
		operationEnd=genomePos-1;
		if(operation==0||operation==2||operation==3||operation==7||operation==8){
			operationEnd+=length;
		}
		if(operation==6){fprintf(stderr,"SAM operations 6 padding undefined in this code but found here");exit(12);}

		if((operation==0||operation==8||operation==7)&&length>0){
			for(jj=0;jj<region->num;jj++){
				if(region->regions[jj].start <= operationEnd && region->regions[jj].end >= genomePos){
					//printf("%d/%d %d-%d:%d\n",ii+1,b->core.n_cigar,genomePos,operationEnd,operation);
					isInsideTarget=1;
					break;
				}
			}
			if(isInsideTarget)break;
		}
		genomePos=operationEnd+1;
	}
	if(isInsideTarget)*regionArrayAndCounter->counter=*regionArrayAndCounter->counter+1;
	if(isInsideTarget)addNameToBuffer(&regionArrayAndCounter->names,bam1_qname(b));
	//free(cigarBuffer); //I think bam1_cigar actually returns a pointer to something inside b and things crash if we free it. let bam take care of it
	return(0);  
}

int getCountsFromFiles(aux_t **data, bam_index_t **indices, int nFiles, region reg, int baseQ, int **counts){
	int base, anyBase;
	int gap;
	int pos;
	int ii, jj;
	int nBases=reg.end-reg.start+1;
	bam_mplp_t mplp;
	int *n_plp;
	const bam_pileup1_t **plp;
	n_plp = calloc(nFiles, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
	plp = calloc(nFiles, sizeof(void*)); // plp[i] points to the array of covering reads (internal in mplp)

	for(ii=0; ii<nFiles; ii++)data[ii]->iter = bam_iter_query(indices[ii], reg.tid, reg.tStart, reg.end); // set the iterator
	mplp = bam_mplp_init(nFiles, read_bam, (void**)data); // initialization
	bam_mplp_set_maxcnt(mplp,INT_MAX); //Otherwise reads/base is silently and magically capped at 8000 thanks to samtools. samtools evaluates with > so INT_MAX should be safe
	//loop through region
	base=0;
	anyBase=0;
	while (bam_mplp_auto(mplp, &reg.tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < reg.tStart || pos >= reg.end) continue; // out of range; skip //tStart is 0-based and end is 1-based
		if(base>=nBases){
			fprintf(stderr,"More bases than region size in %s\n",printRegion(&reg));
			exit(4);
		}
		for (ii = 0; ii < nFiles; ++ii) { // base level filters have to go here
			gap=0;
			for (jj = 0; jj < n_plp[ii]; ++jj) {
				const bam_pileup1_t *p = plp[ii] + jj; // DON'T modfity plp[][] unless you really know
				if (p->is_del || p->is_refskip) ++gap; // having dels or refskips at tid:pos
				else if (bam1_qual(p->b)[p->qpos] < baseQ) ++gap; // low base quality
			}
			counts[ii][pos-reg.tStart]=n_plp[ii] - gap; 
			if(!anyBase && counts[ii][base]>0)anyBase=1;
		}
		//printf("Base %d/%d of %d\n",base,pos,nBases);
		base++;
	}
	//skips uncovered positions
	//if(base!=nBases){
		//fprintf(stderr,"Less bases (%d) than region size (%d) in %s\n",base,nBases,printRegion(&reg));
		//exit(6);
	//}
	bam_mplp_destroy(mplp);
	free(n_plp);
	free(plp);
	for(ii=0; ii<nFiles; ii++)if(data[ii]->iter)bam_iter_destroy(data[ii]->iter);
	return(!anyBase);//0 if any reads found
}

int getUniqueReadsFromFiles(aux_t **data, bam_index_t **indices, const startStopArray *breaks, const region *location, int nFiles, int **exonCounts, int onlyPaired, nameBuffer *nameStores,int reportGlobalUnique){
	int ii,jj;
	//fetchData regionArrayAndCounter;
	NEW_FETCH_DATA(regionArrayAndCounter);

	regionArrayAndCounter.regionArray.regions=(region*)malloc(sizeof(region));
	regionArrayAndCounter.regionArray.regions[0].tid=location->tid;
	regionArrayAndCounter.regionArray.num=1;
	regionArrayAndCounter.onlyPaired=onlyPaired;
	strcpy(regionArrayAndCounter.regionArray.regions[0].chr,location->chr);

	for(ii=0;ii<breaks->num;ii++){ 
		if(breaks->startStop[1][ii]<breaks->startStop[0][ii]){
			fprintf(stderr,"Start greater than end when finding unique reads in breaks:%s\n",printRegion(location));
			exit(9);
		}
		setRegionStart(&regionArrayAndCounter.regionArray.regions[0],location->start+breaks->startStop[0][ii]);
		regionArrayAndCounter.regionArray.regions[0].end=location->start+breaks->startStop[1][ii];
		if(DEBUG)fprintf(stderr,"%d - %d,%s, tid:%d ",breaks->startStop[0][ii],breaks->startStop[1][ii],printRegion(&regionArrayAndCounter.regionArray.regions[0]),regionArrayAndCounter.regionArray.regions[0].tid);
		for(jj=0;jj<nFiles;jj++){
			regionArrayAndCounter.counter=&exonCounts[jj][ii];
			*regionArrayAndCounter.counter=0;
			bam_fetch(data[jj]->fp,indices[jj],regionArrayAndCounter.regionArray.regions[0].tid,regionArrayAndCounter.regionArray.regions[0].tStart,regionArrayAndCounter.regionArray.regions[0].end, &regionArrayAndCounter, fetchFunc);
			if(DEBUG)fprintf(stderr," Naive count%d: %d All names%d: %d Unique names%d: %d\n",jj,exonCounts[jj][ii],jj,regionArrayAndCounter.names.nNames,jj,countUniqueNamesInBuffer(&(regionArrayAndCounter.names)));
			//could put an option here or not calculate naive count
			exonCounts[jj][ii]=countUniqueNamesInBuffer(&(regionArrayAndCounter.names));
			//copy names to global buffer here
			if(reportGlobalUnique)copyBufferToBuffer(&(regionArrayAndCounter.names),&nameStores[jj]);
			//int kk;
			//for(kk=0;kk<regionArrayAndCounter.names.nNames;kk++)printf("%s ",regionArrayAndCounter.names.names[kk]);
			clearBuffer(&(regionArrayAndCounter.names));
		}
		destroyNameBuffer(&regionArrayAndCounter.names);
		free(regionArrayAndCounter.regionArray.regions);
		if(DEBUG)fprintf(stderr,"\n");
	}
	return(0);
}

int storeExonCounts(exonCountArray *exonStoreCell,const startStopArray *breaks, int * const*exonCounts,int nFiles, const region *location){
	int ii,jj;
	sprintRegion(exonStoreCell->parent,location);
	exonStoreCell->countNum=nFiles;
	exonStoreCell->exonNum=breaks->num;
	if(exonStoreCell->exonNum==0)return(0);
	exonStoreCell->exonCounts=(regAndCounts *)malloc(exonStoreCell->exonNum*sizeof(regAndCounts));
	for(ii=0;ii<exonStoreCell->exonNum;ii++){
		exonStoreCell->exonCounts[ii].counts=(int *)malloc(nFiles*sizeof(int));
		sprintf(exonStoreCell->exonCounts[ii].reg,"%s:%d-%d",location->chr,location->start+breaks->startStop[0][ii],location->start+breaks->startStop[1][ii]);
		for(jj=0;jj<nFiles;jj++){
			exonStoreCell->exonCounts[ii].counts[jj]=exonCounts[jj][ii];
		}
	}
	return(0);
}
int destroyExonCountArray(exonCountArray *exonStore, int num){
	int ii,jj;
	for(ii=0;ii<num;ii++){
		if(exonStore[ii].exonNum==0)continue;
		for(jj=0;jj<exonStore[ii].exonNum;jj++){
			free(exonStore[ii].exonCounts[jj].counts);
		}
		free(exonStore[ii].exonCounts);
	}
	return(0);
}

double getCountForRegion(region location, bam_index_t **indices, aux_t **data, int nFiles, int breakPadding, exonCountArray *exonStoreCell, int onlyPaired, nameBuffer *nameStores,int reportGlobalUnique){
	int nBases;
	int ii;
	int isTooShort=0;
	startStopArray breaks;

	if(location.tid<0)return(NOCHROM);

	nBases=location.end-location.start+1;
	if(nBases<=2*breakPadding)isTooShort=1; //don't return here since we want to store zeros

	breaks.num=1;
	for(ii=0;ii<2;ii++)breaks.startStop[ii]=calloc(1,sizeof(int));
	breaks.startStop[0][0]=breakPadding;
	breaks.startStop[1][0]=nBases-breakPadding-1;

	int **exonCounts=calloc(nFiles,sizeof(int*));
	for(ii=0;ii<nFiles;ii++)exonCounts[ii]=calloc(breaks.num,sizeof(int));

	if(!isTooShort)getUniqueReadsFromFiles(data, indices, &breaks, &location, nFiles, exonCounts, onlyPaired, nameStores,reportGlobalUnique);
	
	//Store break counts for later output
	storeExonCounts(exonStoreCell,&breaks,exonCounts,nFiles,&location);

	for(ii=0;ii<nFiles;ii++)free(exonCounts[ii]);
	free(exonCounts);
	//Clean up break finding
	destroyStartStopArray(&breaks);

	//might want to return break coords and counts too 
	if(isTooShort)return(TOOSHORT);
	else return(0);
}


struct scoreArgs{
	struct regionArray locations;
	bam_index_t **indices;
	aux_t **data;
	int nFiles;
	int stepSize;
	int start;
	int breakPadding;
	int onlyPaired;
	exonCountArray *exonStore;
	nameBuffer *nameStore;
	int reportGlobalUnique;
	int vocal;
};

void* getScoreParallel(void *getScoreArgs){
	int ii;
	struct scoreArgs *args=(struct scoreArgs *)getScoreArgs;
	//if(DEBUG)fprintf(stderr,"Thread started. start: %d, stepsize:%d\n",args->start,args->stepSize);
	for(ii=args->start;ii<args->locations.num;ii+=args->stepSize){
		if(args->vocal){
			fprintf(stderr,".");
			fflush(stderr);
		}
		//if(DEBUG)fprintf(stderr,"Thread %d working on region %d (%s)\n",args->start,ii,printRegion(&args->locations.regions[ii]));
		getCountForRegion(args->locations.regions[ii], args->indices, args->data, args->nFiles, args->breakPadding, &args->exonStore[ii],args->onlyPaired,args->nameStore,args->reportGlobalUnique);
		//if(DEBUG)fprintf(stderr,"Thread %d finished region %d. Score: %.3f\n",args->start,ii,args->scores[ii]);
	}
	pthread_exit(NULL);
}

#ifdef _MAIN_BAM2DEPTH
int main(int argc, char *argv[])
#else
int main_depth(int argc, char *argv[])
#endif
{

	int ii, jj, nFiles, mapQ = 20;
	bam_header_t *h=0; // BAM header of the 1st input
	aux_t ***data;
	bam_index_t ***indices;

	char *bed = NULL;
	FILE *bedFile;  

	struct regionArray locations;
	pthread_t *threads;
	struct scoreArgs *args;

	int nThreads=1;
	int breakPadding=15;
	int onlyPaired=1;
	int *globalCounts;
	int reportGlobalUnique=0;
	int vocal=0;
	
	//storing read counts in exons
	exonCountArray *exonStore;
	nameBuffer **nameStore;

	// parse the command line
	while ((ii = getopt(argc, argv, "b:Q:B:t:sGv")) >= 0) {
		switch (ii) {
			case 'b': bed = strdup(optarg); break; // BED or position list file can be parsed now
			case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
			case 'B': breakPadding = atoi(optarg); break; //don't count reads only falling within this number of bases within break
			case 't': nThreads = atoi(optarg); break; //number of threads to use
			case 's': onlyPaired = 0; break; //only report good pairs
			case 'G': reportGlobalUnique = 1; break; //report the total unique reads in all regions
			case 'v': vocal = 1; break; //report progress to stderr
		}
	}

	if(bed==NULL){
		fprintf(stderr,"Must specify a bed file\n");
		return(3);
	}
	
	if (optind == argc) {
		fprintf(stderr, "Usage: bamSplice [-q baseQthres] [-Q mapQthres] -b in.bed -g 1122... <in1.bam> [...]\n");
		return(1);
	}
	nFiles = argc - optind; // the number of BAMs on the command lined

	// initialize the auxiliary data structures
	data = calloc(nThreads, sizeof(aux_t**)); // data[i] for the j-th thread
	for(ii=0;ii<nThreads;ii++)data[ii]=calloc(nFiles, sizeof(aux_t*)); // data[i] for the i-th input in j-th thread
	indices =calloc(nThreads,sizeof(bam_index_t**));
	for(ii=0;ii<nThreads;ii++)indices[ii] = calloc(nFiles,sizeof(bam_index_t*));
	for(jj = 0; jj < nThreads; jj++){ //may not need to duplicate all these
		for (ii = 0; ii < nFiles; ++ii) {
			bam_header_t *htmp;
			data[jj][ii] = calloc(1, sizeof(aux_t));
			data[jj][ii]->fp = bam_open(argv[optind+ii], "r"); // open BAM
			if(data[jj][ii]->fp == NULL){
				fprintf(stderr,"Problem reading file %s\n",argv[optind+ii]);
				return(5);
			}
			data[jj][ii]->min_mapQ = mapQ;                    // set the mapQ filter
			htmp = bam_header_read(data[jj][ii]->fp);         // read the BAM header
			if (ii == 0) {
				h = htmp; // keep the header of the 1st BAM
			} else bam_header_destroy(htmp); // if not the 1st BAM, trash the header
			indices[jj][ii] = bam_index_load(argv[optind+ii]);  // load the index
		}
	}


	//read bed file
	bedFile=fopen(bed,"rt");
	if(bedFile == NULL){
		fprintf(stderr,"Problem reading file %s\n",bed);
		return(2);
	}
	free(bed);
	readBed(bedFile,&locations);
	fclose(bedFile);
	if(locations.num<1){
		fprintf(stderr,"No lines in file %s\n",bed);
		return(3);
	}
	if(assignTidToLocations(&locations,h)!=0){
		fprintf(stderr,"No known chromosomes in file");
		return(4);
	}

	threads=(pthread_t *)malloc(nThreads*sizeof(pthread_t));
	args=(struct scoreArgs *)malloc(nThreads*sizeof(struct scoreArgs));
	exonStore=(exonCountArray *)malloc(locations.num*sizeof(exonCountArray));
	for(ii=0;ii<locations.num;ii++)exonStore[ii].exonNum=0;
	nameStore=(nameBuffer **)malloc(nThreads*sizeof(nameBuffer*));
	for(ii=0;ii<nThreads;ii++){
		nameStore[ii]=(nameBuffer *)malloc(nFiles*sizeof(nameBuffer));
		for(jj=0;jj<nFiles;jj++){
			nameStore[ii][jj].nNames=0;
			nameStore[ii][jj].nBuffers=0;
			if(reportGlobalUnique)addBuffersToNameBuffer(&nameStore[ii][jj],200000); //save a bunch of room
		}
	}

	for(ii=0;ii<nThreads;ii++){
		//set up args
		args[ii].locations=locations;
		args[ii].indices=indices[ii];
		args[ii].data=data[ii];
		args[ii].nFiles=nFiles;
		args[ii].breakPadding=breakPadding;
		args[ii].stepSize=nThreads;
		args[ii].start=ii;
		args[ii].exonStore=exonStore;
		args[ii].onlyPaired=onlyPaired;
		args[ii].nameStore=nameStore[ii];
		args[ii].reportGlobalUnique=reportGlobalUnique;
		args[ii].vocal=vocal;
		if(pthread_create(&threads[ii],NULL,getScoreParallel,&args[ii])){fprintf(stderr,"Couldn't create thread");exit(9);}
	}
	for(ii=0;ii<nThreads;ii++){
		if(pthread_join(threads[ii],NULL)){fprintf(stderr,"Couldn't join thread");exit(10);}
	}
	free(threads);
	fprintExonCountArray(stdout,exonStore,locations.num);
	if(reportGlobalUnique){
		//find global uniques
		globalCounts=(int *)malloc(nFiles*sizeof(int));
		fprintf(stdout, "GLOBAL\tGLOBAL");
		for(jj=0;jj<nFiles;jj++){
			for(ii=1;ii<nThreads;ii++){ //use name[ii][0] as base
				copyBufferToBuffer(&nameStore[ii][jj],&nameStore[0][jj]);
				destroyNameBuffer(&nameStore[ii][jj]);
			}
			globalCounts[jj]=countUniqueNamesInBuffer(&nameStore[0][jj]);
			fprintf(stdout, "\t%d",globalCounts[jj]);
			destroyNameBuffer(&nameStore[0][jj]);
		}
		fprintf(stdout, "\n");
		free(globalCounts);
	}
	
	for(ii=0;ii<nThreads;ii++)free(nameStore[ii]);
	free(nameStore);
	free(args);
	destroyExonCountArray(exonStore,locations.num);
	destroyRegionArray(locations);
	free(exonStore);

	for(jj=0;jj<nThreads;jj++){
		for (ii = 0; ii < nFiles; ++ii) {
			bam_close(data[jj][ii]->fp);
			free(data[jj][ii]);
			bam_index_destroy(indices[jj][ii]);
		}
		free(indices[jj]);
		free(data[jj]);
	}
	bam_header_destroy(h);
	free(indices);
	free(data);
	fprintf(stderr,"All done. Thanks\n----------------\n");
	return 0;
}
