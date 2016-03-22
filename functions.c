#include "bam.h"

//remove last character if + - or * and return else return *
char parseRegionStrand(char *reg){
  int lastPos;
  char lastChar;
  lastPos=strlen(reg)-1;
  if(lastPos<0)return('*');
  lastChar=reg[lastPos];
  switch(lastChar){
    case '+':break;
    case '-':break;
    case '*':break;
    default:
      return('*');
  }
  reg[lastPos]='\0';
  return(lastChar);
}

//return 1 if read correct strand, otherwise 0
int checkStrand(const char strand,const char flag, const int singleEnd){
  if(strand=='*')return(1);
  if(strand == '-'){
    if(flag & BAM_FREVERSE && flag & BAM_FREAD1)return(1);
    if(flag & BAM_FREVERSE && singleEnd)return(1);
    if(!(flag & BAM_FREVERSE) && flag & BAM_FREAD2)return(1);
  }
  if(strand == '+'){
    if(flag & BAM_FREVERSE && flag & BAM_FREAD2)return(1);
    if(!(flag & BAM_FREVERSE) && flag & BAM_FREAD1)return(1);
    if(!(flag & BAM_FREVERSE) && singleEnd)return(1);
  }
  return(0);
}
