/* This program demonstrates how to generate pileup from multiple BAMs
 * simutaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -L. -lbam -lz
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include "bam.h"
#include "functions.h"

typedef struct {     // auxiliary data structure
  bamFile fp;      // the file handler
  bam_iter_t iter; // NULL if a region not specified
  int min_mapQ;    // mapQ filter
} aux_t;

void *bed_read(const char *fn); // read a BED or position list file
void bed_destroy(void *_h);     // destroy the BED data structure
int bed_overlap(const void *_h, const char *chr, int beg, int end); // test if chr:beg-end overlaps

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
  aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
  int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
  if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
  return ret;
}




#ifdef _MAIN_BAM2DEPTH
int main(int argc, char *argv[])
#else
int main_depth(int argc, char *argv[])
#endif
{
  int ii, n, tid, beg, end, pos, *n_plp;
  int baseQ = 0;
  int mapQ = 0;
  int maxDepth=INT_MAX; //I'd rather use a long here but needs to be int to fit into samtools. at least samtools use > so INT_MAX should always pass
  const bam_pileup1_t **plp;
  char *reg = 0; // specified region
  void *bed = 0; // BED data structure
  bam_header_t *h = 0; // BAM header of the 1st input
  aux_t **data;
  bam_mplp_t mplp;
  char strand='*';
  char usage[50000];
  sprintf(usage,"Usage: %.1000s [-r reg] [-q baseQthres] [-Q mapQthres] [-b in.bed] <in1.bam> [...]\n first and additional arguments: bam files to be parsed\n -r: region to get coverage for in samtools format e.g. chr1:1000-1029 (default: all positions in the reference)\n -b: bed file specifying multiple regions\n -q: only count positions with a quality greater than or equal this (default:0)\n -Q: only count reads with a map quality greater than or equal this (default:0) \n -d: approximate maximum depth counted for a base. In samtools version this is 8000. (default: INT_MAX)\n -h: (optional) display this message and exit\n",argv[0]);

  // parse the command line
  while ((n = getopt(argc, argv, "r:b:q:Q:d:h")) >= 0) {
    switch (n) {
      case 'r': reg = strdup(optarg); break;   // parsing a region requires a BAM header
      case 'b': bed = bed_read(optarg); break; // BED or position list file can be parsed
      case 'q': baseQ = atoi(optarg); break;   // base quality threshold
      case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
      case 'd': maxDepth = atoi(optarg); break;    // approximate max number of reads/base (in samtools version this is a hidden magic 8000)
      case 'h':
        fprintf(stdout,"%s",usage);
        return(0);
      case '?':
        if (optopt<33)
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
    }
  }
  if (optind == argc) {
    fprintf(stderr,"%s", usage);
    return 1;
  }

  // initialize the auxiliary data structures
  n = argc - optind; // the number of BAMs on the command line
  data = calloc(n, sizeof(void*)); // data[i] for the i-th input
  beg = 0; end = 1<<30; tid = -1;  // set the default region
  for (ii = 0; ii < n; ++ii) {
    bam_header_t *htmp;
    data[ii] = calloc(1, sizeof(aux_t));
    data[ii]->fp = bam_open(argv[optind+ii], "r"); // open BAM
    data[ii]->min_mapQ = mapQ;                    // set the mapQ filter
    htmp = bam_header_read(data[ii]->fp);         // read the BAM header
    if (ii == 0) {
      h = htmp; // keep the header of the 1st BAM
      if (reg){
			strand=parseRegionStrand(reg);
			bam_parse_region(h, reg, &tid, &beg, &end); // also parse the region
		}
    } else bam_header_destroy(htmp); // if not the 1st BAM, trash the header
    if (tid >= 0) { // if a region is specified and parsed successfully
      bam_index_t *idx = bam_index_load(argv[optind+ii]);  // load the index
      data[ii]->iter = bam_iter_query(idx, tid, beg, end); // set the iterator
      bam_index_destroy(idx); // the index is not needed any more; phase out of the memory
    }
  }

  // the core multi-pileup loop
  mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
  bam_mplp_set_maxcnt(mplp,maxDepth);
  n_plp = calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
  plp = calloc(n, sizeof(void*)); // plp[i] points to the array of covering reads (internal in mplp)
  while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
    if (pos < beg || pos >= end) continue; // out of range; skip
    if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
    fputs(h->target_name[tid], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
    for (ii = 0; ii < n; ++ii) { // base level filters have to go here
      int jj, m = 0;
      for (jj = 0; jj < n_plp[ii]; ++jj) {
        const bam_pileup1_t *p = plp[ii] + jj; // DON'T modfity plp[][] unless you really know
        if (p->is_del || p->is_refskip) ++m; // having dels or refskips at tid:pos
        else if (bam1_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality
		  else if (!checkStrand(strand,p->b->core.flag,1)) ++m;//wrong strand
      }
      printf("\t%d", n_plp[ii] - m); // this the depth to output
    }
    putchar('\n');
  }
  free(n_plp); free(plp);
  bam_mplp_destroy(mplp);

  bam_header_destroy(h);
  for (ii = 0; ii < n; ++ii) {
    bam_close(data[ii]->fp);
    if (data[ii]->iter) bam_iter_destroy(data[ii]->iter);
    free(data[ii]);
  }
  free(data); free(reg);
  if (bed) bed_destroy(bed);
  return 0;
}
