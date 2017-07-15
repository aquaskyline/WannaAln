#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#include "edlib.h"
#include "revcomp.h"

#define DEBUG
#ifdef DEBUG
#define VERBOSE(...) fprintf(stderr, __VA_ARGS__)
#else
#define VERBOSE(...)
#endif

#define CSTRLEN 1024

void help()
{
  fprintf(stderr, "wannaAln\n");
  fprintf(stderr, "Ruibang Luo\n");
  fprintf(stderr, "rluo5@jhu.edu\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "-a [STR]   Input R1\n");
  fprintf(stderr, "-b [STR]   Input R2\n");
  fprintf(stderr, "-x [STR]   Output R1, gzipped\n");
  fprintf(stderr, "-y [STR]   Output R2, gzipped\n");
  fprintf(stderr, "-p [STR]   Matching pattern\n");
  fprintf(stderr, "-m [INT]   Edit distance allowed, default 1\n");
  fprintf(stderr, "-h         Show help information\n");
  fprintf(stderr, "\n");
}

int main(int argc, char** argv)
{
  int c, flag, ed, ql;
  flag = 0;
  ed = 1; // 1 edit distance allowed by default
  char q[CSTRLEN], rcq[CSTRLEN]; // Alignment query pattern
  char *r1ifn, *r2ifn; // Input file names
  char *r1ofn, *r2ofn; // Output file names

  while ((c = getopt (argc, argv, "a:b:x:y:q:m:h")) != -1)
  {
    switch(c)
    {
      case 'a':
        flag |= 0x1;
        r1ifn = optarg;
        VERBOSE("Input R1: %s\n", r1ifn);
        break;
      case 'b':
        flag |= 0x2;
        r2ifn = optarg;
        VERBOSE("Input R2: %s\n", r2ifn);
        break;
      case 'x':
        flag |= 0x4;
        r1ofn = optarg;
        VERBOSE("Output R1: %s\n", r1ofn);
        break;
      case 'y':
        flag |= 0x8;
        r2ofn = optarg;
        VERBOSE("Output R2: %s\n", r2ofn);
        break;
      case 'q':
        flag |= 0x10;
        strcpy(q, optarg);
        VERBOSE("Query pattern: %s\n", q);
        ql = strlen(q);
        RevComp(q, rcq);
        VERBOSE("Query pattern reverse complement: %s\n", rcq);
        break;
      case 'm':
        ed = atoi(optarg);
        VERBOSE("Edit distance: %d\n", ed);
        break;
      case 'h':
        help();
        exit(0);
      case '?':
        if (optopt == 'a' || optopt == 'b' || optopt == 'x' || optopt == 'y' || optopt == 'q')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
      default:
        fprintf(stderr, "Error parsing parameters. Exiting...\n"); exit(255);
    }
  }
  if(flag != 0x1F)
  { help(); fprintf(stderr, "Parameters a, b, x, y and q must be set. Exiting...\n"); exit(255); }

  gzFile r1ifp, r2ifp, r1ofp, r2ofp;
  kseq_t *r1seq, *r2seq;
  int n;
  n = 0;
  r1ifp = gzopen(r1ifn, "r");
  r2ifp = gzopen(r2ifn, "r");
  if(r1ifp == Z_NULL || r2ifp == Z_NULL)
  { fprintf(stderr, "Error opening input files. %s\n"); strerror(errno); }
  r1ofp = gzopen(r1ofn, "w3");
  r2ofp = gzopen(r2ofn, "w3");
  if(r1ofp == Z_NULL || r2ofp == Z_NULL)
  { fprintf(stderr, "Error opening output files. %s\n"); strerror(errno); }
  r1seq = kseq_init(r1ifp);
  r2seq = kseq_init(r2ifp);

  EdlibAlignConfig alnCfg = edlibNewAlignConfig(ed, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE);

  while (kseq_read(r1seq) >= 0)
  {
    ++n;
    kseq_read(r2seq);
    EdlibAlignResult r1 = edlibAlign(q, ql, r1seq->seq.s, r1seq->seq.l, alnCfg);
    EdlibAlignResult r1rc = edlibAlign(rcq, ql, r1seq->seq.s, r1seq->seq.l , alnCfg);
    EdlibAlignResult r2 = edlibAlign(q, ql, r2seq->seq.s, r2seq->seq.l, alnCfg);
    EdlibAlignResult r2rc = edlibAlign(rcq, ql, r2seq->seq.s, r2seq->seq.l , alnCfg);
    //printf("1\t%d\t%s\t%s\t%s\t%d\t%d\n", n, r1seq->name.s, r1seq->seq.s, r1seq->qual.s, r1.editDistance, r1rc.editDistance);
    //printf("2\t%d\t%s\t%s\t%s\t%d\t%d\n", n, r2seq->name.s, r2seq->seq.s, r2seq->qual.s, r2.editDistance, r2rc.editDistance);
    if(r1.editDistance != -1 || r1rc.editDistance != -1 || r2.editDistance != -1 || r2rc.editDistance != -1)
    {
      gzprintf(r1ofp, "@%s\n%s\n+\n%s\n", r1seq->name.s, r1seq->seq.s, r1seq->qual.s);
      gzprintf(r2ofp, "@%s\n%s\n+\n%s\n", r2seq->name.s, r2seq->seq.s, r2seq->qual.s);
    }
    edlibFreeAlignResult(r1);
    edlibFreeAlignResult(r1rc);
    edlibFreeAlignResult(r2);
    edlibFreeAlignResult(r2rc);
    if(n % 100000 == 0)
    { VERBOSE("\b%d reads processed.", n); }
  }
  VERBOSE("\n");

  kseq_destroy(r1seq);
  kseq_destroy(r2seq);
  gzclose(r1ifp);
  gzclose(r2ifp);
  gzclose(r1ofp);
  gzclose(r2ofp);

  exit(0);
}
