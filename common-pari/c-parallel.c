/*
  Command-line options and parallel search for the search programs.

    -a    search all splitting types (a, b, ...) at every prime, no prompts
    -jN   (or -j N) use N processes

  The standard searches (c?-stand) take only -jN (par_jobs_option).

  With -jN the workers walk the same loops in the same order and claim the
  values of A_n (one per iteration of the A_n loop) through a counter in
  shared memory, so each value is searched by exactly one worker.  The
  counters are added up at the end.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>
#include "c-target.h"

#if DEG==3
#define ALL_TYPES "ab"
#elif DEG==4
#define ALL_TYPES "abcd"
#elif DEG==5
#define ALL_TYPES "abcdef"
#elif DEG==6
#define ALL_TYPES "abcdefghij"
#elif DEG==7
#define ALL_TYPES "abcdefghijklmn"
#else
#define ALL_TYPES "abcdefghijklmnopqrstu"
#endif

typedef struct {
  long next;
  long polcount, postdisc, t2count, discfcount;
  long qf_tested, qf_rejected, qf_loops;
  double dmax;
} par_shared_t;

static par_shared_t *par_shared = NULL;
static long par_item = 0, par_want = 0;
static int  par_worker = -1;          /* -1: serial */
static char par_menu[20*32];

static void par_usage(const char *prog)
{
  printf("usage: %s [-a] [-jN] p1 p2 ...\n"
         "  -a   search all splitting types (%s) at every prime, no prompts\n"
         "  -jN  use N processes\n", prog, ALL_TYPES);
  exit(1);
}

/* Remove the options from argv; the remaining arguments (program name
   first) go to args, and their number is returned. */
int par_options(int argc, char **argv, char **args, int *all, int *njobs)
{
  int k, nargs = 0;
  *all = 0;
  *njobs = 1;
  args[nargs++] = argv[0];
  for(k=1; k<argc; k++) {
    if(strcmp(argv[k], "-a") == 0) *all = 1;
    else if(strncmp(argv[k], "-j", 2) == 0) {
      const char *v = argv[k][2] ? argv[k]+2 : (k+1 < argc ? argv[++k] : NULL);
      if(v == NULL || (*njobs = atoi(v)) < 1) par_usage(argv[0]);
    }
    else if(argv[k][0] == '-') par_usage(argv[0]);
    else if(nargs < 19) args[nargs++] = argv[k];
  }
  args[nargs] = NULL;
  if(*all && nargs < 2) {
    printf("-a needs the primes on the command line\n");
    exit(1);
  }
  return nargs;
}

/* For the standard (discriminant bound) searches: remove only -jN / -j N
   from argv, so that negative numbers stay ordinary arguments.  The
   remaining arguments (program name first) go to args, and their number
   is returned. */
int par_jobs_option(int argc, char **argv, char **args, int *njobs)
{
  int k, nargs = 0;
  *njobs = 1;
  args[nargs++] = argv[0];
  for(k=1; k<argc; k++) {
    if(strncmp(argv[k], "-j", 2) == 0) {
      const char *v = argv[k][2] ? argv[k]+2 : (k+1 < argc ? argv[++k] : NULL);
      if(v == NULL || (*njobs = atoi(v)) < 1) {
        printf("-j needs a positive number of processes\n");
        exit(1);
      }
    }
    else if(nargs < 19) args[nargs++] = argv[k];
  }
  args[nargs] = NULL;
  return nargs;
}

/* For -a: make the menu code read "all types" for every prime. */
void par_all_types(void)
{
  long i;
  par_menu[0] = '\0';
  for(i=0; i<totprimes; i++) { strcat(par_menu, ALL_TYPES); strcat(par_menu, "\n"); }
  menuin = fmemopen(par_menu, strlen(par_menu), "r");
  if(menuin == NULL) { perror("fmemopen"); exit(1); }
}

/* Returns 1 if the caller should run the search loops (serial run, or a
   worker), 0 in the parent of a parallel run, after all workers are done
   and their counters have been added up. */
int par_start(int njobs, FILE *afile)
{
  int w;
  setvbuf(afile, NULL, _IOLBF, 0);    /* whole lines, so workers do not interleave */
  if(njobs <= 1) return 1;
  par_shared = mmap(NULL, sizeof(par_shared_t), PROT_READ | PROT_WRITE,
                    MAP_SHARED | MAP_ANONYMOUS, -1, 0);
  if(par_shared == MAP_FAILED) { perror("mmap"); exit(1); }
  memset(par_shared, 0, sizeof(*par_shared));
  fflush(stdout); fflush(afile); fflush(logfile);
  for(w=0; w<njobs; w++) {
    pid_t pid = fork();
    if(pid < 0) { perror("fork"); exit(1); }
    if(pid == 0) {
      par_worker = w;
      if(w > 0) verbosity = 0;        /* only the first worker reports progress */
      par_want = __sync_fetch_and_add(&par_shared->next, 1);
      return 1;
    }
  }
  while(wait(NULL) > 0) ;
  polcount = par_shared->polcount;   postdisc = par_shared->postdisc;
  t2count = par_shared->t2count;     discfcount = par_shared->discfcount;
  qf_tested = par_shared->qf_tested; qf_rejected = par_shared->qf_rejected;
  qf_loops = par_shared->qf_loops;   Dmaxroot = par_shared->dmax;
  return 0;
}

/* 1 if this process should search the current value of A_n. */
int par_claim(void)
{
  if(par_shared == NULL) return 1;
  if(par_item++ == par_want) {
    par_want = __sync_fetch_and_add(&par_shared->next, 1);
    return 1;
  }
  return 0;
}

/* End of the search loops: a worker reports its counters and exits. */
void par_finish(FILE *afile)
{
  if(par_worker < 0) return;
  fflush(afile);
  fflush(stdout);
  __sync_fetch_and_add(&par_shared->polcount, polcount);
  __sync_fetch_and_add(&par_shared->postdisc, postdisc);
  __sync_fetch_and_add(&par_shared->t2count, t2count);
  __sync_fetch_and_add(&par_shared->discfcount, discfcount);
  __sync_fetch_and_add(&par_shared->qf_tested, qf_tested);
  __sync_fetch_and_add(&par_shared->qf_rejected, qf_rejected);
  __sync_fetch_and_add(&par_shared->qf_loops, qf_loops);
  if(par_worker == 0) par_shared->dmax = Dmaxroot;
  _exit(0);
}
