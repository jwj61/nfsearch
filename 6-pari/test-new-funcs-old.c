#include <stdio.h>
#include <stdlib.h>
#include "c-target.h"


long pdata[10][8][7];
long plist[10];
long pflags[10][8];

long long arc[10], shift[10];

int t6init(long, long, long);    long t6next(long, long);
int t51init(long, long, long);   long t51next(long, long);
int t42init(long, long, long);   long t42next(long, long);
int t411init(long, long, long);  long t411next(long, long);
int t33init(long, long, long);  long t33next(long, long);
int t321init(long, long, long);  long t321next(long, long);
int t3111init(long, long, long);  long t3111next(long, long);
int t222init(long, long, long);  long t222next(long, long);
int t2211init(long, long, long);  long t2211next(long, long);
int t21111init(long, long, long);  long t21111next(long, long);
int w_file_init(long, long, long);  long w_file_next(long, long);
long next_all_primes(long , long);
long init_all_primes(long , long);

long (*(ff_next_list[]))() = {t6next, t51next, t42next, t411next, t33next,
                             t321next, t3111next, t222next, t2211next,
                             t21111next, w_file_next};
int (*(ff_init_list[]))() = {t6init, t51init, t42init, t411init, t33init,
                            t321init, t3111init, t222init, t2211init,
                            t21111init, w_file_init};

long (*(pf_next_list[8]))();
int (*(pf_init_list[8]))();
long totprimes;

int main(int argc, char * argv []) {
  /*  int (*f_init)();
      long (*f_next)(); */
  long j,a1,a2,a3,a4,a5,a6, pc;
  long myprime;
  int myfunc;

  getprimelist(argc, argv);

  /* myfunc = atoi(argv[1]);
     myprime = atoi(argv[2]); */
  
  /*  plist[0] = myprime;
  plist[1] = 11;
  plist[2] = 5;
  plist[3] = 2003; */

  /*  totprimes = 1;
      plist[0] = myprime; */

  /*  for(j=0; j<totprimes; j++) {
    pf_next_list[j] = ff_next_list[myfunc];
    pf_init_list[j] = ff_init_list[myfunc];
  }
  */

  getsearchcases();

  /* Here is where we list the functions */
  /*  for(l=0; l<4;l++) {
    f_init = f_init_list[l]; 
    f_next = f_next_list[l];
    printf("\n---------------------------------------------\n"); */

  for(pc=0; pc<totprimes; pc++) {
    printf("\n\n p = %ld\n", plist[pc]);
    init_all_primes(0,0);
  
    init_all_primes(1,0);
    do {
      /* now loop on a_1 */
      for(a1=pdata[pc][1][1]; a1<=3; a1 += pdata[pc][1][0]) {
        init_all_primes(2,a1); do {
          for(a2=pdata[pc][2][1]; a2<plist[pc]; a2 += pdata[pc][2][0]) {
            init_all_primes(6,a2); do {
              for(a6=pdata[pc][6][1]; a6<plist[pc]; a6 += pdata[pc][6][0]) {
                init_all_primes(3,a6); do {
                  for(a3=pdata[pc][3][1]; a3<plist[pc]; a3 += pdata[pc][3][0]) {
                    init_all_primes(4,a3); do {
                      for(a4=pdata[pc][4][1]; a4<plist[pc]; a4 += pdata[pc][4][0]) {
                        init_all_primes(5,a4); do {
                          for(a5=pdata[pc][5][1]; a5<plist[pc]; a5 += pdata[pc][5][0]) {
                            printf("{%ld,%ld,%ld,%ld,%ld,%ld},",
                                   a1,a2,a3,a4,a5,a6);
                        
                          }
                        } while(next_all_primes(5,a4));
                      }
                    } while(next_all_primes(4,a3));
                  }
                } while(next_all_primes(3,a6));
              }
            } while(next_all_primes(6,a2));
          }        
        } while(next_all_primes(2,a1));
      }
    } while (next_all_primes(1,0));
  }
  /*  } */
  
  printf("\n");
  return 0;
}
