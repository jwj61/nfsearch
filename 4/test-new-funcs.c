#include <stdio.h>
#include <stdlib.h>

long pdata[10][8][7];
long plist[10];
long pflags[10][8];

long long arc[10], shift[10];

int t4init(long, long, long);    long t4next(long, long);
int t31init(long, long, long);   long t31next(long, long);
int t22init(long, long, long);   long t22next(long, long);
int t211init(long, long, long);  long t211next(long, long);
long next_all_primes(long , long);
long init_all_primes(long , long);

long (*(f_next_list[]))() = {t4next, t31next, t22next, t211next};
int (*(f_init_list[]))() = {t4init, t31init, t22init, t211init};

long (*(pf_next_list[8]))();
int (*(pf_init_list[8]))();
long totprimes;

long long truemod(long long a, long long m) {
  long long ans = a % m;
  if(ans<0)
    ans += m;
  return ans;
}

/* This assumes that 0 < a < b, and that inverse exists */

long long minv(long long a, long long b) {
  long long c[2], d[2], r, q, m = b;
  int j=0;

  c[0] = 1;
  d[0] = 0;
  c[1] = a/b;
  d[1] = 1;
  r = a % b;
  while(r > 0) {
    q = a/b;
    c[j] = q*c[1-j]+c[j];
    d[j] = q*d[1-j]+d[j];
    a = b;
    b = r;
    r = a % b;
    j = 1-j;
  }
  return (j ? c[1-j] : m-c[1-j]);
}


int main(int argc, char * argv []) {
  /*  int (*f_init)();
      long (*f_next)(); */
  long j,a1,a2,a3,a4,a5,a6, pc;
  long myprime;
  int myfunc;
  myfunc = atoi(argv[1]);
  myprime = atoi(argv[2]);
  
  plist[0] = myprime;
  plist[1] = 11;
  plist[2] = 5;
  plist[3] = 2003;

  totprimes = 1;
  plist[0] = myprime;

  for(j=0; j<totprimes; j++) {
    pf_next_list[j] = f_next_list[myfunc];
    pf_init_list[j] = f_init_list[myfunc];
  }

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
      for(a1=pdata[pc][1][1]; a1<3; a1 += pdata[pc][1][0]) {
        init_all_primes(2,a1); do {
          for(a2=pdata[pc][2][1]; a2<plist[pc]; a2 += pdata[pc][2][0]) {
            init_all_primes(4,a2); do {
              for(a4=pdata[pc][4][1]; a4<plist[pc]; a4 += pdata[pc][4][0]) {
                init_all_primes(3,a4); do {
                  for(a3=pdata[pc][3][1]; a3<plist[pc]; a3 += pdata[pc][3][0]) {
                            printf("{%ld,%ld,%ld,%ld},",
                                   a1,a2,a3,a4);
                        
                  }
                } while(next_all_primes(3,a4));
              }
            } while(next_all_primes(4,a2));
          }        
        } while(next_all_primes(2,a1));
      }
    } while (next_all_primes(1,0));
  }
  /*  } */
  
  printf("\n");
  return 0;
}
