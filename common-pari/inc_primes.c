#include "c-target.h"


void docrt(int ii) {
  long j;
  long long r, s;
  
  arc[ii] = pdata[0][ii][0]; /* This is the modulus */
  shift[ii] = pdata[0][ii][1];
  for(j=1; j<totprimes; j++) {
    if(pdata[j][ii][0]>1) { /* Otherwise, this modulus is 1, and we do nothing */
      /* Thinking, we have the pdata and the arc/shift.  Solve x \equiv a_i mod m_i we
         gcd trick on the moduli r m_1 + s m_2 = 1 and form r m_1 a_2 + s m_2 a_1 */
      r = minv(pdata[j][ii][0] % arc[ii], arc[ii]);
      s = (1-r*pdata[j][ii][0])/arc[ii];
      shift[ii] = truemod(shift[ii]*pdata[j][ii][0]*r
                          +pdata[j][ii][1]*arc[ii]*s, arc[ii]*pdata[j][ii][0]);
      arc[ii] = arc[ii]*pdata[j][ii][0];
      
    }
  }
}

void checkcrt(long a1,long a2,long a3,long m1,long m2,long m3) {
  pdata[0][1][0] = m1;
  pdata[0][1][1] = a1;
  pdata[1][1][0] = m2;
  pdata[1][1][1] = a2;
  pdata[2][1][0] = m3;
  pdata[2][1][1] = a3;
  totprimes=3;
  docrt(1);
  printf("\n");
  printf("Solution is ");
  printf("%lld", shift[1]);
  printf(" modulo ");
  printf("%lld", arc[1]);
  printf("\n");
}



/* i refers to a_i */
long init_all_primes(long i, long what) {
  long ans=1, j;
  for(j=0; ans && j<totprimes;j++)
    if(!(pflags[j][i]>0)) {
      ans *= (int) pf_init_list[j](j,i,(long) truemod(what, plist[j]));
    }
  if(ans) {
    /* This is where we will CRT */
    docrt(i);
  }
  
  return ans;
}

long next_all_primes(long i, long what) {

  long j=0;

  while(j<totprimes && (pflags[j][i] ||
        !((long) pf_next_list[j](j, i)))) {
    j++;
  }

  /* Now, either j = totprimes, or (pflags[j][i]==0 and pf_next returned true) */
  if(j==totprimes) {
    return 0;
  }
  
  /*  j--; */  /* Now we need to backtrack and init some primes */
  while(0<= --j) {
    if(!(pflags[j][i]>0))
      (int) pf_init_list[j](j, i, (long) truemod(what, plist[j]));
  }
  

  /* This is where we CRT */
  docrt(i);
  
  return 1;
  
  /*
    i is as above
    we have to step through nest of cases:
    each time we try to increase prime1.  If that fails increase prime2, and so on until
    we either run out of primes (and return 0), or we succeed, in which case we have to
    reinitialize all remaining primes.
  */
}
