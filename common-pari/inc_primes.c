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
  if(!next_all_primes_range(i, what, 0, totprimes))
    return 0;
  /* This is where we CRT */
  docrt(i);
  return 1;
}

/* One step of the congruence classes for a_i, which are enumerated like an
   odometer with prime 0 the fastest digit: advance the first prime j with
   jlo <= j < jhi that has a next class, and reset the primes before it (all
   of them, also those before jlo) to their first classes.  Returns 0, with
   nothing reset, if none of the primes jlo <= j < jhi has a next class.

   Unlike next_all_primes this does not combine the congruences (docrt(i)),
   so arc[i] and shift[i] are not updated.

   next_all_primes(i, what) is range(i, what, 0, totprimes) followed by
   docrt(i).  For any h (group_prime(i, n) below), the classes come in groups
   with the same classes at the primes j >= h: range(i, what, 0, h) steps
   within the group (0 at its end), and range(i, what, h, totprimes) goes to
   the first class of the next group, skipping the rest of the current one
   without stepping through it. */
long next_all_primes_range(long i, long what, long jlo, long jhi) {

  long j=jlo;

  while(j<jhi && (pflags[j][i] ||
        !((long) pf_next_list[j](j, i)))) {
    j++;
  }

  /* Now, either j = jhi, or (pflags[j][i]==0 and pf_next returned true) */
  if(j==jhi) {
    return 0;
  }
  
  /* Now we need to backtrack and init some primes */
  while(0<= --j) {
    if(!(pflags[j][i]>0))
      (int) pf_init_list[j](j, i, (long) truemod(what, plist[j]));
  }

  return 1;
  
  /*
    i is as above
    we have to step through nest of cases:
    each time we try to increase prime1.  If that fails increase prime2, and so on until
    we either run out of primes (and return 0), or we succeed, in which case we have to
    reinitialize all remaining primes.
  */
}

/* For sharing out the congruence classes for a_i in groups (see
   next_all_primes_range): the largest index h of a prime whose classes are
   stepped through (pflags 0) such that the primes j >= h give at least
   ngroups groups, estimating the number of classes at p as p (exact for the
   21 function in degree 3).  If there are fewer, the first such prime (all
   groups are then single classes), or 0 if there is none. */
long group_prime(long i, long ngroups) {
  long j, h = 0;
  double groups = 1.0;
  for(j=totprimes-1; j>=0; j--)
    if(pflags[j][i] == 0) {
      h = j;
      groups *= plist[j];
      if(groups >= ngroups) break;
    }
  return h;
}
