/*
  Fast necessary test applied before the exact discriminant.

  For the target fields, disc(f) = d_K * (index)^2, and in the current case
  the exponent of every p in S in d_K is known (cases[i][..].discexp).
  So disc(f) lies in the square class of c = prod p^(d_p mod 2), up to sign.
  For a prime q = 1 mod 4 outside S with q not dividing disc(f), this forces
  disc(f)/c to be a square mod q (the sign does not matter since -1 is a
  square mod q).  Each such q rejects about half of the polynomials.

  The test result depends only on the residue of the innermost coefficient
  modulo q once the other coefficients are fixed, so results are cached per
  residue for the duration of one innermost loop (a lazy sieve).

  The filter is turned off in slop mode, where the acceptance test allows
  small extra primes in the field discriminant.  Compile with -DNO_QFILTER
  to turn it off entirely.
*/

#include <stdio.h>
#include <string.h>
#include "c-target.h"

#if DEG==3 || DEG==4
#define QF_INNER 3
#else
#define QF_INNER (DEG-1)
#endif

#define QF_NQ    24
#define QF_MAXQ  400

static const long qf_cand[] = {
  13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137, 149, 157,
  173, 181, 193, 197, 229, 233, 241, 257, 269, 277, 281, 293, 313, 317,
  337, 349, 353, 373, 389, 397 };

static long qf_q[QF_NQ];        /* the primes q */
static long qf_cinv[QF_NQ];     /* inverse of the square-class constant c mod q */
static int  qf_n = 0;           /* number of primes in use */
static int  qf_on = 0;
static long qf_base[QF_NQ][DEG+1];  /* fixed coefficients reduced mod q */
static unsigned long qf_gen = 0;    /* current innermost-loop generation */
static unsigned long qf_stamp[QF_NQ][QF_MAXQ];
static unsigned char qf_val[QF_NQ][QF_MAXQ];

long qf_tested = 0, qf_rejected = 0;   /* statistics */

static long
qf_pow(long b, long e, long m)
{
  long r = 1;
  b %= m; if(b < 0) b += m;
  while(e) { if(e & 1) r = r*b % m; b = b*b % m; e >>= 1; }
  return r;
}

/* Resultant mod q of A (degree da) and B (degree db), coefficients from
   the leading one down.  Uses res(A,B) = (-1)^(da db) lc(B)^(da-dr) res(B,R)
   with R = A mod B of degree dr. */
static long
qf_res(const long *f, int da, const long *g, int db, long q)
{
  long A[16], B[16], R[16], res = 1, inv, c;
  int i, j, dr;

  for(i=0; i<=da; i++) A[i] = f[i];
  for(i=0; i<=db; i++) B[i] = g[i];
  while(db >= 0 && B[0] == 0) { for(i=0; i<db; i++) B[i] = B[i+1]; db--; }
  for(;;) {
    if(db < 0) return 0;
    if(db == 0) return res * qf_pow(B[0], da, q) % q;
    inv = qf_pow(B[0], q-2, q);
    for(i=0; i<=da-db; i++) {
      c = A[i] * inv % q;
      if(c)
        for(j=0; j<=db; j++) {
          A[i+j] = (A[i+j] - c*B[j]) % q;
          if(A[i+j] < 0) A[i+j] += q;
        }
    }
    dr = db-1;
    for(i=0; i<=dr; i++) R[i] = A[da-dr+i];
    while(dr >= 0 && R[0] == 0) { for(i=0; i<dr; i++) R[i] = R[i+1]; dr--; }
    if(dr < 0) return 0;
    if((da*db) & 1) res = (q - res) % q;
    res = res * qf_pow(B[0], da-dr, q) % q;
    for(i=0; i<=db; i++) A[i] = B[i];
    da = db;
    for(i=0; i<=dr; i++) B[i] = R[i];
    db = dr;
  }
}

/* Call once per case, after initcase. */
void
qf_setup(long *casevec)
{
  int k, i;
  long q, c;

  qf_n = 0;
#ifdef NO_QFILTER
  qf_on = 0;
  return;
#endif
  qf_on = !do_slop;
  for(k=0; k < (int) (sizeof(qf_cand)/sizeof(qf_cand[0])) && qf_n < QF_NQ; k++) {
    q = qf_cand[k];
    c = 1;
    for(i=0; i<totprimes; i++) {
      if(plist[i] == q) break;
      if(cases[i][casevec[i]].discexp & 1) c = c * (plist[i] % q) % q;
    }
    if(i < totprimes) continue;   /* q is in S */
    qf_q[qf_n] = q;
    qf_cinv[qf_n] = qf_pow(c, q-2, q);
    qf_n++;
  }
}

static long
qf_mod(long long a, long q)
{
  long r = (long) (a % q);
  return (r < 0) ? r + q : r;
}

/* Call at the start of each innermost loop, when every coefficient except
   A_{QF_INNER} is fixed. */
void
qf_newloop(void)
{
  long long co[9];
  int k, i;

  if(!qf_on) return;
  co[0] = 1; co[1] = A1; co[2] = A2; co[3] = A3; co[4] = A4;
  co[5] = A5; co[6] = A6; co[7] = A7; co[8] = A8;
  for(k=0; k<qf_n; k++)
    for(i=0; i<=DEG; i++) qf_base[k][i] = qf_mod(co[i], qf_q[k]);
  qf_gen++;
}

/* Returns 0 if the current polynomial certainly fails, 1 otherwise. */
int
qf_test(void)
{
  long long inner;
  long f[DEG+1], g[DEG], q, r, v;
  int k, i;

  if(!qf_on) return 1;
  qf_tested++;
  switch(QF_INNER) {
  case 3: inner = A3; break;
  case 4: inner = A4; break;
  case 5: inner = A5; break;
  default: inner = A6; break;
  }
  for(k=0; k<qf_n; k++) {
    q = qf_q[k];
    r = qf_mod(inner, q);
    if(qf_stamp[k][r] != qf_gen) {
      for(i=0; i<=DEG; i++) f[i] = qf_base[k][i];
      f[QF_INNER] = r;
      for(i=0; i<DEG; i++) g[i] = f[i] * (DEG-i) % q;
      v = qf_res(f, DEG, g, DEG-1, q);
      /* q | disc: no information.  Otherwise disc/c must be a square. */
      qf_val[k][r] = (v == 0) ||
        (qf_pow(v * qf_cinv[k] % q, (q-1)/2, q) == 1);
      qf_stamp[k][r] = qf_gen;
    }
    if(!qf_val[k][r]) { qf_rejected++; return 0; }
  }
  return 1;
}
