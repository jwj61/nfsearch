#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][8][0]
#define bb pdata[pnum][8][1]
#define cc pdata[pnum][8][2]
#define i2 pdata[pnum][8][3]
#define i4 pdata[pnum][8][4]
#define i8 pdata[pnum][8][5]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t322init(long pnum, long which, long what) {
  long j;
  long long b2, b3, b4;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 4;
    pdata[pnum][which][1] = 3;
    for(j=1; j<=DEG; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][DEG] = CALL_BOTH;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1; /* modulus */
    pdata[pnum][2][1] = 0; /* starting value */
    for(j=3; j<=DEG; j++)
      pdata[pnum][j][0] = pp;
    i2 = (pp+1)/2;
    i4 = (i2*i2) % pp;
    i8 = (i2*i4) % pp;
    break;
  case 2:
    bb = (i2*what) % pp;
    break;
  case DEG: /* Capture b, loop through values of a */
    cc = (what*i2) % pp;  /* Store this for future reference */
    aa = 0;
    /* Now [2*b, 2*c, -b^3 + 2*c*b, 1/4*b^4 - c*b^2 + c^2, 0, 0, 0] */
    b2 = (bb*bb) % pp;
    b3 = (b2*bb) % pp;
    b4 = (b2*b2) % pp;
    pdata[pnum][3][1] = truemod(-b3 + 2*cc*bb,pp);
    pdata[pnum][4][1] = truemod(i4*b4 - cc*b2 + cc*cc,pp);
    pdata[pnum][5][1] = 0;
    pdata[pnum][6][1] = 0;
    pdata[pnum][7][1] = 0; 
   break;
  default:
    fprintf(stderr,"Error in t322 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t322next(long pnum, long which) { /* This should never be called */
  long long a2, a4, a5, a3, b2, b3, a6, c2, c3, tmp1;
  long long tmp2, tmp3, tmp4, tmp5, b4;
  
  if(which != DEG) {
    fprintf(stderr,"Error in t322 next.\n");
    return 0;
  }
  aa++;
  if(aa==pp)
    return 0;
  a2 = (aa*aa) % pp; /* a^2 */
  a3 = (aa*a2) % pp;
  b2 = (bb*bb) % pp;
  b3 = (bb*b2) % pp;
  b4 = (bb*b3) % pp;
  c2 = (cc*cc) % pp;
  c3 = (cc*c2) % pp;
  a4 = (a2*a2) % pp; /* a^4 */
  a5 = (a4*aa) % pp; /* a^5 */
  a6 = (a5*aa) % pp; /* a^5 */

  tmp1 = (bb*a2) % pp;
  tmp2 = (-3*i2*b2 + 3*cc) % pp;
  pdata[pnum][3][1] = truemod(35*i8*a3 - 15*i4*tmp1 + tmp2*aa
                              -b3 +2*cc*bb, pp);
  tmp1 = (-33*i8*b2 + 3*i4*cc) % pp;
  tmp2 = (-3*i2*b3+ 3*cc*bb) % pp;
  tmp3 = (bb*a3) % pp;
  tmp4 = (a2*i8) % pp;
  tmp4 = (tmp4*tmp4) % pp;
  pdata[pnum][4][1] = truemod(105*tmp4 + 25*i8*tmp3 +
                              tmp1*a2 + tmp2*aa +
                              i4*b4 - cc*b2 + c2, pp);
  tmp1 = (-27*i8*b2 + 17*i4*cc) % pp;
  tmp2 = (3*i2*b3 - 3*cc*bb) % pp;
  tmp3 = (3*i4*b4 - 3*cc*b2 + 3*c2) % pp;
  tmp4 = (i8*bb) % pp;
  tmp5 = (a5*i8) % pp;
  tmp5 = (tmp5*i8) % pp;
  pdata[pnum][5][1] = truemod(-21*tmp5 + 15*tmp4*a4 +
                              tmp1*a3 + tmp2*a2 + tmp3*aa, pp);
  tmp1 = (-3*i8*b2 + 33*i4*cc) % pp;
  tmp2 = (7*i2*b3 - 7*cc*bb) % pp;
  tmp3 = (3*i4*b4 - 3*cc*b2 + 3*c2) % pp;
  tmp4 = (i8*bb) % pp;
  tmp5 = (a3*i8) % pp;
  tmp5 = (tmp5*tmp5) % pp;
  
  pdata[pnum][6][1] = truemod(315*tmp5 - 69*tmp4*a5 + tmp1*a4 +
                              tmp2*a3 + tmp3*a2, pp);
  tmp1 = (3*i8*b2 + 15*i4*cc) % pp;
  tmp2 = (3*i2*b3 - 3*cc*bb) % pp;
  tmp3 = (i4*b4 - cc*b2 + c2) % pp;
  tmp4 = (i8*a6) % pp;
  tmp5 = (aa*i8) % pp;
  pdata[pnum][7][1] = truemod(225*tmp5*tmp4 - 45*bb*tmp4 +
                              tmp1*a5 + tmp2*a4 + tmp3*a3, pp);
  return 1;
}
