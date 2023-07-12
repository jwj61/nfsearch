#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][8][0]
#define bb pdata[pnum][8][1]
#define cc pdata[pnum][8][2]
#define i3 pdata[pnum][8][3]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t331init(long pnum, long which, long what) {
  long j;
  long long b2, b3;

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
    i3 = minv(3, pp);
    break;
  case 2:
    cc = what;
    break;
  case DEG: /* Capture b, loop through values of a */
    bb = (what*i3) % pp;  /* Store this for future reference */
    aa = 0;
    b2 = (bb*bb) % pp;
    b3 = (b2*bb) % pp;
    pdata[pnum][3][1] = truemod(3*bb*cc,pp);
    pdata[pnum][4][1] = truemod(3*b2,pp);
    pdata[pnum][5][1] = truemod(3*b2*cc,pp);
    pdata[pnum][6][1] = truemod(b3,pp);
    pdata[pnum][7][1] = truemod(b3*cc,pp); 
   break;
  default:
    fprintf(stderr,"Error in t331 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t331next(long pnum, long which) { /* This should never be called */
  long long a2, ab, a4, a5, a3, b2, b3, a6, a7, c2, c3, cb, tmp1;
  long long cb2, c3b, c2b2, c4;
  
  if(which != DEG) {
    fprintf(stderr,"Error in t331 next.\n");
    return 0;
  }
  aa++;
  if(aa==pp)
    return 0;
  a2 = (aa*aa) % pp; /* a^2 */
  a3 = (aa*a2) % pp;
  b2 = (bb*bb) % pp;
  b3 = (bb*b2) % pp;
  c2 = (cc*cc) % pp;
  c3 = (cc*c2) % pp;
  c4 = (c2*c2) % pp;
  ab = (bb*aa) % pp; /* ab */
  a4 = (a2*a2) % pp; /* a^4 */
  a5 = (a4*aa) % pp; /* a^5 */
  a6 = (a5*aa) % pp; /* a^5 */
  a7 = (a6*aa) % pp; /* a^5 */
  cb = (cc*bb) % pp;
  
  pdata[pnum][3][1] = truemod(-14*a3 + 12*cc*a2 + (-3*bb - 3*c2)*aa +
                              3*cc*bb, pp);
  pdata[pnum][4][1] = truemod(-21*a4 +16*cc*a3 + (-3*bb - 3*c2)*a2
                              + 3*b2, pp);
  tmp1 = (c2*bb) % pp;
  pdata[pnum][5][1] = truemod(-42*a5 + 51*cc*a4 + (-33*bb -21*c2)*a3
                              + (27*cb + 3*c3)*a2 +(-6*b2 - 6*tmp1)*aa
                              + 3*cc*b2, pp);
  pdata[pnum][6][1] = truemod(-28*a6 + 36*cc*a5 + (-24*bb - 15*c2)*a4
                              + (18*cb + 2*c3)*a3 +
                              (-3*b2 - 3*tmp1)*a2 + b3, pp);
  cb2 = (cc*b2) % pp;
  c3b = (c3*bb) % pp;
  c2b2 = (c2*b2) % pp;
  pdata[pnum][7][1] = truemod(-24*a7 + 44*cc*a6 + (-36*bb -30*c2)*a5
                              + (48*cb + 9*c3)*a4 + (-18*b2 - 21*tmp1 - c4)*a3
                              +(15*cb2 + 3*c3b)*a2
                              + (-3*b3 - 3*c2b2)*aa + cc*b3, pp);
  return 1;
}
