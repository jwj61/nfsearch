#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][7][0]
#define i6 pdata[pnum][7][1]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t6init(long pnum, long which, long what) {
  long j, i2, i3, i66, ww;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 5;
    pdata[pnum][which][1] = 1;
    i2 = (pp+1)/2;
    i3 = ((pp%3) == 1) ?  (2*pp+1)/3 : (pp+1)/3;
    i66 = (i2*i3) % pp;
    i6 = i66;
    pflags[pnum][1] = -1;
    for(j=2; j<=6; j++)
      pflags[pnum][j] = 0;
    break;
  case 1:
    pdata[pnum][which][0] = 1; /* modulus */
    pdata[pnum][which][1] = 0; /* starting value */
    break;
  case 2:
    pdata[pnum][which][0] = pp;
    aa = (what*i6) % pp;  /* Store this for future reference */
    pdata[pnum][which][1] = (15*aa*aa) % pp;
    ww = 3;
    pdata[pnum][ww][0] = pp;
    j = (aa*aa) % pp;
    pdata[pnum][ww][1] = (20*aa) % pp; /* 20 a^3 */
    pdata[pnum][ww][1] = (pdata[pnum][ww][1]*j) % pp;
    ww = 4;
    pdata[pnum][ww][0] = pp;
    j = (aa*aa) % pp;
    j = (j*j) % pp;
    pdata[pnum][ww][1] = (15*j) % pp;
    ww = 5;
    pdata[pnum][ww][0] = pp;
    i2 = (aa*aa) % pp;
    i2 = (6*i2*i2) % pp;
    pdata[pnum][ww][1] = (i2*aa) % pp;
    ww = 6;
    pdata[pnum][ww][0] = pp;
    j = (aa*aa) % pp;
    i2 = (j*j) % pp;
    pdata[pnum][ww][1] = (j*i2) % pp;
    break;
  default:
    pdata[pnum][which][0] = pp;
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t6next(long pnum, long which) {
  return 0;
}

