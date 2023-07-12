#include "c-target.h"
#include <stdio.h>

#define pp plist[pnum]
#define bb pdata[pnum][7][0]
#define aa pdata[pnum][7][1]
#define cc pdata[pnum][7][2]
#define dd pdata[pnum][7][3]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* Let's START *************************************************************/
int t3111init(long pnum, long which, long what) {
  long long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][0][0] = 2;
    pdata[pnum][0][1] = 1;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    pdata[pnum][6][0] = 1;
    pdata[pnum][6][1] = 0;
    pdata[pnum][3][0] = pp;
    pdata[pnum][4][0] = pp;
    pdata[pnum][5][0] = pp;

    for(j=1; j<=5; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][6] = CALL_INIT; /* This start like this, but can change */
    pflags[pnum][3] = CALL_BOTH;
    break;
    /*
      [b, c, d,
  -15*a^4 + 10*b*a^3 - 6*c*a^2 + 3*d*a,
  -24*a^5 + 15*b*a^4 - 8*c*a^3 + 3*d*a^2,
  -10*a^6 + 6*b*a^5 - 3*c*a^4 + d*a^3]

if a == 0, loop on d (rest are 0)
if a != 0, any constant term is ok, but this determines d
    */
  case 2:  /* c */
    bb = what; /* for future reference */
    break;
  case 6:
    cc = what;
    break;
  case 3:
    dd = what; /* for future reference */
    if(dd==0) {
      aa=0;
      pdata[pnum][3][0] = 1; /* Sets up the a=0 part */
      pdata[pnum][3][1] = 0;
      pdata[pnum][4][1] = 0;
      pdata[pnum][5][1] = 0;
    } else {
      aa=1;
      pdata[pnum][3][0] = pp; /* Sets up the a=1 part */
      pdata[pnum][3][1] = truemod(10-6*bb+3*cc+dd, pp);
      pdata[pnum][4][1] = truemod(3*cc -8*bb + 3*dd + 15, pp);
      pdata[pnum][5][1] = truemod(cc -3*bb + 3*dd + 6, pp);
    }
    break;
  default:
    /*    pdata[pnum][which][0] = pp; */
    fprintf(stderr,"Error in t3111 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t3111next(long pnum, long which) {
  long long a2, a3, ainv, ainv2, ainv3, t1;

  if(which != 3) {
    fprintf(stderr,"Error in t3111 next.\n");
    return 0;
  }

  if(aa==0) {
    pdata[pnum][3][0] = pp;
  }

  aa++;
  if(aa==pp)
    return 0;

  a2 = (aa*aa) % pp;
  a3 = (a2*aa) % pp;
  t1 = (3*aa*cc + 10*a3 - 6*bb*a2) % pp;
  ainv = minv(aa, pp);
  ainv2 = (ainv*ainv) % pp;
  ainv3 = (ainv2*ainv) % pp;
    
  pdata[pnum][3][1] = truemod(t1+dd*ainv3, pp);
  t1 = (3*a2*cc + 15*a2*a2 - 8*bb*a3) % pp;
  pdata[pnum][4][1] = truemod(t1+3*dd*ainv2, pp);
  t1 = (aa*cc + 6*a3 - 3*bb*a2) % pp;
  pdata[pnum][5][1] = truemod(a2*t1+3*dd*ainv, pp);
  
  return 1;
}


/* 


   
3111

 [b, c, (10*a^6 - 6*b*a^5 + 3*c*a^4 + d)/a^3,
   (15*a^6 - 8*b*a^5 + 3*c*a^4 + 3*d)/a^2,
   (6*a^6 - 3*b*a^5 + c*a^4 + 3*d)/a, d]

 and

 [b, c, d, 0, 0, 0]
 
  --------------------------------------
222

 [a, b, -1/8*a^3 + 1/2*b*a + 2*c, 1/64*a^4 - 1/8*b*a^2 + c*a + 1/4*b^2,
   -1/4*c*a^2 + c*b, c^2]

  --------------------------------------
2211

  [c, d, 4*a^3 - 3*c*a^2 + (-6*b + 2*d)*a + 2*c*b,
    3*a^4 - 2*c*a^3 + d*a^2 - 2*c*b*a + (-3*b^2 + 2*d*b),
    6*b*a^3 - 4*c*b*a^2 + (-6*b^2 + 2*d*b)*a + c*b^2,
    3*b^2*a^2 - 2*c*b^2*a + (-2*b^3 + d*b^2)]

  --------------------------------------
21111

 

[b, c, d, e,
  6*a^5 - 5*b*a^4 + 4*c*a^3 - 3*d*a^2 + 2*e*a,
  5*a^6 - 4*b*a^5 + 3*c*a^4 - 2*d*a^3 + e*a^2]


From a,b,c,d,e:
 [b, e, c, d, 0, 0]

 or
 
 [b, e, c,
   (-5*a^6 + 4*b*a^5 - 3*e*a^4 + 2*c*a^3 + d)/a^2,
   (-4*a^6 + 3*b*a^5 - 2*e*a^4 + c*a^3 + 2*d)/a,
   d]
*/
