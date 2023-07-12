#include "c-target.h"

#define pp plist[pnum]
#define bb pdata[pnum][7][0]
#define aa pdata[pnum][7][1]
#define cc pdata[pnum][7][2]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* Let's START *************************************************************/
int t321init(long pnum, long which, long what) {
  long long j;
  long jj;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  
  switch(which) {
  case 0:
    pdata[pnum][0][0] = 3;
    pdata[pnum][0][1] = 1;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */

    pflags[pnum][1] = 1;
    pflags[pnum][2] = 0;
    pflags[pnum][6] = 1;
    pflags[pnum][3] = 1;
    for(j=4; j<=5; j++)
      pflags[pnum][j] = 1;
    break;
    
  case 2:  /* c */
    /* a_2 == -6*a^2 + (-6*b + 3*c)*a + (-3*b^2 + 2*c*b)  */
    cc = what; /* for future reference */
    aa = 0; /* ready for double loop */
    bb = 0;
    for(jj=2; jj<=6; jj++) {
      pdata[pnum][jj][0] = pp;
      pdata[pnum][jj][1] = 0;
    }
    
    break;
  default:
    /*    pdata[pnum][which][0] = pp; */
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t321next(long pnum, long which) {
  long long a2, a3, b2, b3, a4, bc, b2c;
  
  /* Should only be called with which ==2 */
  if(which != 2)
    return 0;
  (aa)++;
  if(aa == pp) {
    (bb)++;
    if(bb == pp) {
      return 0; /* We are done with the double loop */
    } /* Really an else now */
    aa = 0;
  }
  /* We now have a new aa and a new bb, time to set everyone:
     [c, -6*a^2 + (-6*b + 3*c)*a + (-3*b^2 + 2*c*b),
  -8*a^3 + (-18*b + 3*c)*a^2 + (-12*b^2 + 6*c*b)*a + (-2*b^3 + c*b^2),
  -3*a^4 + (-18*b + c)*a^3 + (-18*b^2 + 6*c*b)*a^2 + (-6*b^3 + 3*c*b^2)*a,
  -6*b*a^4 + (-12*b^2 + 2*c*b)*a^3 + (-6*b^3 + 3*c*b^2)*a^2,
  -3*b^2*a^4 + (-2*b^3 + c*b^2)*a^3] */
  a2 = (aa*aa) % pp;  a3 = (a2 * aa) % pp; a4 = (a3*aa) % pp;
  b2 = (bb*bb) % pp;  b3 = (b2 * bb) % pp;
  bc = (bb*cc) % pp;
  b2c = (b2*cc) % pp;
  
  pdata[pnum][2][1] = truemod(-6*a2 + (3*cc-6*bb)*aa -3*b2 + 2*bc,pp);
  pdata[pnum][3][1] = truemod(-8*a3 + (-18*bb + 3*cc)*a2 +
                              (-12*b2 + 6*bc)*aa -2*b3 + b2c,pp);
  pdata[pnum][4][1] = truemod(-3*a4 + (-18*bb + cc)*a3 +
                              (-18*b2 + 6*bc)*a2 + (-6*b3 + 3*b2c)*aa,pp);
  pdata[pnum][5][1] = truemod(-6*bb*a4 + (-12*b2 + 2*bc)*a3
                              + (-6*b3 + 3*b2c)*a2,pp);
  pdata[pnum][6][1] = truemod(-3*b2*a4 + (-2*b3 + b2c)*a3,pp);
  
  return 1;
}


/* 

321
 
[c, -6*a^2 + (-6*b + 3*c)*a + (-3*b^2 + 2*c*b),
  -8*a^3 + (-18*b + 3*c)*a^2 + (-12*b^2 + 6*c*b)*a + (-2*b^3 + c*b^2),
  -3*a^4 + (-18*b + c)*a^3 + (-18*b^2 + 6*c*b)*a^2 + (-6*b^3 + 3*c*b^2)*a,
  -6*b*a^4 + (-12*b^2 + 2*c*b)*a^3 + (-6*b^3 + 3*c*b^2)*a^2,
  -3*b^2*a^4 + (-2*b^3 + c*b^2)*a^3]

  --------------------------------------
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
