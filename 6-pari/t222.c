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
int t222init(long pnum, long which, long what) {
  long long j, a2, i2, a22, b2;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 3;
    pdata[pnum][which][1] = 1;
    i2 = (pp+1)/2;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    
    pflags[pnum][1] = 1;
    pflags[pnum][2] = -1;
    pflags[pnum][6] = 0;
    for(j=3; j<=5; j++)
      pflags[pnum][j] = 1;
    break;
  case 2:  
    aa = what; /* Store this for future reference, it is b */
    break;
  case 6:/* first for 3: -4*a^3 - 2*b*a^2 + b^2*a */
    cc = 0;
    bb = what;
    pdata[pnum][6][0] = pp;
    pdata[pnum][6][1] = 0;
    /* a_3: -1/8*a^3 + 1/2*b*a */
    i2 = (pp+1)/2;
    a2 = (aa*i2) % pp;
    b2 = (bb*i2) % pp;
    a22 = (a2*a2) % pp;
    pdata[pnum][3][0] = pp;
    pdata[pnum][3][1] = truemod(bb*a2-a22*a2, pp);
    /* a_4: 1/64*a^4 - 1/8*b*a^2 + 1/4*b^2 */
    pdata[pnum][4][0] = pp;
    j = (a22*i2) % pp;
    pdata[pnum][4][1] =truemod(j*j - a22*b2 + b2*b2, pp);
    
    /* a_5: 0 */
    pdata[pnum][5][0] = pp;
    pdata[pnum][5][1] = 0;
    break;
  default:
    pdata[pnum][which][0] = pp;
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t222next(long pnum, long which) {
  long long j, i2, a2, b2, a22;

  if(which==6) {
    ++(cc);
    if(cc==pp)
      return 0; /* Done with this case! */
    pdata[pnum][6][1] = (cc*cc) % pp;
    i2 = (pp+1)/2;
    a2 = (aa*i2) % pp;
    b2 = (bb*i2) % pp;
    a22 = (a2*a2) % pp;
    pdata[pnum][5][1] = truemod(cc*bb-cc*a22, pp);
    j = (a22*i2) % pp;
    pdata[pnum][4][1] = truemod(j*j-b2*a22+cc*aa+b2*b2, pp);
    pdata[pnum][3][1] = truemod(bb*a2-a22*a2+2*cc, pp);
    return 1;
  }
  /* in all other cases, there is only one case to consider */
  return 0;
}


/* 

222

 [a, b,
 -1/8*a^3 + 1/2*b*a + 2*c,
 1/64*a^4 - 1/8*b*a^2 + c*a + 1/4*b^2,
   -1/4*c*a^2 + c*b, c^2]

  --------------------------------------
2211

  [c, d, 4*a^3 - 3*c*a^2 + (-6*b + 2*d)*a + 2*c*b,
    3*a^4 - 2*c*a^3 + d*a^2 - 2*c*b*a + (-3*b^2 + 2*d*b),
    6*b*a^3 - 4*c*b*a^2 + (-6*b^2 + 2*d*b)*a + c*b^2,
    3*b^2*a^2 - 2*c*b^2*a + (-2*b^3 + d*b^2)]

  --------------------------------------
21111

 From a,b,c,d,e

 [b, e, c, d, 0, 0]

 or
 
 [b, e, c,
   (-5*a^6 + 4*b*a^5 - 3*e*a^4 + 2*c*a^3 + d)/a^2,
   (-4*a^6 + 3*b*a^5 - 2*e*a^4 + c*a^3 + 2*d)/a,
   d]
*/
