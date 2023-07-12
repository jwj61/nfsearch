#include "c-target.h"

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
int t2211init(long pnum, long which, long what) {
  long long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 2;
    pdata[pnum][which][1] = 1;
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
    cc = what;
    break;
  case 6:/* first for 3: -4*a^3 - 2*b*a^2 + b^2*a */
    dd = what;
    bb = 0;
    aa = 0;
    pdata[pnum][6][0] = pp;
    pdata[pnum][6][1] = 0;
    pdata[pnum][5][0] = pp;
    pdata[pnum][5][1] = 0;
    pdata[pnum][4][0] = pp;
    pdata[pnum][4][1] = 0;
    pdata[pnum][3][0] = pp;
    pdata[pnum][3][1] = 0;
    break;
  default:
    pdata[pnum][which][0] = pp;
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t2211next(long pnum, long which) {
  long long a2, b2, ba;

  if(which==6) {
    ++(aa);
    if(aa==pp) {
      (bb)++;
      if(bb==pp)
        return 0; /* Done with this case! */
      aa=0;
    }
    b2=(bb*bb) % pp;
    a2=(aa*aa) % pp;
    ba = (aa*bb) % pp;
    pdata[pnum][6][1] = truemod(3*b2*a2 - 2*cc*b2*aa + -2*bb*b2 + dd*b2,pp);
    pdata[pnum][5][1] = truemod(6*ba*a2- 4*cc*bb*a2-6*b2*aa+2*dd*ba+ cc*b2,
                                pp);
    pdata[pnum][4][1] = truemod(3*a2*a2 - 2*cc*a2*aa + dd*a2
                                - 2*cc*ba + -3*b2 + 2*dd*bb, pp);
    pdata[pnum][3][1] = truemod(4*a2*aa - 3*cc*a2  -6*ba + 2*dd*aa + 2*cc*bb,
                                pp);
    return 1;
  }
  /* in all other cases, there is only one case to consider */
  return 0;
}


/* 

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
