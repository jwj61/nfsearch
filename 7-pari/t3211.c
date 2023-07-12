#include "c-target.h"

#define pp plist[pnum]
#define bb pdata[pnum][8][0]
#define aa pdata[pnum][8][1]
#define cc pdata[pnum][8][2]
#define dd pdata[pnum][8][3]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* Let's START *************************************************************/
int t3211init(long pnum, long which, long what) {
  long jj;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  
  switch(which) {
  case 0:
    pdata[pnum][0][0] = 3;
    pdata[pnum][0][1] = 4;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1; /* modulus */
    pdata[pnum][2][1] = 0; /* starting value */

    pflags[pnum][1] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    for(jj=3;jj<=DEG; jj++)
      pflags[pnum][jj] = CALL_NONE;
    pflags[pnum][DEG] = CALL_BOTH;
    for(jj=3;jj<=DEG; jj++)
      pdata[pnum][jj][0] = pp;
    break;
    
  case 2:  /* c */
    cc = what; /* for future reference */
    break;
  case DEG:
    dd = what;
    
    aa = 0; /* ready for double loop */
    bb = 0;
    for(jj=3; jj<=DEG; jj++) {
      pdata[pnum][jj][1] = 0;
    }
    
    break;
  default:
    /*    pdata[pnum][which][0] = pp; */
    fprintf(stderr,"Error in t3211 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t3211next(long pnum, long which) {
  long long a2, a3, b2, b3, a4, a5, bc, b2c, bd, b4, b3c, b2d;
  
  /* Should only be called with which ==2 */
  if(which != DEG) {
    fprintf(stderr,"Error in t3211 next.\n");
    return 0;
  }
  
  (aa)++;
  if(aa == pp) {
    (bb)++;
    if(bb == pp) {
      return 0; /* We are done with the double loop */
    } /* Really an else now */
    aa = 0;
  }
  /* We now have a new aa and a new bb, time to set everyone: */
  a2 = (aa*aa) % pp;  a3 = (a2 * aa) % pp; a4 = (a3*aa) % pp;
  a5 = (a4*aa) % pp;
  b2 = (bb*bb) % pp;  b3 = (b2 * bb) % pp;
  b4 = (bb*b3) % pp;
  bc = (bb*cc) % pp; 
  b2c = (b2*cc) % pp; b3c = (b3*cc) % pp;
  bd = (bb*dd) % pp; b2d = (b2*dd) % pp;
  
  pdata[pnum][3][1] = truemod(10*a3 + (12*bb - 6*cc)*a2 +
                              (9*b2 - 6*bc + 3*dd)*aa +
                              (4*b3 -3*b2c + 2*dd*bb),pp);
  pdata[pnum][4][1] = truemod(15*a4 + (36*bb - 8*cc)*a3 +
                              (33*b2 - 18*bc +3*dd)*a2 +
                              (18*b3 -  12*b2c + 6*bd)*aa +
                              (3*b4 - 2*cc*b3 + dd*b2),pp);
  pdata[pnum][5][1] = truemod(6*a5 + (36*bb - 3*cc)*a4 +
                              (45*b2 - 18*bc + dd)*a3 +
                              (30*b3 -18*b2c + 6*bd)*a2 +
                              (9*b4 - 6*b3c + 3*b2d)*aa,pp);
  pdata[pnum][6][1] = truemod(12*bb*a5 +
                              (27*b2 - 6*bc)*a4 +
                              (22*b3 - 12*b2c + 2*bd)*a3 +
                              (9*b4 -6*b3c + 3*b2d)*a2,pp);
  pdata[pnum][7][1] = truemod(6*b2*a5 +
                              (6*b3 - 3*b2c)*a4 + (3*b4 -2*b3c + b2d)*a3,pp);
  
  return 1;
}

