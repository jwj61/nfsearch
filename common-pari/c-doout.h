#define DUMPIT Dumpit dUmpit

#if DEG==3
#define Dumpit fprintf(afile, \
 "x^3+"); \
 plintf(afile, A1); \
 fprintf(afile, "*x^2+"); \
 plintf(afile, A2); \
 fprintf(afile, "*x+"); \
 plintf(afile, A3);
#else
#if DEG==4
#define Dumpit fprintf(afile, \
 "x^4+"); \
 plintf(afile, A1); \
 fprintf(afile, "*x^3+"); \
 plintf(afile, A2); \
 fprintf(afile, "*x^2+"); \
 plintf(afile, A3); \
 fprintf(afile, "*x+"); \
 plintf(afile, A4);
#else
#if DEG==5
#define Dumpit fprintf(afile, \
 "x^5+"); \
 plintf(afile, A1); \
 fprintf(afile, "*x^4+"); \
 plintf(afile, A2); \
 fprintf(afile, "*x^3+"); \
 plintf(afile, A3); \
 fprintf(afile, "*x^2+"); \
 plintf(afile, A4); \
 fprintf(afile, "*x+"); \
 plintf(afile, A5);
#else
#if DEG==6
#define Dumpit fprintf(afile, \
 "x^6+"); \
 plintf(afile, A1); \
 fprintf(afile, "*x^5+"); \
 plintf(afile, A2); \
 fprintf(afile, "*x^4+"); \
 plintf(afile, A3); \
 fprintf(afile, "*x^3+"); \
 plintf(afile, A4); \
 fprintf(afile, "*x^2+"); \
 plintf(afile, A5); \
 fprintf(afile, "*x+"); \
 plintf(afile, A6);
#else
#if DEG==7
#define Dumpit fprintf(afile, \
 "x^7+"); \
 plintf(afile, A1); \
 fprintf(afile, "*x^6+"); \
 plintf(afile, A2); \
 fprintf(afile, "*x^5+"); \
 plintf(afile, A3); \
 fprintf(afile, "*x^4+"); \
 plintf(afile, A4); \
 fprintf(afile, "*x^3+"); \
 plintf(afile, A5); \
 fprintf(afile, "*x^2+"); \
 plintf(afile, A6); \
 fprintf(afile, "*x+"); \
 plintf(afile, A7);
#else
#if DEG==8
#define Dumpit fprintf(afile, \
 "x^8+"); \
 plintf(afile, A1); \
 fprintf(afile, "*x^7+"); \
 plintf(afile, A2); \
 fprintf(afile, "*x^6+"); \
 plintf(afile, A3); \
 fprintf(afile, "*x^5+"); \
 plintf(afile, A4); \
 fprintf(afile, "*x^4+"); \
 plintf(afile, A5); \
 fprintf(afile, "*x^3+"); \
 plintf(afile, A6); \
 fprintf(afile, "*x^2+"); \
 plintf(afile, A7); \
 fprintf(afile, "*x+"); \
 plintf(afile, A8);
#endif
#endif
#endif
#endif
#endif
#endif

#define dUmpit fprintf(afile, "\n");

#define dUmpitold \
 fprintf(afile, ",");fprintf(afile,"[%1ld",plist[0]); \
 for(ii=1; ii<totprimes; ii++) fprintf(afile, ",%1ld",plist[ii]); \
 fprintf(afile, "],[%1ld", j_p[0]); \
 for(ii=1; ii<totprimes; ii++) fprintf(afile, ",%1ld",j_p[ii]); \
 fprintf(afile, "]], \\\n"); \
 fflush(afile);
