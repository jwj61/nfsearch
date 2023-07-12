#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
/*
  Inputs a search program file of congruences, and sorts them.
*/

#define cmp(a,b)   ((a<b)?-1:((a>b)? 1:0))

#define NUM  500000  /* Maximum number of congruences */
#define COMMENT_CHAR '%'

void showai(int);

int * (ai[NUM]);
long degree;

void getmemory() {
  int i;

  for(i=0; i<NUM; i++) {
    if((ai[i] = (int *) malloc(degree*sizeof(int))) == NULL) {
      printf("Cannot allocate memory for ai\n");
      exit(1);
    }
  }
}

int mycompare3 (const void *xx, const void *yy) {
  int *x = *((int *) xx);
  int *y = *((int *) yy);
  printf("Cmp (%d,%d,%d) ? (%d,%d,%d) and %d, %d, %d\n",
         x[0],x[1],x[2],y[0],y[1],y[2],
         cmp(x[0],y[0]),cmp(x[1],y[1]),cmp(x[2],y[2]));
  if(cmp(x[0],y[0]))
    return cmp(x[0],y[0]);
  if(cmp(x[1],y[1]))
    return cmp(x[1],y[1]);
  return cmp(x[2],y[2]);
  
}

int mycompare6 (const void *xx, const void *yy) {
  int *x = *((int *) xx);
  int *y = *((int *) yy);
  
  if(cmp(x[0],y[0]))
    return cmp(x[0],y[0]);
  if(cmp(x[1],y[1]))
    return cmp(x[1],y[1]);
  if(cmp(x[5],y[5]))
    return cmp(x[5],y[5]);
  if(cmp(x[2],y[2]))
    return cmp(x[2],y[2]);
  if(cmp(x[3],y[3]))
    return cmp(x[3],y[3]);
  return cmp(x[4],y[4]);
}

int mycompare7 (const void *xx, const void *yy) {
  int *x = *((int *) xx);
  int *y = *((int *) yy);
  
  if(cmp(x[0],y[0]))
    return cmp(x[0],y[0]);
  if(cmp(x[1],y[1]))
    return cmp(x[1],y[1]);
  if(cmp(x[6],y[6]))
    return cmp(x[6],y[6]);
  if(cmp(x[2],y[2]))
    return cmp(x[2],y[2]);
  if(cmp(x[3],y[3]))
    return cmp(x[3],y[3]);
  if(cmp(x[4],y[4]))
    return cmp(x[4],y[4]);
  return cmp(x[5],y[5]);
}

void showai(int deg){
  int j, k;
  
  for(j=0; j<5; j++) {
    for(k=0; k<deg; k++)
      printf("%d ", ai[j][k]);
    printf("\n");
  }
}

void mysort(int num) {
  qsort(ai,num,sizeof(int*), mycompare7);
}

void dump_cases(int num, int deg, int c, int d, FILE *outfile) {
  int j, k, cnt;

  fprintf(outfile, "\n\nc = %d\nd = %d\n\n", c, d);  
  cnt=0;
  
  for(j=0; j<num; j++) {
    if(cnt==4) {
      cnt=0;
      fprintf(outfile, "\n");
    }
    cnt++;
    fprintf(outfile, "{");
    for(k=0; k<deg-1; k++) {
      fprintf(outfile, "%d,", ai[j][k]);  
    }
    fprintf(outfile, "%d}, ", ai[j][deg-1]);
  }
}

long fmyreadint(FILE *inf)
{
  long num=0, tmpc;

  while(!isdigit(tmpc=getc(inf))) ; /* Find a digit */
  num = tmpc - '0';
  while(isdigit(tmpc=getc(inf)))
    num = 10*num + tmpc - '0';
  ungetc(tmpc, inf);
  return(num);
}

int main(int argc, char * argv []) {

  char * infilename;
  char outfilename[128];
  FILE *infile, *outfile;
  long  c, d, old_c, old_d;
  int k, nc, num_cases;

  
  infilename = argv[1];
  degree = atoi(argv[2]);
  if((infile = fopen(infilename, "r")) == NULL) {
    printf("Cannot open input file %s.\n", infilename);
    exit(1);
  }
  
  outfilename[0] = '\0';
  sprintf(outfilename, "%s_sorted", infilename);
  if((outfile = fopen(outfilename, "w")) == NULL) {
    printf("Cannot open output file %s\n", outfilename);
    exit(1);
  }
  
  getmemory();

/***********************************************************************/
/****************      Main Part     ***********************************/
/***********************************************************************/

  
  c = 1;       old_c = 1;
  d = 0; old_d = 0;
  num_cases=0;
  nc=getc(infile);
  while(nc != EOF) {
    switch (nc) {
    case COMMENT_CHAR: /* Skip comment lines. */
        fprintf(outfile, "%c", nc);
        while(((nc=getc(infile)) != EOF) && (nc != '\n')) 
          fprintf(outfile, "%c", nc);
        fprintf(outfile, "\n");
        break;
    case EOF: break;
    case 'c':
      old_c = c;
      c=fmyreadint(infile);
      if(num_cases && (old_c != c)) { /* We have some cases */
        mysort(num_cases);
        dump_cases(num_cases, degree, old_c, d, outfile);
        num_cases = 0;
      }
      break;
    case 'd':
      old_d = d;
      d=fmyreadint(infile);
      if(num_cases && (d != old_d)) {
        mysort(num_cases);
        dump_cases(num_cases, degree, c, old_d, outfile);
        num_cases = 0;
      }
      break;
    case '\n': break;
    default:
      break;
    }
    if(isdigit(nc)) { /* line starting with a digit, load it */
      ungetc(nc, infile);
      for(k=0; k<degree; k++)
        ai[num_cases][k] = fmyreadint(infile);
      num_cases++;
    }
    if(nc != EOF) nc = getc(infile);
  }
  /* End up input file, clean up now */
  if(num_cases) {
    mysort(num_cases);
    dump_cases(num_cases, degree, c, d, outfile);
  }
  fprintf(outfile, "\n");
  
  

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

  fclose(infile);
  fclose(outfile);
  return 0;
}

