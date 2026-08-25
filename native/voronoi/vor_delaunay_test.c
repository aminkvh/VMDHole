#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "vor_delaunay.h"
static double now(void){struct timeval t;gettimeofday(&t,0);return t.tv_sec+1e-6*t.tv_usec;}
static int fails=0;
static void ck(int c,const char*w){ if(!c){printf("  FAIL  %s\n",w);fails++;} }
int main(void){
    dt_mesh m; int i;
    /* 1. a single tetrahedron's worth of points */
    { double p[]={0,0,0, 10,0,0, 0,10,0, 0,0,10};
      ck(dt_build(&m,p,4)==0,"build 4 points");
      printf("  4 pts: live=%d finite=%d verify=%d\n",dt_count_live(&m),dt_count_finite(&m),dt_verify(&m));
      ck(dt_count_finite(&m)==1,"4 general-position points give exactly 1 finite tetra");
      ck(dt_verify(&m)==0,"verify clean"); dt_free(&m); }
    /* 2. a cube: 8 points, classic case, must triangulate into 5 or 6 tetra */
    { double p[]={0,0,0, 10,0,0, 0,10,0, 10,10,0, 0,0,10, 10,0,10, 0,10,10, 10,10,10};
      ck(dt_build(&m,p,8)==0,"build cube");
      printf("  cube: live=%d finite=%d verify=%d\n",dt_count_live(&m),dt_count_finite(&m),dt_verify(&m));
      ck(dt_verify(&m)==0,"cube verify clean (cospherical corners!)"); dt_free(&m); }
    /* 3. random points: the Delaunay property must hold exactly */
    for (int n=20;n<=200;n*=10){
      double *p=malloc(3*n*sizeof(double)); srand(7);
      for(i=0;i<3*n;i++) p[i]=(rand()%200000)/1000.0-100.0;
      ck(dt_build(&m,p,n)==0,"build random");
      int v=dt_verify(&m);
      printf("  n=%d: live=%d finite=%d verify=%d\n",n,dt_count_live(&m),dt_count_finite(&m),v);
      ck(v==0,"random point set satisfies the empty-circumsphere property");
      /* Euler-ish sanity: finite tetra should be O(n) */
      ck(dt_count_finite(&m)>=n-3,"enough finite tetrahedra");
      dt_free(&m); free(p); }
    /* 4. degenerate: points on a regular grid (many cospherical/coplanar sets) */
    { int N=5,k=0; double *p=malloc(3*N*N*N*sizeof(double));
      for(int a=0;a<N;a++)for(int b=0;b<N;b++)for(int c=0;c<N;c++){p[k++]=a*10;p[k++]=b*10;p[k++]=c*10;}
      ck(dt_build(&m,p,N*N*N)==0,"build regular grid (heavily degenerate)");
      int v=dt_verify(&m);
      printf("  grid %d^3: live=%d finite=%d verify=%d\n",N,dt_count_live(&m),dt_count_finite(&m),v);
      ck(v==0,"regular grid satisfies the empty-circumsphere property");
      dt_free(&m); free(p); }
    /* 5. circumcentre of a known tetrahedron */
    { double p[]={0,0,0, 10,0,0, 0,10,0, 0,0,10}; double cc[3],cl;
      dt_build(&m,p,4);
      for(i=0;i<m.nt;i++) if(!m.t[i].dead&&dt_is_finite(&m,i)) break;
      ck(dt_circumcentre(&m,i,NULL,cc,&cl)==0,"circumcentre computed");
      printf("  circumcentre (%.3f %.3f %.3f) R=%.4f  (expect 5,5,5 R=8.6603)\n",cc[0],cc[1],cc[2],cl);
      ck(fabs(cc[0]-5)<1e-6&&fabs(cc[1]-5)<1e-6&&fabs(cc[2]-5)<1e-6,"circumcentre correct");
      ck(fabs(cl-sqrt(75.0))<1e-6,"circumradius correct"); dt_free(&m); }
    /* 6. timing on a protein-sized set */
    { int n=13658; double *p=malloc(3*n*sizeof(double)); srand(3);
      for(i=0;i<3*n;i++) p[i]=(rand()%120000)/1000.0-60.0;
      double t0=now(); int rc=dt_build(&m,p,n); double t1=now();
      ck(rc==0,"build 13658 points");
      printf("  n=13658: %.2f s, live=%d finite=%d\n",t1-t0,dt_count_live(&m),dt_count_finite(&m));
      dt_free(&m); free(p); }
    printf("\n%s\n", fails?"FAILURES":"all checks passed");
    return fails?1:0;
}
