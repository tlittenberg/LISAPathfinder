// LPF.c Routines for defining basic LISA Pathfinder structural properties
// Written by John Baker 2016
#include "LPF.h"
#include "Spacecraft.h"
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//***********************************************
//   Routines setting Spacecraft data
//***********************************************


void initialize_spacecraft(struct Spacecraft *spacecraft)
{
  int i,j;

  /* Moment of inertia tensor */
  spacecraft->I    = malloc(2*sizeof(double **));
  spacecraft->invI = malloc(2*sizeof(double **));

  for(j=0; j<2; j++) spacecraft->I[j]    = malloc(3*sizeof(double *));
  for(j=0; j<2; j++) spacecraft->invI[j] = malloc(3*sizeof(double *));
  for(j=0; j<2; j++) for(i=0; i<3; i++) spacecraft->I[j][i]    = malloc(3*sizeof(double));
  for(j=0; j<2; j++) for(i=0; i<3; i++) spacecraft->invI[j][i] = malloc(3*sizeof(double));

  /* Position of proof masses */
  spacecraft->R = malloc(2*sizeof(double *));
  for(i=0; i<2; i++) spacecraft->R[i] = malloc(3*sizeof(double));

  /* Spacecraft M-frame corners */  //disused?
  spacecraft->x = malloc(9*sizeof(double *));
  for(j=0; j<9; j++) spacecraft->x[j] = malloc(2*sizeof(double));

  /* Info on the spacecraft faces and their connections */
  int nfaces=10;
  spacecraft->nfaces=nfaces;
  spacecraft->faces = malloc(nfaces*sizeof(struct FaceData*));
  for(i=0;i<nfaces;i++)spacecraft->faces[i]=malloc(sizeof(struct FaceData));
  
  double xmin=SC_BOT_CORNER_1_X,xmax=SC_BOT_CORNER_5_X,ymin=SC_BOT_CORNER_7_Y,ymax=SC_BOT_CORNER_3_Y;
  //face 0 extends on side from corner 1 to corner 2.
  initialize_face(spacecraft->faces[0],0,SC_BOT_CORNER_1_X,SC_BOT_CORNER_1_Y,SC_BOT_CORNER_2_X,SC_BOT_CORNER_2_Y,
		  7, 0, 1, 0,     1, 0, 0, 0,    8, 0, 0, 0,     9, 0, 0, 0 );   
  printf("r0=(%g,%g,%g)\n",spacecraft->faces[0]->r0[0],spacecraft->faces[0]->r0[1],spacecraft->faces[0]->r0[2]);
  printf("rmax=(%g,%g)\n",spacecraft->faces[0]->rmax[0],spacecraft->faces[0]->rmax[1]);
  //face 1
  initialize_face(spacecraft->faces[1],0,SC_BOT_CORNER_2_X,SC_BOT_CORNER_2_Y,SC_BOT_CORNER_3_X,SC_BOT_CORNER_3_Y,
                  0, 0, 1, 0,     2, 0, 0, 0,   -1, 0, 0, 0,    -1, 0, 0, 0 );//no wrap over diagonal edge.   
  printf("r0=(%g,%g,%g)\n",spacecraft->faces[1]->r0[0],spacecraft->faces[1]->r0[1],spacecraft->faces[1]->r0[2]);
  printf("rmax=(%g,%g)\n",spacecraft->faces[1]->rmax[0],spacecraft->faces[1]->rmax[1]);
  //face 2
  initialize_face(spacecraft->faces[2],0,SC_BOT_CORNER_3_X,SC_BOT_CORNER_3_Y,SC_BOT_CORNER_4_X,SC_BOT_CORNER_4_Y,
		  1, 0, 1, 0,    3, 0, 0, 0,    8, 1, 1, 0,    9, 1, 1, 0 );   
  printf("r0=(%g,%g,%g)\n",spacecraft->faces[2]->r0[0],spacecraft->faces[2]->r0[1],spacecraft->faces[2]->r0[2]);
  printf("rmax=(%g,%g)\n",spacecraft->faces[2]->rmax[0],spacecraft->faces[2]->rmax[1]);
  //face 3
  initialize_face(spacecraft->faces[3],0,SC_BOT_CORNER_4_X,SC_BOT_CORNER_4_Y,SC_BOT_CORNER_5_X,SC_BOT_CORNER_5_Y,
                  2, 0, 1, 0,    4, 0, 0, 0,   -1, 0, 0, 0,   -1, 0, 0, 0 );//no wrap over diagonal edge.   
  printf("r0=(%g,%g,%g)\n",spacecraft->faces[3]->r0[0],spacecraft->faces[3]->r0[1],spacecraft->faces[3]->r0[2]);
  printf("rmax=(%g,%g)\n",spacecraft->faces[3]->rmax[0],spacecraft->faces[3]->rmax[1]);
  //face 4
  initialize_face(spacecraft->faces[4],0,SC_BOT_CORNER_5_X,SC_BOT_CORNER_5_Y,SC_BOT_CORNER_6_X,SC_BOT_CORNER_6_Y,
		  3, 0, 1, 0,    5, 0, 0, 0,    8, 0, 1, 1,    9, 0, 1, 1 );   
  printf("r0=(%g,%g,%g)\n",spacecraft->faces[4]->r0[0],spacecraft->faces[4]->r0[1],spacecraft->faces[4]->r0[2]);
  printf("rmax=(%g,%g)\n",spacecraft->faces[4]->rmax[0],spacecraft->faces[4]->rmax[1]);
  //face 5
  initialize_face(spacecraft->faces[5],0,SC_BOT_CORNER_6_X,SC_BOT_CORNER_6_Y,SC_BOT_CORNER_7_X,SC_BOT_CORNER_7_Y,
                  4, 0, 1, 0,    6, 0, 0, 0,   -1, 0, 0, 0,   -1, 0, 0, 0 );//no wrap over diagonal edge.   
  printf("r0=(%g,%g,%g)\n",spacecraft->faces[5]->r0[0],spacecraft->faces[5]->r0[1],spacecraft->faces[5]->r0[2]);
  printf("rmax=(%g,%g)\n",spacecraft->faces[5]->rmax[0],spacecraft->faces[5]->rmax[1]);
  //face 6
  initialize_face(spacecraft->faces[6],0,SC_BOT_CORNER_7_X,SC_BOT_CORNER_7_Y,SC_BOT_CORNER_8_X,SC_BOT_CORNER_8_Y,
		  5, 0, 1, 0,    7, 0, 0, 0,    8, 1, 0, 1,    9, 1, 0, 1 );   
  printf("r0=(%g,%g,%g)\n",spacecraft->faces[6]->r0[0],spacecraft->faces[6]->r0[1],spacecraft->faces[6]->r0[2]);
  printf("rmax=(%g,%g)\n",spacecraft->faces[6]->rmax[0],spacecraft->faces[6]->rmax[1]);
  //face 7
  initialize_face(spacecraft->faces[7],0,SC_BOT_CORNER_8_X,SC_BOT_CORNER_8_Y,SC_BOT_CORNER_1_X,SC_BOT_CORNER_1_Y,
		  6, 0, 1, 0,    0, 0, 0, 0,   -1, 0, 0, 0,   -1, 0, 0, 0 );//no wrap over diagonal edge.   
  printf("r0=(%g,%g,%g)\n",spacecraft->faces[7]->r0[0],spacecraft->faces[7]->r0[1],spacecraft->faces[7]->r0[2]);
  printf("rmax=(%g,%g)\n",spacecraft->faces[7]->rmax[0],spacecraft->faces[7]->rmax[1]);
  //face 8, Bottom face:
  initialize_face(spacecraft->faces[8],4,xmin,ymin,xmax,ymin,
		  0, 1, 0, 0,    4, 1, 0, 1,    6, 1, 0, 1,    2, 1, 0, 0 );
  spacecraft->faces[8]->rmax[1]=ymax-ymin;
  set_basis(spacecraft->faces[8]->basis,1,0,0,0,1,0);
  spacecraft->faces[8]->area=2.77331;
  printf("r0=(%g,%g,%g)\n",spacecraft->faces[8]->r0[0],spacecraft->faces[8]->r0[1],spacecraft->faces[8]->r0[2]);
  printf("rmax=(%g,%g)\n",spacecraft->faces[8]->rmax[0],spacecraft->faces[8]->rmax[1]);
  //set cut-outs for bottom
  set_cut( spacecraft->faces[8], 0, SC_BOT_CORNER_8_X, SC_BOT_CORNER_1_Y, 0, 1, 0);
  set_cut( spacecraft->faces[8], 1, SC_BOT_CORNER_3_X, SC_BOT_CORNER_2_Y, 1, 1, 1);
  set_cut( spacecraft->faces[8], 2, SC_BOT_CORNER_4_X, SC_BOT_CORNER_5_Y, 1, 0, 1);
  set_cut( spacecraft->faces[8], 3, SC_BOT_CORNER_7_X, SC_BOT_CORNER_6_Y, 0, 0, 0);  
  //face 9, Top face:
  initialize_face(spacecraft->faces[9],4,xmin,ymin,xmax,ymin,
		  0, 1, 1, 0,    4, 1, 1, 1,    6, 1, 1, 1,    2, 1, 1, 0 );
  spacecraft->faces[9]->r0[2]=SC_H;
  spacecraft->faces[9]->rmax[1]=ymax-ymin;
  set_basis(spacecraft->faces[9]->basis,1,0,0,0,1,0);
  spacecraft->faces[9]->area=2.77331;
  spacecraft->faces[9]->invert=0;
  printf("r0=(%g,%g,%g)\n",spacecraft->faces[9]->r0[0],spacecraft->faces[9]->r0[1],spacecraft->faces[9]->r0[2]);
  printf("rmax=(%g,%g)\n",spacecraft->faces[9]->rmax[0],spacecraft->faces[9]->rmax[1]);
  //set cut-outs for top
  set_cut( spacecraft->faces[9], 0, SC_BOT_CORNER_8_X, SC_BOT_CORNER_1_Y, 0, 1, 0);
  set_cut( spacecraft->faces[9], 1, SC_BOT_CORNER_3_X, SC_BOT_CORNER_2_Y, 1, 1, 1);
  set_cut( spacecraft->faces[9], 2, SC_BOT_CORNER_4_X, SC_BOT_CORNER_5_Y, 1, 0, 1);
  set_cut( spacecraft->faces[9], 3, SC_BOT_CORNER_7_X, SC_BOT_CORNER_6_Y, 0, 0, 0);  
  spacecraft->lpf_area=0;
  for(i=0;i<nfaces;i++)spacecraft->lpf_area+=spacecraft->faces[i]->area;
  printf("lpf_area=%g\n",spacecraft->lpf_area);
}

void set_basis(double ** basis, double e1x,double e1y,double e1z,double e2x,double e2y,double e2z){
  basis[0][0]=e1x;
  basis[0][1]=e1y;
  basis[0][2]=e1z;
  basis[1][0]=e2x;
  basis[1][1]=e2y;
  basis[1][2]=e2z;
  basis[2][0]=basis[0][1]*basis[1][2]-basis[0][2]*basis[1][1];
  basis[2][1]=basis[0][2]*basis[1][0]-basis[0][0]*basis[1][2];
  basis[2][2]=basis[0][0]*basis[1][1]-basis[0][1]*basis[1][0];
  if(1){
    double delta[2][2];
    delta[0][0]=basis[0][0]*basis[0][0]+basis[0][1]*basis[0][1]+basis[0][2]*basis[0][2];
    delta[1][0]=basis[1][0]*basis[0][0]+basis[1][1]*basis[0][1]+basis[1][2]*basis[0][2];
    delta[0][1]=basis[0][0]*basis[1][0]+basis[0][1]*basis[1][1]+basis[0][2]*basis[1][2];
    delta[1][1]=basis[1][0]*basis[1][0]+basis[1][1]*basis[1][1]+basis[1][2]*basis[1][2];
    printf("basis test: delta = { { %g, %g } , { %g, %g } }\n",delta[0][0],delta[0][1],delta[1][0],delta[1][1]);
  }
  printf("normal=(%g,%g,%g)\n",basis[2][0],basis[2][1],basis[2][2]);
}

void initialize_face(struct FaceData *face, int ncuts,double x0,double y0,double x1, double y1,
		     int xdn_face, int xdn_dir, int xdn_end, int xdn_flip,
		     int xup_face, int xup_dir, int xup_end, int xup_flip,
		     int ydn_face, int ydn_dir, int ydn_end, int ydn_flip,
		     int yup_face, int yup_dir, int yup_end, int yup_flip){
  double dx=x1-x0,dy=y1-y0;
  face->rmax=malloc(2*sizeof(double));
  face->rmax[0]=sqrt(dx*dx+dy*dy);
  face->rmax[1]=SC_H;//default ymax=spacecraft height
  face->basis=malloc(3*sizeof(double*));
  face->basis[0]=malloc(3*sizeof(double));
  face->basis[1]=malloc(3*sizeof(double));
  face->basis[2]=malloc(3*sizeof(double));
  set_basis(face->basis,dx/face->rmax[0],dy/face->rmax[0],0,0,0,1);
  face->invert=1; //We are setting outward normals
  face->r0=malloc(3*sizeof(double));
  face->r0[0]=x0;
  face->r0[1]=y0;
  face->r0[2]=0;
  face->area=face->rmax[0]*face->rmax[1];

  face->xdn_face=xdn_face;
  face->xdn_dir=xdn_dir;
  face->xdn_end=xdn_end;
  face->xdn_flip=xdn_flip;

  face->xup_face=xup_face;
  face->xup_dir=xup_dir;
  face->xup_end=xup_end;
  face->xup_flip=xup_flip;

  face->ydn_face=ydn_face;
  face->ydn_dir=ydn_dir;
  face->ydn_end=ydn_end;
  face->ydn_flip=ydn_flip;

  face->yup_face=yup_face;
  face->yup_dir=yup_dir;
  face->yup_end=yup_end;
  face->yup_flip=yup_flip;

  face->ncuts=ncuts;
  face->cut_xc=malloc(ncuts*sizeof(double));
  face->cut_yc=malloc(ncuts*sizeof(double));
  face->cut_left=malloc(ncuts*sizeof(int));
  face->cut_top=malloc(ncuts*sizeof(int));
  face->cut_above=malloc(ncuts*sizeof(int));
} 

void set_cut(struct FaceData *face, int idx, double xc, double yc, int above, int left, int top){
  face->cut_xc[idx]    = xc - face->r0[0];
  face->cut_yc[idx]    = yc - face->r0[1];
  printf("setting cut: xc=%g, yc=%g, a,l,t=%i,%i,%i\n",face->cut_xc[idx],face->cut_yc[idx],above,left,top);
  face->cut_above[idx] = above;
  face->cut_left[idx]  = left;
  face->cut_top[idx]   = top;
}

void MomentOfInertia(double ***I)
{

  I[0][0][0] = EOM_SC_IH1_XX;
  I[0][0][1] = EOM_SC_IH1_XY;
  I[0][0][2] = EOM_SC_IH1_XZ;
  
  I[0][1][0] = EOM_SC_IH1_YX;
  I[0][1][1] = EOM_SC_IH1_YY;
  I[0][1][2] = EOM_SC_IH1_YZ;
  
  I[0][2][0] = EOM_SC_IH1_ZX;
  I[0][2][1] = EOM_SC_IH1_ZY;
  I[0][2][2] = EOM_SC_IH1_ZZ;

  I[1][0][0] = EOM_SC_IH2_XX;
  I[1][0][1] = EOM_SC_IH2_XY;
  I[1][0][2] = EOM_SC_IH2_XZ;

  I[1][1][0] = EOM_SC_IH2_YX;
  I[1][1][1] = EOM_SC_IH2_YY;
  I[1][1][2] = EOM_SC_IH2_YZ;

  I[1][2][0] = EOM_SC_IH2_ZX;
  I[1][2][1] = EOM_SC_IH2_ZY;
  I[1][2][2] = EOM_SC_IH2_ZZ;

}


//***********************************************
//   Routines handling Spacecraft faces
//***********************************************

void adjust_face(struct Spacecraft *lpf, int *iface, double *r, gsl_rng *seed)
{
  if(*iface<0)return;
  int i;
  //spacecraft dimensions
  int last=-1;
  int nextlast=-1;
  int ixy=gsl_rng_uniform(seed)*2;//randomly choose 0/1 to avoid asymmetry in jumps that wrap around 2 edges near a cut.

  //First we must check if the face coordinate point falls in (or beyond) one of the cut regions.
  struct FaceData *fd=lpf->faces[*iface];
  //printf("adjust: [%i](%g,%g), rmax=(%g,%g)\n",*iface,r[0],r[1],fd->rmax[0],fd->rmax[1]);
  double xmax=fd->rmax[0];
  double ymax=fd->rmax[1];
  double x=r[0],y=r[1]; //define the coordinate point (x,y) confined to the rectangle for the cut test.
  if(x<0)x=0; 
  if(x>xmax)x=xmax;
  if(y<0)y=0;
  if(y>ymax)y=ymax;
  for(i=0;i<fd->ncuts;i++){ //loop over the specified cuts.
    if(test_cut(fd,i,x,y)==1){
      *iface=-1;
      return;
    }
  }//end of cuts
  
  //Next roll over to another face if needed.
  while(*iface>=0&&(*iface!=last||last!=nextlast)){//require no change twice consequtively
    nextlast=last;
    last=*iface;
    wrap_face(lpf,ixy,iface,r);
    ixy=1-ixy;
  }
}


void wrap_face(struct Spacecraft *lpf, int dir, int *iface, double *r){
  struct FaceData *face=lpf->faces[*iface];
  struct FaceData *oface;
  //printf("wrap?: [%i](%g,%g), rmax[%i]=%g\n",*iface,r[0],r[1],dir,face->rmax[dir]);
  double newr[2];
  int odir=-1,end=0,ioface=-1,flip=0;
  if(dir==0){//x
    if(r[dir]<0){
      ioface=face->xdn_face;
      odir=face->xdn_dir;
      flip=face->xdn_flip;
      end=face->xdn_end;
      r[dir]*=-1;
    } else if (r[dir]>face->rmax[dir]) {
      ioface=face->xup_face;
      odir=face->xup_dir;
      flip=face->xup_flip;
      end=face->xup_end;
      r[dir]-=face->rmax[dir];
    }
  } else { //y
    if(r[dir]<0){
      ioface=face->ydn_face;
      odir=face->ydn_dir;
      flip=face->ydn_flip;
      end=face->ydn_end;
      r[dir]*=-1;
    } else if (r[dir]>face->rmax[dir]) {
      ioface=face->yup_face;
      odir=face->yup_dir;
      flip=face->yup_flip;
      end=face->yup_end;
      r[dir]-=face->rmax[dir];
    }
  }
  if(odir<0)return;
  //printf("wrapping from face %i to face %i\n",*iface, ioface);
  if(ioface>=0){
    //printf("start: [%i](%g,%g), rmax[%i]=%g\n",*iface,r[0],r[1],dir,face->rmax[dir]);
    int i;
    oface=lpf->faces[ioface];
    //compute the zero-offset for the transf to the new face's frame
    double r0i=0;for(i=0;i<3;i++)r0i+=face->basis[1-dir][i]*face->r0[i];
    double r0o=0;for(i=0;i<3;i++)r0o+=oface->basis[1-odir][i]*oface->r0[i];
    if(flip==1)r0o*=-1;
    double dr0=r0i-r0o;;
    //tranfsform to new face frame
    newr[odir]=r[dir];
    if(end)newr[odir]=oface->rmax[odir]-newr[odir];
    newr[1-odir]=r[1-dir];
    newr[1-odir]+=dr0;
    if(flip==1)newr[1-odir]*=-1;
    r[0]=newr[0];
    r[1]=newr[1];
    //printf("   to: [%i](%g,%g), rmax[%i]=%g, end=%i, flip=%i\n",ioface,newr[0],newr[1],odir,oface->rmax[odir],end,flip);
  } 
  *iface=ioface;    
}

//Returns 3d body-frame vector from in-face coordinate
void face2body(struct Spacecraft *lpf, int iface, double *rface, double *rbody){
  int i;
  if(iface<0)return;
  struct FaceData *face=lpf->faces[iface];
  for(i=0;i<3;i++) rbody[i] = face->r0[i] + rface[0]*face->basis[0][i] + rface[1]*face->basis[1][i];
  //printf("r0=(%g,%g,%g)\n",face->r0[0],face->r0[1],face->r0[2]);
  //printf(" converted face -> body: [%i](%g,%g) -> (%g,%g,%g)\n",iface,rface[0],rface[1],rbody[0],rbody[1],rbody[2]);
  //double rref[3];
  //printf("r0=(%g,%g,%g)\n",face->r0[0],face->r0[1],face->r0[2]);
  //for(i=0;i<3;i++) rref[i] = rbody[i] - face->r0[i];
  //double zmiss=face->basis[2][0]*rref[0]+face->basis[2][1]*rref[1]+face->basis[2][2]*rref[2];
  //if(fabs(zmiss)>1e-12)printf("* * * * * Off-face! zmiss=%g\n",zmiss);
}

//Returns 3d body-frame vector from in-face coordinate
void body2face(struct Spacecraft *lpf, int iface, double *rbody, double *rface){
  int i,j;
  if(iface<0)return;
  struct FaceData *face=lpf->faces[iface];
  double rref[3];
  //First reference face corner
  //printf("r0=(%g,%g,%g)\n",face->r0[0],face->r0[1],face->r0[2]);
  for(i=0;i<3;i++) rref[i] = rbody[i] - face->r0[i];
  //double zmiss=face->basis[2][0]*rref[0]+face->basis[2][1]*rref[1]+face->basis[2][2]*rref[2];
  //printf("e3.rref=%g\n",zmiss);
  //if(fabs(zmiss)>1e-12)printf("***** Off-face! zmiss=%g\n",zmiss);
  for(j=0;j<2;j++){
    rface[j]=0;
    for(i=0;i<3;i++) rface[j] += face->basis[j][i]*rref[i];
  }
  //printf(" converted body -> face: (%g,%g,%g) -> [%i](%g,%g)\n",rbody[0],rbody[1],rbody[2],iface,rface[0],rface[1]); 
}


///JBG: Tyson probably already has something that could be made to work, but this is indep. and optimized for this computation
void face_sky_to_body_sky(struct Spacecraft *lpf, int face, double costhetaF, double phiF, double* costheta, double* phi)
{
  double *n=malloc(3*sizeof(double));
  get_normal_lpf(lpf,n,face);//We assume this is unit-normalized

  //There are two conventions that are important here:
  //We assume either:  sky-costheta and sky-phi point *toward* the source and the face normals are *outward*
  //              or:  sky-costheta and sky-phi point *away from* the source and the face normals are *inward*
  
  //This computation should realize, beginning with some unit vector in the face frame, (with the "x-side" assumed to
  //be moving toward the bottom [if that matters]) first a rotation to shift the face normal to the right inclination
  //in the x-z plane, then a rotation of the x-y plane.
  double rhoF=sqrt(1-costhetaF*costhetaF);
  double rhoN=sqrt(1-n[2]*n[2]);
  double xF=cos(phiF)*rhoF;
  double yF=sin(phiF)*rhoF;
  double coeff=n[2]*xF+rhoN*costhetaF;
  if(fabs(n[2])<1)
    *phi      = atan2(  n[1]*coeff + n[0]*yF , n[0]*coeff - n[1]*yF );
  else
    *phi=phiF*n[2];
  while(*phi<0)*phi+=2*M_PI;
  *costheta = -rhoN*xF + n[2]*costhetaF;
  free(n);
}

void get_normal_lpf(struct Spacecraft *lpf, double *n, int face){
  //Whether this is inward-facing or outward facing depends on how inverts are set in initialization
  if(face<0)return;
  struct FaceData *fc = lpf->faces[face];
  int i;
  int sign=1-2*fc->invert;
  for(i=0;i<3;i++)n[i]=fc->basis[2][i]*sign;
}

double incidence(struct Spacecraft *lpf,int iface, double cth, double phi){
  double normvec[3];
  get_normal_lpf(lpf,normvec,iface);
  //printf("normal[%i] = (%g,%g,%g)\n",iface,normvec[0],normvec[1],normvec[2]);
  double rho=sqrt(1-cth*cth);
  double nk = rho*(normvec[0]*cos(phi) + normvec[1]*sin(phi)) + normvec[2]*cth;
  return nk;
}  

//Handling face cuts

void write_faces(FILE *out, struct Spacecraft *lpf){
  int n,i,j,k;
  for(n=0;n<lpf->nfaces;n++){//loop over faces
    struct FaceData *f=lpf->faces[n];
    int np=4;
    double p[lpf->nfaces+2*f->ncuts][2];
    //first add in the corner points (in face frame)
    p[0][0]=0;
    p[0][1]=0;
    p[1][0]=f->rmax[0];
    p[1][1]=0;
    p[2][0]=f->rmax[0];
    p[2][1]=f->rmax[1];
    p[3][0]=0;
    p[3][1]=f->rmax[1];
    
    //The trick is handling the cuts, each of which adds two points, but may eliminate other points
    for(k=0;k<f->ncuts;k++){
      //Now loop over the corners to see which (if any are cut);
      //For a convex shape the cut-corners must be consequtive
      int iex0=-1,nex=0,cutlast=0; //exclude nex corners beginning with corner iex0;
      for(i=0;i<np;i++){
	if(test_cut(f,k,p[i][0],p[i][1])==1){//exclude corner
	  if(cutlast==0)iex0=i;
	  nex++;
	  cutlast=1;
	  //printf("face %i, cut %i: cutting point[%i]=(%g,%g)\n",n,k,i,p[i][0],p[i][1]);
	} else cutlast=0;
      }
      double newp1[2],newp2[2];
      int i0,i1;
      if(nex==0)continue;//nothing to cut.
      //The cut intersects the line from p[iex0-1] to p[iex0] to form a new point
      i1=iex0;
      i0=i1-1; if(i0<0)i0+=np;//wrap
      get_cut_intersection(f,k,p[i0][0],p[i0][1],p[i1][0],p[i1][1],&newp1[0],&newp1[1]);
      double x0,y0,x1,y1;
      get_cut_line(f,k,&x0,&y0,&x1,&y1);
      //printf("face %i, cut %i: first intersection cut=(%g,%g)-(%g,%g)\n",n,k,x0,y0,x1,y1);
      //printf("face %i, cut %i: first intersection edge=[%i](%g,%g)-[%i](%g,%g)\n",n,k,i0,p[i0][0],p[i0][1],i1,p[i1][0],p[i1][1]);
      //printf("face %i, cut %i: first intersection point=(%g,%g)\n",n,k,newp1[0],newp1[1]);
      //The cut intersects the line from p[iex0+nex-1] to p[iex0+nex] to form a new point
      i1=iex0+nex; if(i1>=np)i1-=np;//wrap
      i0=i1-1;if(i0<0)i0+=np;//wrap
      get_cut_intersection(f,k,p[i0][0],p[i0][1],p[i1][0],p[i1][1],&newp2[0],&newp2[1]);
      if(nex==1){//need to insert a point
	for(i=np;i>iex0+nex;i--)for(j=0;j<2;j++)p[i][j]=p[i-1][j];
	np++;
      } else if(nex>2){//need to delete points
	for(i=iex0+nex;i<np;i++)for(j=0;j<2;j++)p[i+2-nex][j]=p[i][j];
	np-=(nex-2);
      }
      //insert new points
      for(j=0;j<2;j++)p[iex0][j]=newp1[j];
      for(j=0;j<2;j++)p[iex0+1][j]=newp2[j];
    }//end of cuts

    //Write the results out

    //First write a comment with basic face info:
    double rn[3];
    get_normal_lpf(lpf,rn,n);
    fprintf(out,"#Face %i: normal=(%g,%g,%g)\n",n,rn[0],rn[1],rn[2]);
    /*
    double r[3];
    for(i=0;i<np;i++){
      face2body(lpf,n,p[i],r);
      fprintf(out,"%g, %g, %g\n",r[0],r[1],r[2]);
    }
    face2body(lpf,n,p[0],r);
    */
    for(i=0;i<(np-1)/2;i++){
      double r00[3];
      double r01[3];
      double r10[3];
      double r11[3];
      int i0=i;
      int i0opp=np-1-i0;
      int i1=i0+1;
      int i1opp=np-1-i1;
      face2body(lpf,n,p[i0],r00);
      face2body(lpf,n,p[i0opp],r01);
      face2body(lpf,n,p[i1],r10);
      face2body(lpf,n,p[i1opp],r11);
      int ngrid=5;
      for(j=0;j<ngrid+1;j++){
	for(k=0;k<ngrid+1;k++){
	  double dj=j/(double)ngrid;
	  double dk=k/(double)ngrid;
	  fprintf(out,"%g, %g, %g\n",
		  (r00[0]*(1-dj)+r01[0]*dj)*(1-dk)+(r10[0]*(1-dj)+r11[0]*dj)*dk,
		  (r00[1]*(1-dj)+r01[1]*dj)*(1-dk)+(r10[1]*(1-dj)+r11[1]*dj)*dk,
		  (r00[2]*(1-dj)+r01[2]*dj)*(1-dk)+(r10[2]*(1-dj)+r11[2]*dj)*dk);
	}
	fprintf(out,"\n");
      }
      fprintf(out,"\n");
    }
  }//end of loop over faces
}

//***********************************************
//   Routines handling Spacecraft face cuts
//***********************************************


//test whether point (x,y) is excluded by the cut. Returns 1 if excluded; Result undef. if (x,y) is outside rect.
int test_cut(struct FaceData *fd,int icut, double x, double y){
  double x0=0,x1=0,y0=0,y1=0;
  int sign=1-2*fd->cut_above[icut]; // i.e. -1 if cut_above==1, 1 if cut_above==0
  get_cut_line(fd,icut,&x0,&y0,&x1,&y1);
  //cut line is 0=x0+(x-x0)*(y1-y0)/(x1-x0)-y
  if( sign * ( y0+(x-x0)*(y1-y0)/(x1-x0)-y ) > 0 )return 1; //do cut!
  return 0;
}

void get_cut_line(struct FaceData *fd,int icut, double *x0, double *y0, double *x1, double *y1){
  double xmax=fd->rmax[0];
  double ymax=fd->rmax[1];
  double xc=fd->cut_xc[icut];
  double yc=fd->cut_yc[icut];
  int left=fd->cut_left[icut];
  int top=fd->cut_top[icut];
  //printf("rmax=(%g,%g), xc=(%g,%g), left=%i, top=%i\n",xmax,ymax,xc,yc,left,top);
  if(left==1){//cut defined from ( 0, yc ) to ( xc, 0 ) or ( xc, ymax ) if top==1
    *x0=0;
    *y0=yc;
    *x1=xc;
    *y1=ymax*top;
    //printf(" leftt: pc0=(%g,%g), pc1=(%g,%g)\n",0.0,yc,xc,ymax*top);
  }  else {              //cut defined from ( xc, 0 ) or ( xc, ymax ) [if top==1] to ( xmax, yc ).
    *x0=xc;
    *y0=ymax*top;
    *x1=xmax;
    *y1=yc;
    //printf("right: pc0=(%g,%g), pc1=(%g,%g)\n",xc,xmax*top,xmax,yc);
  }
  //printf("cut: pc0=(%g,%g), pc1=(%g,%g)\n",*x0,*y0,*x1,*y1);
  //cut line is 0=x0+(x-x0)*(y1-y0)/(x1-x0)-y
}

//find where cut meets a line
void get_cut_intersection(struct FaceData *fd,int icut, double x0, double y0, double x1, double y1, double *x, double *y){
  double xc0=0,yc0=0,xc1=0,yc1=0;
  get_cut_line(fd,icut,&xc0,&yc0,&xc1,&yc1);
  //printf("p0=(%g,%g), p1=(%g,%g)\n",x0,y0,x1,y1);
  //printf("pc0=(%g,%g), pc1=(%g,%g)\n",xc0,yc0,xc1,yc1);
  double dx=x1-x0;
  double dy=y1-y0;
  double dxc=xc1-xc0;
  double dyc=yc1-yc0;
  double Cxx=dx*dxc;
  double Cxy=dx*dyc;
  double Cyx=dy*dxc;
  double Cyy=dy*dyc;
  //printf("dx=%g, dy=%g, dxc=%g, dyc=%g\n",dx,dy,dxc,dyc);
  //printf("Cxx=%g, Cxy=%g, Cyx=%g, Cyy=%g\n",Cxx,Cxy,Cyx,Cyy);
  *x = ( (yc0-y0)*Cxx - xc0*Cxy + x0*Cyx ) / ( Cyx - Cxy );
  *y = ( (xc0-x0)*Cyy - yc0*Cyx + y0*Cxy ) / ( Cxy - Cyx );

  //check:
  //  Delta = Cxy-Cyx
  // (x-x0)*-Delta = (yc0-y0)*Cxx - xc0*Cxy + x0*Cyx + x0*(Cxy - Cyx)
  //               = (yc0-y0)*Cxx - (xc0-x0)*Cxy 
  // (y-y0)*Delta =  (xc0-x0)*Cyy - yc0*Cyx + y0*Cxy - y0*(Cxy - Cyx)
  // (y-y0)*Delta =  (xc0-x0)*Cyy - (yc0-y0)*Cyx
  // (y-y0)*dx0*Delta = (xc0-x0)*dx0*Cyy - (yc0-y0)*dx0*Cyx
  //                  = (xc0-x0)*dy0*Cxy - (yc0-y0)*dy0*Cxx
  //                  = dy0*((xc0-x0)*Cxy - (yc0-y0)*Cxx )
  //                  = dy0*(x-x0)*Delta                               check!
}





