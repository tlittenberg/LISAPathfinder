// LPF.c Routines for defining basic LISA Pathfinder structural properties
// Written by John Baker 2015

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>


//*******************************************
//             Data Structures
//*******************************************

struct FaceData
{
  //Spacecraft face parameters
  double *rmax;       //(x,y) extent of rectangular panel
  double **basis;
  double *r0;
  double area;
  int invert;        //set to 1 if XxY is ingoing.
  int xdn_face;           //index of the connected face in the neg-x direction
  int xdn_dir;       //the idirection (0=x,1=y) the movement in x dir changes to on next panel.
  int xdn_end;       //connected to low (0) or high end (1) side of "dir" direction on other panel 
  int xdn_flip;      //if 1 the tranverse direction increases in opposite direction.
  double xdn0;       //y=0 maps to not-dir coord has specified value;
  int xup_face;           //index of the connected face in the pos-x direction
  int xup_dir;       //the idirection (0=x,1=y) the movement in x dir changes to on next panel.
  int xup_end;       //connected to low (0) or high end (1) side of "dir" direction on other panel 
  int xup_flip;      //if 1 the tranverse direction increases in opposite direction.
  double xup0;       //y=0 maps to not-dir coord has specified value;
  int ydn_face;           //index of the connected face in the neg-y direction
  int ydn_dir;       //the idirection (0=x,1=y) the movement in y dir changes to on next panel.
  int ydn_end;       //connected to low (0) or high end (1) side of "dir" direction on other panel 
  int ydn_flip;      //if 1 the tranverse direction increases in opposite direction.
  double ydn0;       //y=0 maps to not-dir coord has specified value;
  int yup_face;           //index of the connected face in the pos-y direction
  int yup_dir;       //the idirection (0=x,1=y) the movement in y dir changes to on next panel.
  int yup_end;       //connected to low (0) or high end (1) side of "dir" direction on other panel 
  int yup_flip;      //if 1 the tranverse direction increases in opposite direction.
  double yup0;       //y=0 maps to not-dir coord has specified value;
  int ncuts;         //number of specified cutouts for this panel
  double *cut_xc;    //one intercept at (xc,0)
  double *cut_yc;    //other intercept at (xmax,yc)
  int *cut_above;    //1 to exclude the region above the cut line, 0 to exclude below.
  int *cut_left;     //1 to have yc define the left intercept
  int *cut_top ;     //1 to have xc define the top intercept
};

struct Spacecraft
{
  /* Spacecraft parameters */
  double **R;        //position of proof masses
  double **x;        //position of spacecraft corners
  double ***I;       //moment of inertia (ultimately a tensor)
  double ***invI;    //inverse moment of inertia (ultimately a tensor)
  double M;          //mass
  double W;          //spacecraft width  (x-direction, m)
  double D;          //spacecraft depth  (y-direction, m)
  double H;          //spacecraft height (z-direciton, m)
  int nfaces;
  struct FaceData **faces;   //info about the spaces.
  double lpf_area;
};

//*******************************************
//           Function Prototypes
//*******************************************


//Spacecraft struct initialization

void initialize_spacecraft(struct Spacecraft *spacecraft);

void set_basis(double ** basis, double e1x,double e1y,double e1z,double e2x,double e2y,double e2z);

void initialize_face(struct FaceData *face, int ncuts, double x0, double y0, double x1, double y1,
		     int xdn_face, int xdn_dir, int xdn_end, int xdn_flip,
		     int xup_face, int xup_dir, int xup_end, int xup_flip,
		     int ydn_face, int ydn_dir, int ydn_end, int ydn_flip,
		     int yup_face, int yup_dir, int yup_end, int yup_flip);

void set_cut(struct FaceData *face, int idx, double xc, double yc, int above, int left, int top);

void MomentOfInertia(double ***I);
  
//Spacecraft struct utilization

void adjust_face(struct Spacecraft *lpf, int *face, double *r, gsl_rng *seed);

void wrap_face(struct Spacecraft *lpf, int dir, int *iface, double *r);

void face2body(struct Spacecraft *lpf, int iface, double *rface, double *rbody);

void body2face(struct Spacecraft *lpf, int iface, double *rbody, double *rface);

void face_sky_to_body_sky(struct Spacecraft *lpf, int face, double costhetaF, double phiF, double* costheta, double* phi);

void get_normal_lpf(struct Spacecraft *lpf, double *n, int face);

double incidence(struct Spacecraft *lpf,int iface, double cth, double phi);

void write_faces(FILE *out, struct Spacecraft *lpf);

int test_cut(struct FaceData *fd,int icut, double x, double y);

void get_cut_line(struct FaceData *fd,int icut, double *x0, double *y0, double *x1, double *y1);

void get_cut_intersection(struct FaceData *fd,int icut, double x0, double y0, double x1, double y1, double *x, double *y);


