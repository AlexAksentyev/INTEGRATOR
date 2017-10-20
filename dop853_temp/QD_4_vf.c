/*  Vector field function and events for Dopri853 integrator.
  This code was automatically generated by PyDSTool, but may be modified by hand. */

#include <math.h>
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "events.h"
#include "maxmin.h"
#include "signum.h"
#include "vfield.h"

extern double *gICs;
extern double **gBds;
extern double globalt0;

static double pi = 3.1415926535897931;

double signum(double x)
{
  if (x<0) {
    return -1;
  }
  else if (x==0) {
    return 0;
  }
  else if (x>0) {
    return 1;
  }
  else {
    /* must be that x is Not-a-Number */
    return x;
  }
}


/* Variable, aux variable, parameter, and input definitions: */ 
#define Curve	p_[0]
#define G	p_[1]
#define KinEn0	p_[2]
#define Length	p_[3]
#define Mass0	p_[4]
#define clight	p_[5]
#define m0	p_[6]
#define offset	p_[7]
#define q	p_[8]
#define H	Y_[0]
#define Ss	Y_[1]
#define Sx	Y_[2]
#define Sy	Y_[3]
#define dK	Y_[4]
#define px	Y_[5]
#define py	Y_[6]
#define s	Y_[7]
#define start	Y_[8]
#define ts	Y_[9]
#define x	Y_[10]
#define y	Y_[11]

double NaN_event(unsigned n_, double t, double *Y_, double *p_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_);
double passto19(unsigned n_, double t, double *Y_, double *p_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_);

double Bs(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_);
double Bx(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_);
double By(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_);
double Es(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_);
double Ex(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_);
double Ey(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_);
double KinEn(double __dK__, double *p_, double *wk_, double *xv_);
double Lbeta(double __dK__, double *p_, double *wk_, double *xv_);
double Lgamma(double __dK__, double *p_, double *wk_, double *xv_);
double Pc(double __dK__, double *p_, double *wk_, double *xv_);
double __maxof2(double e1_, double e2_, double *p_, double *wk_, double *xv_);
double __maxof3(double e1_, double e2_, double e3_, double *p_, double *wk_, double *xv_);
double __maxof4(double e1_, double e2_, double e3_, double e4_, double *p_, double *wk_, double *xv_);
double __minof2(double e1_, double e2_, double *p_, double *wk_, double *xv_);
double __minof3(double e1_, double e2_, double e3_, double *p_, double *wk_, double *xv_);
double __minof4(double e1_, double e2_, double e3_, double e4_, double *p_, double *wk_, double *xv_);
double __rhs_if(int cond_, double e1_, double e2_, double *p_, double *wk_, double *xv_);
double getbound(char *name, int which_bd, double *p_, double *wk_, double *xv_);
double globalindepvar(double t, double *p_, double *wk_, double *xv_);
double initcond(char *varname, double *p_, double *wk_, double *xv_);
int getindex(char *name, double *p_, double *wk_, double *xv_);
int heav(double x_, double *p_, double *wk_, double *xv_);

int N_EVENTS = 2;
void assignEvents(EvFunType *events){
 events[0] = &NaN_event;
events[1] = &passto19;

}

void auxvars(unsigned, unsigned, double, double*, double*, double*, unsigned, double*, unsigned, double*);
void jacobian(unsigned, unsigned, double, double*, double*, double**, unsigned, double*, unsigned, double*);
void jacobianParam(unsigned, unsigned, double, double*, double*, double**, unsigned, double*, unsigned, double*);
int N_AUXVARS = 0;


int N_EXTINPUTS = 0;


void vfieldfunc(unsigned n_, unsigned np_, double t, double *Y_, double *p_, double *f_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_){
/* reused term definitions */
double v0_Bs = Bs(x,y,ts,px,py,dK,H,s,start,Sx,Sy,Ss, p_, wk_, xv_);
double v0_Bx = Bx(x,y,ts,px,py,dK,H,s,start,Sx,Sy,Ss, p_, wk_, xv_);
double v0_By = By(x,y,ts,px,py,dK,H,s,start,Sx,Sy,Ss, p_, wk_, xv_);
double v0_Es = Es(x,y,ts,px,py,dK,H,s,start,Sx,Sy,Ss, p_, wk_, xv_);
double v0_Ex = Ex(x,y,ts,px,py,dK,H,s,start,Sx,Sy,Ss, p_, wk_, xv_);
double v0_Ey = Ey(x,y,ts,px,py,dK,H,s,start,Sx,Sy,Ss, p_, wk_, xv_);
double v0_Lbeta = Lbeta(dK, p_, wk_, xv_);
double v0_Lgamma = Lgamma(dK, p_, wk_, xv_);
double v0_P0c = Pc(0, p_, wk_, xv_);
double v0_Pc = Pc(dK, p_, wk_, xv_);
double v0_hs = (1+Curve*x);
double v1_V = v0_Lbeta*clight;
double v1_Px = v0_P0c*px;
double v1_Py = v0_P0c*py;
double v2_Ps = sqrt(pow(v0_Pc,2)-pow(v1_Px,2)-pow(v1_Py,2));
double v2_Vx = v1_V*v1_Px/v0_Pc;
double v2_Vy = v1_V*v1_Py/v0_Pc;
double v3_Vs = v1_V*v2_Ps/v0_Pc;
double v3_Xp = v0_hs*v0_P0c*px/v2_Ps;
double v3_Yp = v0_hs*v0_P0c*py/v2_Ps;
double v3_Hp = v0_hs*v0_Pc/v2_Ps;
double v4_Tp = v3_Hp/v1_V;
double v5_t6 = v4_Tp*(q/(v0_Lgamma*m0*m0*clight*clight))*(G+1/(1+v0_Lgamma));
double v5_sp1 = v4_Tp*(-q/(v0_Lgamma*m0))*(1+G*v0_Lgamma);
double v5_sp2 = v4_Tp*(q/(v0_Lgamma*m0*m0*m0*clight*clight))*(G/(1+v0_Lgamma))*(v1_Px*v0_Bx+v1_Py*v0_By+v2_Ps*v0_Bs);
double v6_Sxp = Curve*Ss+v5_t6*((v2_Ps*v0_Ex-v1_Px*v0_Es)*Ss-(v1_Px*v0_Ey-v1_Py*v0_Ex)*Sy)+(v5_sp1*v0_By+v5_sp2*v1_Py)*Ss-(v5_sp1*v0_Bs+v5_sp2*v2_Ps)*Sy;
double v6_Syp = v5_t6*((v1_Px*v0_Ey-v1_Py*v0_Ex)*Sx-(v1_Py*v0_Es-v2_Ps*v0_Ey)*Ss)+(v5_sp1*v0_Bs+v5_sp2*v2_Ps)*Sx-(v5_sp1*v0_Bx+v5_sp2*v1_Px)*Ss;
double v6_Ssp = (-1)*Curve*Sx+v5_t6*((v1_Py*v0_Es-v2_Ps*v0_Ey)*Sy-(v2_Ps*v0_Ex-v1_Px*v0_Es)*Sx)+(v5_sp1*v0_Bx+v5_sp2*v1_Px)*Sy-(v5_sp1*v0_By+v5_sp2*v1_Py)*Sx;

f_[0] = v3_Hp;
f_[1] = v6_Ssp;
f_[2] = v6_Sxp;
f_[3] = v6_Syp;
f_[4] = ((v0_Ex*v0_hs*v1_Px/v2_Ps)+(v0_Ey*v0_hs*v1_Py/v2_Ps)+v0_Es)*1e-6/KinEn0;
f_[5] = ((q*(v0_Ex+((v2_Vy*v0_Bs)-(v0_By*v3_Vs)))*v4_Tp)+Curve*v2_Ps)/v0_P0c;
f_[6] = (q*(v0_Ey+((v3_Vs*v0_Bx)-(v0_Bs*v2_Vx)))*v4_Tp)/v0_P0c;
f_[7] = 1;
f_[8] = 0;
f_[9] = v4_Tp;
f_[10] = v3_Xp;
f_[11] = v3_Yp;

}




double Bs(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_) {


return 0 ;

}


double Bx(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_) {


return -8.2*(-__y__);

}


double By(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_) {


return -8.2*(-__x__);

}


double Es(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_) {


return 0 ;

}


double Ex(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_) {


return 0 ;

}


double Ey(double __x__, double __y__, double __ts__, double __px__, double __py__, double __dK__, double __H__, double __s__, double __start__, double __Sx__, double __Sy__, double __Ss__, double *p_, double *wk_, double *xv_) {


return 0 ;

}


double KinEn(double __dK__, double *p_, double *wk_, double *xv_) {


return KinEn0*(1+__dK__);

}


double Lbeta(double __dK__, double *p_, double *wk_, double *xv_) {


return sqrt(pow(Lgamma(__dK__, p_, wk_, xv_),2)-1)/Lgamma(__dK__, p_, wk_, xv_);

}


double Lgamma(double __dK__, double *p_, double *wk_, double *xv_) {


return KinEn(__dK__, p_, wk_, xv_)/Mass0+1 ;

}


double Pc(double __dK__, double *p_, double *wk_, double *xv_) {


return sqrt(pow(Mass0+KinEn(__dK__, p_, wk_, xv_),2)-pow(Mass0,2));

}


double __maxof2(double e1_, double e2_, double *p_, double *wk_, double *xv_) {
if (e1_ > e2_) {return e1_;} else {return e2_;};
}


double __maxof3(double e1_, double e2_, double e3_, double *p_, double *wk_, double *xv_) {
double temp_;
if (e1_ > e2_) {temp_ = e1_;} else {temp_ = e2_;};
if (e3_ > temp_) {return e3_;} else {return temp_;};
}


double __maxof4(double e1_, double e2_, double e3_, double e4_, double *p_, double *wk_, double *xv_) {
double temp_;
if (e1_ > e2_) {temp_ = e1_;} else {temp_ = e2_;};
if (e3_ > temp_) {temp_ = e3_;};
if (e4_ > temp_) {return e4_;} else {return temp_;};
}


double __minof2(double e1_, double e2_, double *p_, double *wk_, double *xv_) {
if (e1_ < e2_) {return e1_;} else {return e2_;};
}


double __minof3(double e1_, double e2_, double e3_, double *p_, double *wk_, double *xv_) {
double temp_;
if (e1_ < e2_) {temp_ = e1_;} else {temp_ = e2_;};
if (e3_ < temp_) {return e3_;} else {return temp_;};
}


double __minof4(double e1_, double e2_, double e3_, double e4_, double *p_, double *wk_, double *xv_) {
double temp_;
if (e1_ < e2_) {temp_ = e1_;} else {temp_ = e2_;};
if (e3_ < temp_) {temp_ = e3_;};
if (e4_ < temp_) {return e4_;} else {return temp_;};
}


double __rhs_if(int cond_, double e1_, double e2_, double *p_, double *wk_, double *xv_) {
  if (cond_) {return e1_;} else {return e2_;};
}


double getbound(char *name, int which_bd, double *p_, double *wk_, double *xv_) {
  return gBds[which_bd][getindex(name, p_, wk_, xv_)];
}


double globalindepvar(double t, double *p_, double *wk_, double *xv_) {
  return globalt0+t;
}


double initcond(char *varname, double *p_, double *wk_, double *xv_) {

  if (strcmp(varname, "H")==0)
	return gICs[0];
  else if (strcmp(varname, "Ss")==0)
	return gICs[1];
  else if (strcmp(varname, "Sx")==0)
	return gICs[2];
  else if (strcmp(varname, "Sy")==0)
	return gICs[3];
  else if (strcmp(varname, "dK")==0)
	return gICs[4];
  else if (strcmp(varname, "px")==0)
	return gICs[5];
  else if (strcmp(varname, "py")==0)
	return gICs[6];
  else if (strcmp(varname, "s")==0)
	return gICs[7];
  else if (strcmp(varname, "start")==0)
	return gICs[8];
  else if (strcmp(varname, "ts")==0)
	return gICs[9];
  else if (strcmp(varname, "x")==0)
	return gICs[10];
  else if (strcmp(varname, "y")==0)
	return gICs[11];
  else {
	fprintf(stderr, "Invalid variable name %s for initcond call\n", varname);
	return 0.0/0.0;
	}
}


int getindex(char *name, double *p_, double *wk_, double *xv_) {

  if (strcmp(name, "H")==0)
	return 0;
  else if (strcmp(name, "Ss")==0)
	return 1;
  else if (strcmp(name, "Sx")==0)
	return 2;
  else if (strcmp(name, "Sy")==0)
	return 3;
  else if (strcmp(name, "dK")==0)
	return 4;
  else if (strcmp(name, "px")==0)
	return 5;
  else if (strcmp(name, "py")==0)
	return 6;
  else if (strcmp(name, "s")==0)
	return 7;
  else if (strcmp(name, "start")==0)
	return 8;
  else if (strcmp(name, "ts")==0)
	return 9;
  else if (strcmp(name, "x")==0)
	return 10;
  else if (strcmp(name, "y")==0)
	return 11;
  else if (strcmp(name, "Curve")==0)
	return 12;
  else if (strcmp(name, "G")==0)
	return 13;
  else if (strcmp(name, "KinEn0")==0)
	return 14;
  else if (strcmp(name, "Length")==0)
	return 15;
  else if (strcmp(name, "Mass0")==0)
	return 16;
  else if (strcmp(name, "clight")==0)
	return 17;
  else if (strcmp(name, "m0")==0)
	return 18;
  else if (strcmp(name, "offset")==0)
	return 19;
  else if (strcmp(name, "q")==0)
	return 20;
  else {
	fprintf(stderr, "Invalid name %s for getindex call\n", name);
	return 0.0/0.0;
	}
}


int heav(double x_, double *p_, double *wk_, double *xv_) {
  if (x_>0.0) {return 1;} else {return 0;}
}

void auxvars(unsigned n_, unsigned np_, double t, double *Y_, double *p_, double *f_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_){


}

double NaN_event(unsigned n_, double t, double *Y_, double *p_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_) {
return  pow(Pc(dK, p_, wk_, xv_),2)-pow(Pc(0, p_, wk_, xv_),2)*(pow(px,2)+pow(py,2))-offset; 
}

double passto19(unsigned n_, double t, double *Y_, double *p_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_) {
return  s-0.05; 
}


void massMatrix(unsigned n_, unsigned np_, double t, double *Y_, double *p_, double **f_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_) {
}

void jacobian(unsigned n_, unsigned np_, double t, double *Y_, double *p_, double **f_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_) {
}

void jacobianParam(unsigned n_, unsigned np_, double t, double *Y_, double *p_, double **f_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_) {
}
