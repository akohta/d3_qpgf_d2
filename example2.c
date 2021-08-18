// example of second derivative 
#include "d3_qpgf_d2.h"

// verification function using central difference
void central_diff(double *r,QPD2 *qd);

int main()
{
  QPD2 qd;
  double complex qG,dqG[3],d2qG[6];
  double lambda,eps,theta,phi,r[3];
  int err;

  eps=1.0e-15; // requested relative error
  lambda=1.0;  // wave length
  theta=0.3;   // parameter of wave number vector 
  phi=0.5;     // kx=k*sin(theta)*cos(phi), ky=k*sin(theta)*sin(phi), kz=k*cos(theta)

  qd.k=2.0*M_PI/lambda;              // wave number 
  qd.vk[0]=qd.k*sin(theta)*cos(phi); // x component of wave number vector
  qd.vk[1]=qd.k*sin(theta)*sin(phi); // y component of wave number vector
  qd.vk[2]=qd.k*cos(theta);          // z component of wave number vector

  qd.vd1[0]=1.0;  qd.vd1[1]=0.5; // 1st lattice vector, vd1[0]:x component, vd1[1]:y component, z component is zero
  qd.vd2[0]=0.5;  qd.vd2[1]=1.0; // 2nd lattice vector, vd2[0]:x component, vd2[1]:y component, z component is zero

  printf("eps = %2.1e\n",eps);
  printf("|k| = %g\n",qd.k);
  printf("k   = (%8g, %8g, %8g)\n",qd.vk[0],qd.vk[1],qd.vk[2]);
  printf("d1  = (%8g, %8g)\n",qd.vd1[0],qd.vd1[1]);
  printf("d2  = (%8g, %8g)\n\n",qd.vd2[0],qd.vd2[1]);

  r[0]=0.1;
  r[1]=0.4;
  r[2]=10.1;
  printf("r = (%8g, %8g, %8g)\n",r[0],r[1],r[2]);
  printf("----- second derivative of quasi-periodic Green's function -----\n");
  err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,r,eps,&qd);
  printf("d^2qG/dx^2 =% 15.14e %+15.14e I, err=%d\n",creal(d2qG[0]),cimag(d2qG[0]),err);
  printf("d^2qG/dy^2 =% 15.14e %+15.14e I\n",creal(d2qG[1]),cimag(d2qG[1]));
  printf("d^2qG/dz^2 =% 15.14e %+15.14e I\n",creal(d2qG[2]),cimag(d2qG[2]));
  printf("d^2qG/dxdy =% 15.14e %+15.14e I\n",creal(d2qG[3]),cimag(d2qG[3]));
  printf("d^2qG/dydz =% 15.14e %+15.14e I\n",creal(d2qG[4]),cimag(d2qG[4]));
  printf("d^2qG/dzdx =% 15.14e %+15.14e I\n",creal(d2qG[5]),cimag(d2qG[5]));
  // verification 
  central_diff(r,&qd);

  r[0]=0.1;
  r[1]=0.4;
  r[2]=0.07;
  printf("\n\nr = (%8g, %8g, %8g)\n",r[0],r[1],r[2]);
  printf("----- second derivative of quasi-periodic Green's function -----\n");
  err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,r,eps,&qd);
  printf("d^2qG/dx^2 =% 15.14e %+15.14e I, err=%d\n",creal(d2qG[0]),cimag(d2qG[0]),err);
  printf("d^2qG/dy^2 =% 15.14e %+15.14e I\n",creal(d2qG[1]),cimag(d2qG[1]));
  printf("d^2qG/dz^2 =% 15.14e %+15.14e I\n",creal(d2qG[2]),cimag(d2qG[2]));
  printf("d^2qG/dxdy =% 15.14e %+15.14e I\n",creal(d2qG[3]),cimag(d2qG[3]));
  printf("d^2qG/dydz =% 15.14e %+15.14e I\n",creal(d2qG[4]),cimag(d2qG[4]));
  printf("d^2qG/dzdx =% 15.14e %+15.14e I\n",creal(d2qG[5]),cimag(d2qG[5]));
  // verification 
  central_diff(r,&qd);
  
  return 0;
  
  printf("test\n");
  return 0;
}

void central_diff(double *r,QPD2 *qd)
{
  double complex dqGp[3],dqGm[3],fv;
  double h,rt[3],eps;
  
  h=1.0e-7;
  eps=1.0e-15;
  
  int d3hm_qpgf_d2_dqG(double complex *dqG,double *r,double eps,QPD2 *qd);
  printf("----- verification of second derivative using central difference -----\n");
  // x 
  rt[0]=r[0]+h;
  rt[1]=r[1];
  rt[2]=r[2];
  d3hm_qpgf_d2_dqG(dqGp,rt,eps,qd);
  rt[0]=r[0]-h;
  d3hm_qpgf_d2_dqG(dqGm,rt,eps,qd);
  fv=(dqGp[0]-dqGm[0])/(2.0*h);
  printf("d^2qG/dx^2 =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));
  fv=(dqGp[1]-dqGm[1])/(2.0*h);
  printf("d^2qG/dydx =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));
  fv=(dqGp[2]-dqGm[2])/(2.0*h);
  printf("d^2qG/dzdx =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));
  // y 
  rt[0]=r[0];
  rt[1]=r[1]+h;
  rt[2]=r[2];
  d3hm_qpgf_d2_dqG(dqGp,rt,eps,qd);
  rt[1]=r[1]-h;
  d3hm_qpgf_d2_dqG(dqGm,rt,eps,qd);
  fv=(dqGp[0]-dqGm[0])/(2.0*h);
  printf("d^2qG/dxdy =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));
  fv=(dqGp[1]-dqGm[1])/(2.0*h);
  printf("d^2qG/dy^2 =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));
  fv=(dqGp[2]-dqGm[2])/(2.0*h);
  printf("d^2qG/dzdy =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));
  // z
  rt[0]=r[0];
  rt[1]=r[1];
  rt[2]=r[2]+h;
  d3hm_qpgf_d2_dqG(dqGp,rt,eps,qd);
  rt[2]=r[2]-h;
  d3hm_qpgf_d2_dqG(dqGm,rt,eps,qd);
  fv=(dqGp[0]-dqGm[0])/(2.0*h);
  printf("d^2qG/dxdz =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));
  fv=(dqGp[1]-dqGm[1])/(2.0*h);
  printf("d^2qG/dydz =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));
  fv=(dqGp[2]-dqGm[2])/(2.0*h);
  printf("d^2qG/dz^2 =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));  
}
