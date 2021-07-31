#include "d3_qpgf_d2.h"

void normal_sum(double complex *qG,double complex *dqG,double k,double *vk,double *vd1,double *vd2,double *r);

int main()
{
	QPD2 qd;
	double complex qG,dqG[3];
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

	qd.vd1[0]=1.0;	qd.vd1[1]=0.5; // 1st lattice vector, vd1[0]:x component, vd1[1]:y component, z component is zero
	qd.vd2[0]=0.5;	qd.vd2[1]=1.0; // 2nd lattice vector, vd2[0]:x component, vd2[1]:y component, z component is zero

	printf("eps = %2.1e\n",eps);
	printf("|k| = %g\n",qd.k);
	printf("k   = (%8g, %8g, %8g)\n",qd.vk[0],qd.vk[1],qd.vk[2]);
	printf("d1  = (%8g, %8g)\n",qd.vd1[0],qd.vd1[1]);
	printf("d2  = (%8g, %8g)\n\n",qd.vd2[0],qd.vd2[1]);

	r[0]=0.1;
	r[1]=0.4;
	r[2]=10.1;
	printf("r = (%8g, %8g, %8g)\n",r[0],r[1],r[2]);

	printf("----- normal_sum -----\n");
	normal_sum(&qG,dqG,qd.k,qd.vk,qd.vd1,qd.vd2,r);

	// Fourier domain
	printf("----- Fourier domain summation -----\n");
	err=d3hm_qpgf_d2_fd(&qG,dqG,r,eps,&qd);
	printf("qG    =% 15.14e %+15.14e I, err=%d\n",creal(qG),cimag(qG),err);
	printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
	printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
	printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));

	printf("----- Ewald method -----\n");
	err=d3hm_qpgf_d2_ew(&qG,dqG,r,eps,&qd);
	printf("qG    =% 15.14e %+15.14e I, err=%d\n",creal(qG),cimag(qG),err);
	printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
	printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
	printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));
	// auto
	printf("----- Auto -----\n");
	err=d3hm_qpgf_d2(&qG,dqG,r,eps,&qd);
	printf("qG    =% 15.14e %+15.14e I, err=%d\n",creal(qG),cimag(qG),err);
	printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
	printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
	printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));


	r[0]=0.1;
	r[1]=0.4;
	r[2]=0.07;
	printf("\n\nr = (%8g, %8g, %8g)\n",r[0],r[1],r[2]);

	printf("----- normal_sum -----\n");
	normal_sum(&qG,dqG,qd.k,qd.vk,qd.vd1,qd.vd2,r);

	// Fourier domain
	printf("----- Fourier domain summation -----\n");
	err=d3hm_qpgf_d2_fd(&qG,dqG,r,eps,&qd);
	printf("qG    =% 15.14e %+15.14e I, err=%d\n",creal(qG),cimag(qG),err);
	printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
	printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
	printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));
	// Ewald method
	printf("----- Ewald method -----\n");
	err=d3hm_qpgf_d2_ew(&qG,dqG,r,eps,&qd);
	printf("qG    =% 15.14e %+15.14e I, err=%d\n",creal(qG),cimag(qG),err);
	printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
	printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
	printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));
	// auto
	printf("----- Auto -----\n");
	err=d3hm_qpgf_d2(&qG,dqG,r,eps,&qd);
	printf("qG    =% 15.14e %+15.14e I, err=%d\n",creal(qG),cimag(qG),err);
	printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
	printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
	printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));

	return 0;
}

void normal_sum(double complex *qG,double complex *dqG,double k,double *vk,double *vd1,double *vd2,double *r)
{
	double complex tf;
	double ar,arg,Rl[2],tr[3];
	int lmax,l,i;

	lmax=4000;
	printf("sum_{l1=-l}^{l} sum_{l2=-l}^{l} f_{l1,l2}(r)\n");

	// l=0 (l1=l2=0)
	*qG=0.0;
	dqG[0]=0.0;	dqG[1]=0.0;	dqG[2]=0.0;
	ar=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	arg=k*ar;
	tf=(cos(arg)+I*sin(arg))/ar;
	*qG=tf;
	dqG[0]=tf*r[0]*(I*k*ar-1.0)/(ar*ar);
	dqG[1]=tf*r[1]*(I*k*ar-1.0)/(ar*ar);
	dqG[2]=tf*r[2]*(I*k*ar-1.0)/(ar*ar);

	for(l=1;l<=lmax;l++){ // sum_{l1=-l}^{l1=l} sum_{l2=-l}^{l2=l}
		// side of l1=+-l
		for(i=-l;i<=l;i++){
			// l1=+l
			Rl[0]=(double)l*vd1[0]+(double)i*vd2[0];
			Rl[1]=(double)l*vd1[1]+(double)i*vd2[1];
			tr[0]=r[0]+Rl[0];
			tr[1]=r[1]+Rl[1];
			tr[2]=r[2];
			ar=sqrt(tr[0]*tr[0]+tr[1]*tr[1]+tr[2]*tr[2]);
			arg=k*ar+vk[0]*Rl[0]+vk[1]*Rl[1];
			tf=(cos(arg)+I*sin(arg))/ar;
			*qG+=tf;
			dqG[0]+=tf*tr[0]*(I*k*ar-1.0)/(ar*ar);
			dqG[1]+=tf*tr[1]*(I*k*ar-1.0)/(ar*ar);
			dqG[2]+=tf*tr[2]*(I*k*ar-1.0)/(ar*ar);
			// l1=-l
			Rl[0]=-(double)l*vd1[0]+(double)i*vd2[0];
			Rl[1]=-(double)l*vd1[1]+(double)i*vd2[1];
			tr[0]=r[0]+Rl[0];
			tr[1]=r[1]+Rl[1];
			tr[2]=r[2];
			ar=sqrt(tr[0]*tr[0]+tr[1]*tr[1]+tr[2]*tr[2]);
			arg=k*ar+vk[0]*Rl[0]+vk[1]*Rl[1];
			tf=(cos(arg)+I*sin(arg))/ar;
			*qG+=tf;
			dqG[0]+=tf*tr[0]*(I*k*ar-1.0)/(ar*ar);
			dqG[1]+=tf*tr[1]*(I*k*ar-1.0)/(ar*ar);
			dqG[2]+=tf*tr[2]*(I*k*ar-1.0)/(ar*ar);
		}

		// side of l2=+-l, except corner
		for(i=-l+1;i<l;i++){
			// l2=+l
			Rl[0]=(double)i*vd1[0]+(double)l*vd2[0];
			Rl[1]=(double)i*vd1[1]+(double)l*vd2[1];
			tr[0]=r[0]+Rl[0];
			tr[1]=r[1]+Rl[1];
			tr[2]=r[2];
			ar=sqrt(tr[0]*tr[0]+tr[1]*tr[1]+tr[2]*tr[2]);
			arg=k*ar+vk[0]*Rl[0]+vk[1]*Rl[1];
			tf=(cos(arg)+I*sin(arg))/ar;
			*qG+=tf;
			dqG[0]+=tf*tr[0]*(I*k*ar-1.0)/(ar*ar);
			dqG[1]+=tf*tr[1]*(I*k*ar-1.0)/(ar*ar);
			dqG[2]+=tf*tr[2]*(I*k*ar-1.0)/(ar*ar);
			// l2=-l
			Rl[0]=(double)i*vd1[0]-(double)l*vd2[0];
			Rl[1]=(double)i*vd1[1]-(double)l*vd2[1];
			tr[0]=r[0]+Rl[0];
			tr[1]=r[1]+Rl[1];
			tr[2]=r[2];
			ar=sqrt(tr[0]*tr[0]+tr[1]*tr[1]+tr[2]*tr[2]);
			arg=k*ar+vk[0]*Rl[0]+vk[1]*Rl[1];
			tf=(cos(arg)+I*sin(arg))/ar;
			*qG+=tf;
			dqG[0]+=tf*tr[0]*(I*k*ar-1.0)/(ar*ar);
			dqG[1]+=tf*tr[1]*(I*k*ar-1.0)/(ar*ar);
			dqG[2]+=tf*tr[2]*(I*k*ar-1.0)/(ar*ar);
		}

		 if(l%(lmax/2)==0){
			 printf("l=%5d\n",l);
			 printf("qG    =% 15.14e %+15.14e I\n",creal(*qG)/(4.0*M_PI),cimag(*qG)/(4.0*M_PI));
			 printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0])/(4.0*M_PI),cimag(dqG[0])/(4.0*M_PI));
			 printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1])/(4.0*M_PI),cimag(dqG[1])/(4.0*M_PI));
			 printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2])/(4.0*M_PI),cimag(dqG[2])/(4.0*M_PI));
		 }

	}

	*qG/=4.0*M_PI;
	dqG[0]/=4.0*M_PI;
	dqG[1]/=4.0*M_PI;
	dqG[2]/=4.0*M_PI;
}

