/*
 * d3_qpgf_d2.c
 *
 *  Created on: Jun 11, 2019
 *      Author: ohta
 */

#include "d3_qpgf_d2.h"

int d3hm_qpgf_d2_fd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd)
{
	int d3hm_qpgf_d2_fd_sd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd);

	double complex cf;
	double rt[3],i_detd,arg;
	int l1,l2,err;

	i_detd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
	l1=(int)floor( (r[0]*qd->vd2[1]-r[1]*qd->vd2[0])*i_detd+0.5 );
	l2=(int)floor(-(r[0]*qd->vd1[1]-r[1]*qd->vd1[0])*i_detd+0.5 );

	rt[0]=r[0]-(double)l1*qd->vd1[0]-(double)l2*qd->vd2[0];
	rt[1]=r[1]-(double)l1*qd->vd1[1]-(double)l2*qd->vd2[1];
	rt[2]=r[2];

	err=d3hm_qpgf_d2_fd_sd(qG,dqG,rt,eps,qd);

	if(l1!=0 || l2!=0){
		arg=-(qd->vk[0]*((double)l1*qd->vd1[0]+(double)l2*qd->vd2[0])+qd->vk[1]*((double)l1*qd->vd1[1]+(double)l2*qd->vd2[1]));
		cf=cos(arg)+I*sin(arg);

		*qG*=cf;
		dqG[0]*=cf;
		dqG[1]*=cf;
		dqG[2]*=cf;
	}

	return err;
}

int d3hm_qpgf_d2_ew(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd)
{
	int d3hm_qpgf_d2_ew_sd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd);

	double complex cf;
	double rt[3],i_detd,arg;
	int l1,l2,err;

	i_detd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
	l1=(int)floor( (r[0]*qd->vd2[1]-r[1]*qd->vd2[0])*i_detd+0.5 );
	l2=(int)floor(-(r[0]*qd->vd1[1]-r[1]*qd->vd1[0])*i_detd+0.5 );

	rt[0]=r[0]-(double)l1*qd->vd1[0]-(double)l2*qd->vd2[0];
	rt[1]=r[1]-(double)l1*qd->vd1[1]-(double)l2*qd->vd2[1];
	rt[2]=r[2];

	err=d3hm_qpgf_d2_ew_sd(qG,dqG,rt,eps,qd);

	if(l1!=0 || l2!=0){
		arg=-(qd->vk[0]*((double)l1*qd->vd1[0]+(double)l2*qd->vd2[0])+qd->vk[1]*((double)l1*qd->vd1[1]+(double)l2*qd->vd2[1]));
		cf=cos(arg)+I*sin(arg);

		*qG*=cf;
		dqG[0]*=cf;
		dqG[1]*=cf;
		dqG[2]*=cf;
	}

	return err;
}

int d3hm_qpgf_d2(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd)
{
	double bN,aNx,aNy,k2,idd,aNpx,aNmx,aNpy,aNmy,ad12,ad22,cpd,bNp,bNm;
	int n1,n2,err;

	k2=qd->k*qd->k;
	idd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
	cpd=2.0*M_PI*idd;
	ad12=qd->vd1[0]*qd->vd1[0]+qd->vd1[1]*qd->vd1[1];
	ad22=qd->vd2[0]*qd->vd2[0]+qd->vd2[1]*qd->vd2[1];

	n1=(int)round( (qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI) );
	n2=(int)round( (qd->vk[0]*qd->vd2[0]+qd->vk[1]*qd->vd2[1])/(2.0*M_PI) );

	aNx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
	aNy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
	bN=sqrt(k2-aNx*aNx-aNy*aNy);

	if(ad22<ad12){ // |g_1| < |g_2|
		aNpx=-qd->vk[0]+(double)(n1+S_LIMIT)*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
		aNpy=-qd->vk[1]-(double)(n1+S_LIMIT)*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
		aNmx=-qd->vk[0]+(double)(n1-S_LIMIT)*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
		aNmy=-qd->vk[1]-(double)(n1-S_LIMIT)*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
	}
	else { // |g_1| >= |g_2|
		aNpx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)(n2+S_LIMIT)*cpd*qd->vd1[1];
		aNpy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)(n2+S_LIMIT)*cpd*qd->vd1[0];
		aNmx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)(n2-S_LIMIT)*cpd*qd->vd1[1];
		aNmy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)(n2-S_LIMIT)*cpd*qd->vd1[0];
	}
	bNp=k2-aNpx*aNpx-aNpy*aNpy;
	bNm=k2-aNmx*aNmx-aNmy*aNmy;

	if(bNp>0.0 || bNm>0.0){ // Ewald method
		err=d3hm_qpgf_d2_ew(qG,dqG,r,eps,qd);
		if(err<0) err-=20;
		return err;
	}
	else {
		bNp=sqrt(-bNp);
		bNm=sqrt(-bNm);
		if( bN/bNp*exp(-bNp*fabs(r[2]))<eps && bN/bNm*exp(-bNm*fabs(r[2]))<eps ){ // Fourier domain sum
			err=d3hm_qpgf_d2_fd(qG,dqG,r,eps,qd);
			if(err<0) err-=10;
			return err;
		}
		else { // Ewald method
			err=d3hm_qpgf_d2_ew(qG,dqG,r,eps,qd);
			if(err<0) err-=20;
			return err;
		}
	}
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////
int d3hm_qpgf_d2_fd_sd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd)
{
	//function prototype
	int fd_sum_n2(double complex *tqG,double complex *tdqG,int n1,double *vg1,double *vg2,double *r,double eps,QPD2 *qd,int *err);

	double complex oqG,odqG[3];
	double vg1[2],vg2[2],detd,idd;
	int n1,err,cc,nc;

	if(r[2]==0.0) return -4; // z==0.0

	detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
	if(detd==0.0) return -3; // determinant d == 0
	idd=1.0/detd;

	vg1[0]= 2.0*M_PI*idd*qd->vd2[1];
	vg1[1]=-2.0*M_PI*idd*qd->vd2[0];
	vg2[0]=-2.0*M_PI*idd*qd->vd1[1];
	vg2[1]= 2.0*M_PI*idd*qd->vd1[0];

	nc=(int)round((qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI));

	cc=0;
	*qG=0.0;
	dqG[0]=0.0;
	dqG[1]=0.0;
	dqG[2]=0.0;

	n1=nc;
	cc+=fd_sum_n2(qG,dqG,n1,vg1,vg2,r,eps,qd,&err);	if(err<0) return err;

	// n1 > nc
	for(n1=nc+1;n1<nc+FD_LIMIT;n1++){
		oqG=*qG;
		odqG[0]=dqG[0];
		odqG[1]=dqG[1];
		odqG[2]=dqG[2];
		cc+=fd_sum_n2(qG,dqG,n1,vg1,vg2,r,eps,qd,&err);		if(err<0) return err;

		if(	(  2.0*fabs(creal(*qG)-creal(oqG))<eps*(fabs(creal(*qG))+fabs(creal(oqG)))
				&& 2.0*fabs(cimag(*qG)-cimag(oqG))<eps*(fabs(cimag(*qG))+fabs(cimag(oqG))) )
			&& ( 2.0*fabs(creal(dqG[0])-creal(odqG[0]))<eps*(fabs(creal(dqG[0]))+fabs(creal(odqG[0])))
				&& 2.0*fabs(cimag(dqG[0])-cimag(odqG[0]))<eps*(fabs(cimag(dqG[0]))+fabs(cimag(odqG[0]))) )
			&& ( 2.0*fabs(creal(dqG[1])-creal(odqG[1]))<eps*(fabs(creal(dqG[1]))+fabs(creal(odqG[1])))
				&& 2.0*fabs(cimag(dqG[1])-cimag(odqG[1]))<eps*(fabs(cimag(dqG[1]))+fabs(cimag(odqG[1]))) )
			&& ( 2.0*fabs(creal(dqG[2])-creal(odqG[2]))<eps*(fabs(creal(dqG[2]))+fabs(creal(odqG[2])))
				&& 2.0*fabs(cimag(dqG[2])-cimag(odqG[2]))<eps*(fabs(cimag(dqG[2]))+fabs(cimag(odqG[2]))) )
		) break;
	}
	if(n1==nc+FD_LIMIT) return -1; // summation limit

	// n1 < nc
	for(n1=nc-1;n1>nc-FD_LIMIT;n1--){
		oqG=*qG;
		odqG[0]=dqG[0];
		odqG[1]=dqG[1];
		odqG[2]=dqG[2];
		cc+=fd_sum_n2(qG,dqG,n1,vg1,vg2,r,eps,qd,&err);		if(err<0) return err;

		if(	(  2.0*fabs(creal(*qG)-creal(oqG))<eps*(fabs(creal(*qG))+fabs(creal(oqG)))
				&& 2.0*fabs(cimag(*qG)-cimag(oqG))<eps*(fabs(cimag(*qG))+fabs(cimag(oqG))) )
			&& ( 2.0*fabs(creal(dqG[0])-creal(odqG[0]))<eps*(fabs(creal(dqG[0]))+fabs(creal(odqG[0])))
				&& 2.0*fabs(cimag(dqG[0])-cimag(odqG[0]))<eps*(fabs(cimag(dqG[0]))+fabs(cimag(odqG[0]))) )
			&& ( 2.0*fabs(creal(dqG[1])-creal(odqG[1]))<eps*(fabs(creal(dqG[1]))+fabs(creal(odqG[1])))
				&& 2.0*fabs(cimag(dqG[1])-cimag(odqG[1]))<eps*(fabs(cimag(dqG[1]))+fabs(cimag(odqG[1]))) )
			&& ( 2.0*fabs(creal(dqG[2])-creal(odqG[2]))<eps*(fabs(creal(dqG[2]))+fabs(creal(odqG[2])))
				&& 2.0*fabs(cimag(dqG[2])-cimag(odqG[2]))<eps*(fabs(cimag(dqG[2]))+fabs(cimag(odqG[2]))) )
		) break;
	}
	if(n1==nc-FD_LIMIT) return -1; // summation limit

	*qG*=0.5*I*fabs(idd);
	dqG[0]*=-0.5*fabs(idd);
	dqG[1]*=-0.5*fabs(idd);
	dqG[2]*=-0.5*fabs(idd)*r[2]/fabs(r[2]);

	return cc;
}

int fd_sum_n2(double complex *tqG,double complex *tdqG,int n1,double *vg1,double *vg2,double *r,double eps,QPD2 *qd,int *err)
{
	// function prototype
	int fd_cFn(double complex *Fn,double complex *dFn,int n1,int n2,double *vg1,double *vg2,double *r,QPD2 *qd);

	double complex Fn,dFn[3];
	int n2,cc,nc;

	nc=(int)round( (qd->vk[0]*vg2[0]+qd->vk[1]*vg2[1]-(double)n1*(vg1[0]*vg2[0]+vg1[1]*vg2[1]))/(vg2[0]*vg2[0]+vg2[1]*vg2[1]));

	cc=0;
	n2=nc;
	*err=fd_cFn(&Fn,dFn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
	*tqG+=Fn;
	tdqG[0]+=dFn[0];
	tdqG[1]+=dFn[1];
	tdqG[2]+=dFn[2];
	cc+=1;

	// n2 > nc
	for(n2=nc+1;n2<nc+FD_LIMIT;n2++){
		*err=fd_cFn(&Fn,dFn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
		*tqG+=Fn;
		tdqG[0]+=dFn[0];
		tdqG[1]+=dFn[1];
		tdqG[2]+=dFn[2];
		cc+=1;

		if(*err>0){
			if( (fabs(creal(Fn))<fabs(creal(*tqG))*eps && fabs(cimag(Fn))<fabs(cimag(*tqG))*eps)
					&& (fabs(creal(dFn[0]))<fabs(creal(tdqG[0]))*eps && fabs(cimag(dFn[0]))<fabs(cimag(tdqG[0]))*eps)
					&& (fabs(creal(dFn[1]))<fabs(creal(tdqG[1]))*eps && fabs(cimag(dFn[1]))<fabs(cimag(tdqG[1]))*eps)
					&& (fabs(creal(dFn[2]))<fabs(creal(tdqG[2]))*eps && fabs(cimag(dFn[2]))<fabs(cimag(tdqG[2]))*eps)	)break;
		}
	}
	if(n2==nc+FD_LIMIT){
		*err=-1; // summation limit
		return 0;
	}

	// n2 < nc
	for(n2=nc-1;n2>nc-FD_LIMIT;n2--){
		*err=fd_cFn(&Fn,dFn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
		*tqG+=Fn;
		tdqG[0]+=dFn[0];
		tdqG[1]+=dFn[1];
		tdqG[2]+=dFn[2];
		cc+=1;

		if(*err>0){
			if( (fabs(creal(Fn))<fabs(creal(*tqG))*eps && fabs(cimag(Fn))<fabs(cimag(*tqG))*eps)
					&& (fabs(creal(dFn[0]))<fabs(creal(tdqG[0]))*eps && fabs(cimag(dFn[0]))<fabs(cimag(tdqG[0]))*eps)
					&& (fabs(creal(dFn[1]))<fabs(creal(tdqG[1]))*eps && fabs(cimag(dFn[1]))<fabs(cimag(tdqG[1]))*eps)
					&& (fabs(creal(dFn[2]))<fabs(creal(tdqG[2]))*eps && fabs(cimag(dFn[2]))<fabs(cimag(tdqG[2]))*eps)	)break;
		}
	}
	if(n2==nc-FD_LIMIT){
		*err=-1; // summation limit
		return 0;
	}

	return cc;
}

int fd_cFn(double complex *Fn,double complex *dFn,int n1,int n2,double *vg1,double *vg2,double *r,QPD2 *qd)
{
	double an[2],bn2,bn,arg;

	an[0]=-qd->vk[0]+(double)n1*vg1[0]+(double)n2*vg2[0];
	an[1]=-qd->vk[1]+(double)n1*vg1[1]+(double)n2*vg2[1];
	bn2=qd->k*qd->k-an[0]*an[0]-an[1]*an[1];

	if(bn2>0.0){
		bn=sqrt(bn2);
		arg=an[0]*r[0]+an[1]*r[1]+bn*fabs(r[2]);
		dFn[2]=cos(arg)+I*sin(arg);
		*Fn=dFn[2]/bn;
		dFn[0]=*Fn*an[0];
		dFn[1]=*Fn*an[1];
		return 0;
	}
	else if(bn2<0.0){
		bn=sqrt(-bn2);
		arg=an[0]*r[0]+an[1]*r[1];
		dFn[2]=(cos(arg)+I*sin(arg))*exp(-bn*fabs(r[2]));
		*Fn=-I*dFn[2]/bn;
		dFn[0]=*Fn*an[0];
		dFn[1]=*Fn*an[1];
		return 1;
	}
	else {
		return -2; // bn==0.0
	}
}


int d3hm_qpgf_d2_ew_sd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd)
{
	// function prototype
	int ew_qG1(double complex *qG1,double complex *dqG1,double *r,double eps,double veps,QPD2 *qd);
	int ew_qG2(double complex *qG2,double complex *dqG2,double *r,double eps,double veps,QPD2 *qd);

	double complex qG1,dqG1[3],qG2,dqG2[3];
	double detd,veps,adp,adm,r2_min;
	int c1,c2,l1,l2;

	detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
	if(detd==0.0) return -4; // determinant d == 0

	adp=sqrt(pow(qd->vd1[0]+qd->vd2[0],2)+pow(qd->vd1[1]+qd->vd2[1],2));
	adm=sqrt(pow(qd->vd1[0]-qd->vd2[0],2)+pow(qd->vd1[1]-qd->vd2[1],2));
	veps=sqrt(fabs(detd)/(4.0*M_PI)*adp/adm);

	//l1=(int)round(-idd*(r[0]*qd->vd2[1]-r[1]*qd->vd2[0]));
	l1=0;
	//l2=(int)round( idd*(r[0]*qd->vd1[1]-r[1]*qd->vd1[0]));
	l2=0;
	r2_min=pow(r[0]+(double)l1*qd->vd1[0]+(double)l2*qd->vd2[0],2)
			+  pow(r[1]+(double)l1*qd->vd1[1]+(double)l2*qd->vd2[1],2)+r[2]*r[2];

	if( qd->k*qd->k*veps*veps-r2_min/(4.0*veps*veps)>EW_EXPM )
		veps=sqrt(EW_EXPM+sqrt(qd->k*qd->k*r2_min+EW_EXPM*EW_EXPM))/(sqrt(2.0)*qd->k);

	// qG1
	c1=ew_qG1(&qG1,dqG1,r,eps,veps,qd);
	*qG=qG1;
	dqG[0]=dqG1[0];
	dqG[1]=dqG1[1];
	dqG[2]=dqG1[2];
	if(c1<0) return c1;
/*
	printf("c1=%d\n",c1);
	printf("qG1    =% 15.14e %+15.14e\n",creal(qG1),cimag(qG1));
	printf("dqG1[0]=% 15.14e %+15.14e\n",creal(dqG1[0]),cimag(dqG1[0]));
	printf("dqG1[1]=% 15.14e %+15.14e\n",creal(dqG1[1]),cimag(dqG1[1]));
	printf("dqG1[2]=% 15.14e %+15.14e\n",creal(dqG1[2]),cimag(dqG1[2]));
*/
	// qG2
	c2=ew_qG2(&qG2,dqG2,r,eps,veps,qd);
	*qG+=qG2;
	dqG[0]+=dqG2[0];
	dqG[1]+=dqG2[1];
	dqG[2]+=dqG2[2];
	if(c2<0) return c2;
/*
	printf("c2=%d\n",c2);
	printf("qG2    =% 15.14e %+15.14e\n",creal(qG2),cimag(qG2));
	printf("dqG2[0]=% 15.14e %+15.14e\n",creal(dqG2[0]),cimag(dqG2[0]));
	printf("dqG2[1]=% 15.14e %+15.14e\n",creal(dqG2[1]),cimag(dqG2[1]));
	printf("dqG2[2]=% 15.14e %+15.14e\n",creal(dqG2[2]),cimag(dqG2[2]));
*/
	return c1+c2;
}

int ew_qG1(double complex *qG1,double complex *dqG1,double *r,double eps,double veps,QPD2 *qd)
{
	// function prototype
	int ew_qG1_sum_l2(double complex *tqG,double complex *tdqG,int l1,double *r,double eps,double veps,QPD2 *qd);

	double complex oqG,odqG[3];
	int l1,cc,lc;

	//lc=(int)round(-(r[0]*qd->vd2[1]-r[1]*qd->vd2[0])/(qd->vd1[0]*qd->vd2[1]-qd->vd1[1]*qd->vd2[0]));
	lc=0;

	cc=0;
	*qG1=0.0;
	dqG1[0]=0.0;
	dqG1[1]=0.0;
	dqG1[2]=0.0;

	l1=lc;
	cc+=ew_qG1_sum_l2(qG1,dqG1,l1,r,eps,veps,qd);	if(cc<0) return cc;

	// l1 > lc
	for(l1=lc+1;l1<lc+EW_LIMIT;l1++){
		oqG=*qG1;
		odqG[0]=dqG1[0];
		odqG[1]=dqG1[1];
		odqG[2]=dqG1[2];
		cc+=ew_qG1_sum_l2(qG1,dqG1,l1,r,eps,veps,qd); if(cc<0) return cc;

		if(	(  2.0*fabs(creal(*qG1)-creal(oqG))<eps*(fabs(creal(*qG1))+fabs(creal(oqG)))
				&& 2.0*fabs(cimag(*qG1)-cimag(oqG))<eps*(fabs(cimag(*qG1))+fabs(cimag(oqG))) )
					&& ( 2.0*fabs(creal(dqG1[0])-creal(odqG[0]))<eps*(fabs(creal(dqG1[0]))+fabs(creal(odqG[0])))
						&& 2.0*fabs(cimag(dqG1[0])-cimag(odqG[0]))<eps*(fabs(cimag(dqG1[0]))+fabs(cimag(odqG[0]))) )
					&& ( 2.0*fabs(creal(dqG1[1])-creal(odqG[1]))<eps*(fabs(creal(dqG1[1]))+fabs(creal(odqG[1])))
						&& 2.0*fabs(cimag(dqG1[1])-cimag(odqG[1]))<eps*(fabs(cimag(dqG1[1]))+fabs(cimag(odqG[1]))) )
					&& ( 2.0*fabs(creal(dqG1[2])-creal(odqG[2]))<eps*(fabs(creal(dqG1[2]))+fabs(creal(odqG[2])))
						&& 2.0*fabs(cimag(dqG1[2])-cimag(odqG[2]))<eps*(fabs(cimag(dqG1[2]))+fabs(cimag(odqG[2]))) )
		) break;
	}
	if(l1==lc+EW_LIMIT){
		return -1;
	}

	// l1 < lc
	for(l1=lc-1;l1>lc-EW_LIMIT;l1--){
		oqG=*qG1;
		odqG[0]=dqG1[0];
		odqG[1]=dqG1[1];
		odqG[2]=dqG1[2];
		cc+=ew_qG1_sum_l2(qG1,dqG1,l1,r,eps,veps,qd); if(cc<0) return cc;

		if(	(  2.0*fabs(creal(*qG1)-creal(oqG))<eps*(fabs(creal(*qG1))+fabs(creal(oqG)))
				&& 2.0*fabs(cimag(*qG1)-cimag(oqG))<eps*(fabs(cimag(*qG1))+fabs(cimag(oqG))) )
					&& ( 2.0*fabs(creal(dqG1[0])-creal(odqG[0]))<eps*(fabs(creal(dqG1[0]))+fabs(creal(odqG[0])))
						&& 2.0*fabs(cimag(dqG1[0])-cimag(odqG[0]))<eps*(fabs(cimag(dqG1[0]))+fabs(cimag(odqG[0]))) )
					&& ( 2.0*fabs(creal(dqG1[1])-creal(odqG[1]))<eps*(fabs(creal(dqG1[1]))+fabs(creal(odqG[1])))
						&& 2.0*fabs(cimag(dqG1[1])-cimag(odqG[1]))<eps*(fabs(cimag(dqG1[1]))+fabs(cimag(odqG[1]))) )
					&& ( 2.0*fabs(creal(dqG1[2])-creal(odqG[2]))<eps*(fabs(creal(dqG1[2]))+fabs(creal(odqG[2])))
						&& 2.0*fabs(cimag(dqG1[2])-cimag(odqG[2]))<eps*(fabs(cimag(dqG1[2]))+fabs(cimag(odqG[2]))) )
		) break;
	}
	if(l1==lc-EW_LIMIT){
		return -1;
	}

	*qG1/=4.0*M_PI;
	dqG1[0]/=-4.0*M_PI;
	dqG1[1]/=-4.0*M_PI;
	dqG1[2]*=-r[2]/(4.0*M_PI);

	return cc;
}

int ew_qG1_sum_l2(double complex *tqG,double complex *tdqG,int l1,double *r,double eps,double veps,QPD2 *qd)
{
	// function prototype
	void ew_qG1_cFl(double complex *Fl,double complex *dFl,int l1,int l2,double *r,double eps,double veps,QPD2 *qd);

	double complex Fl,dFl[3];
	int cc,l2,lc;

	lc=-(int)round( ((r[0]+(double)l1*qd->vd1[0])*qd->vd2[0]+(r[1]+(double)l1*qd->vd1[1])*qd->vd2[1])
										/(qd->vd2[0]*qd->vd2[0]+qd->vd2[1]*qd->vd2[1]) );

	cc=0;
	l2=lc;
	cc+=1;
	ew_qG1_cFl(&Fl,dFl,l1,l2,r,eps,veps,qd);
	*tqG+=Fl;
	tdqG[0]+=dFl[0];
	tdqG[1]+=dFl[1];
	tdqG[2]+=dFl[2];

	// l2 > lc
	for(l2=lc+1;l2<lc+EW_LIMIT;l2++){
		cc+=1;
		ew_qG1_cFl(&Fl,dFl,l1,l2,r,eps,veps,qd);
		*tqG+=Fl;
		tdqG[0]+=dFl[0];
		tdqG[1]+=dFl[1];
		tdqG[2]+=dFl[2];

		if( (fabs(creal(Fl))<fabs(creal(*tqG))*eps && fabs(cimag(Fl))<fabs(cimag(*tqG))*eps)
				&& (fabs(creal(dFl[0]))<fabs(creal(tdqG[0]))*eps && fabs(cimag(dFl[0]))<fabs(cimag(tdqG[0]))*eps)
				&& (fabs(creal(dFl[1]))<fabs(creal(tdqG[1]))*eps && fabs(cimag(dFl[1]))<fabs(cimag(tdqG[1]))*eps)
				&& (fabs(creal(dFl[2]))<fabs(creal(tdqG[2]))*eps && fabs(cimag(dFl[2]))<fabs(cimag(tdqG[2]))*eps)	)break;
	}
	if(l2==lc+EW_LIMIT){
		return -1;
	}

	// l2 < 0
	for(l2=lc-1;l2>lc-EW_LIMIT;l2--){
		cc+=1;
		ew_qG1_cFl(&Fl,dFl,l1,l2,r,eps,veps,qd);
		*tqG+=Fl;
		tdqG[0]+=dFl[0];
		tdqG[1]+=dFl[1];
		tdqG[2]+=dFl[2];

		if( (fabs(creal(Fl))<fabs(creal(*tqG))*eps && fabs(cimag(Fl))<fabs(cimag(*tqG))*eps)
				&& (fabs(creal(dFl[0]))<fabs(creal(tdqG[0]))*eps && fabs(cimag(dFl[0]))<fabs(cimag(tdqG[0]))*eps)
				&& (fabs(creal(dFl[1]))<fabs(creal(tdqG[1]))*eps && fabs(cimag(dFl[1]))<fabs(cimag(tdqG[1]))*eps)
				&& (fabs(creal(dFl[2]))<fabs(creal(tdqG[2]))*eps && fabs(cimag(dFl[2]))<fabs(cimag(tdqG[2]))*eps)	)break;
	}
	if(l2==lc-EW_LIMIT){
		return -1;
	}

	return cc;
}

void ew_qG1_cFl(double complex *Fl,double complex *dFl,int l1,int l2,double *r,double eps,double veps,QPD2 *qd)
{
	double complex w,ce;
	double rl,i_rl,Rl[2],ke,i_ve,rl_2e,arg,ef;

	Rl[0]=(double)l1*qd->vd1[0]+(double)l2*qd->vd2[0];
	Rl[1]=(double)l1*qd->vd1[1]+(double)l2*qd->vd2[1];
	rl=sqrt((r[0]+Rl[0])*(r[0]+Rl[0])+(r[1]+Rl[1])*(r[1]+Rl[1])+r[2]*r[2]);
	i_rl=1.0/rl;

	ke=qd->k*veps;
	i_ve=1.0/veps;
	rl_2e=rl*0.5*i_ve;

	arg=qd->vk[0]*Rl[0]+qd->vk[1]*Rl[1];
	ce=cos(arg)+I*sin(arg);
	ef=exp(ke*ke-rl_2e*rl_2e);
	w=Faddeeva_w(ke+I*rl_2e,eps);

	*Fl=ce*i_rl*ef*creal(w);
	dFl[0]=(r[0]+Rl[0])*ce*i_rl*i_rl*ef*(i_rl*creal(w)-qd->k*cimag(w)+i_ve/sqrt(M_PI));
	dFl[1]=(r[1]+Rl[1])*ce*i_rl*i_rl*ef*(i_rl*creal(w)-qd->k*cimag(w)+i_ve/sqrt(M_PI));
	dFl[2]=ce*i_rl*i_rl*ef*(i_rl*creal(w)-qd->k*cimag(w)+i_ve/sqrt(M_PI));
}

int ew_qG2(double complex *qG2,double complex *dqG2,double *r,double eps,double veps,QPD2 *qd)
{
	// function prototype
	int ew_qG2_sum_n2(double complex *tqG,double complex *tdqG,int n1,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd,int *err);

	double complex oqG,odqG[3];
	double vg1[2],vg2[2],detd,idd;
	int n1,cc,err,nc;

	detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
	if(detd==0.0) return -4; // determinant d == 0
	idd=1.0/detd;

	nc=(int)round( (qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI) );

	vg1[0]= 2.0*M_PI*idd*qd->vd2[1];
	vg1[1]=-2.0*M_PI*idd*qd->vd2[0];
	vg2[0]=-2.0*M_PI*idd*qd->vd1[1];
	vg2[1]= 2.0*M_PI*idd*qd->vd1[0];

	cc=0;
	*qG2=0.0;
	dqG2[0]=0.0;
	dqG2[1]=0.0;
	dqG2[2]=0.0;

	n1=nc;
	cc+=ew_qG2_sum_n2(qG2,dqG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

	// n1 > nc
	for(n1=nc+1;n1<nc+EW_LIMIT;n1++){
		oqG=*qG2;
		odqG[0]=dqG2[0];
		odqG[1]=dqG2[1];
		odqG[2]=dqG2[2];
		cc+=ew_qG2_sum_n2(qG2,dqG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

		if(	(  2.0*fabs(creal(*qG2)-creal(oqG))<eps*(fabs(creal(*qG2))+fabs(creal(oqG)))
				&& 2.0*fabs(cimag(*qG2)-cimag(oqG))<eps*(fabs(cimag(*qG2))+fabs(cimag(oqG))) )
			&& ( 2.0*fabs(creal(dqG2[0])-creal(odqG[0]))<eps*(fabs(creal(dqG2[0]))+fabs(creal(odqG[0])))
				&& 2.0*fabs(cimag(dqG2[0])-cimag(odqG[0]))<eps*(fabs(cimag(dqG2[0]))+fabs(cimag(odqG[0]))) )
			&& ( 2.0*fabs(creal(dqG2[1])-creal(odqG[1]))<eps*(fabs(creal(dqG2[1]))+fabs(creal(odqG[1])))
				&& 2.0*fabs(cimag(dqG2[1])-cimag(odqG[1]))<eps*(fabs(cimag(dqG2[1]))+fabs(cimag(odqG[1]))) )
			&& ( 2.0*fabs(creal(dqG2[2])-creal(odqG[2]))<eps*(fabs(creal(dqG2[2]))+fabs(creal(odqG[2])))
				&& 2.0*fabs(cimag(dqG2[2])-cimag(odqG[2]))<eps*(fabs(cimag(dqG2[2]))+fabs(cimag(odqG[2]))) )
		) break;
	}
	if(n1==nc+EW_LIMIT){
		return -2;
	}

	// n1 < nc
	for(n1=nc-1;n1>nc-EW_LIMIT;n1--){
		oqG=*qG2;
		odqG[0]=dqG2[0];
		odqG[1]=dqG2[1];
		odqG[2]=dqG2[2];
		cc+=ew_qG2_sum_n2(qG2,dqG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

		if(	(  2.0*fabs(creal(*qG2)-creal(oqG))<eps*(fabs(creal(*qG2))+fabs(creal(oqG)))
				&& 2.0*fabs(cimag(*qG2)-cimag(oqG))<eps*(fabs(cimag(*qG2))+fabs(cimag(oqG))) )
			&& ( 2.0*fabs(creal(dqG2[0])-creal(odqG[0]))<eps*(fabs(creal(dqG2[0]))+fabs(creal(odqG[0])))
				&& 2.0*fabs(cimag(dqG2[0])-cimag(odqG[0]))<eps*(fabs(cimag(dqG2[0]))+fabs(cimag(odqG[0]))) )
			&& ( 2.0*fabs(creal(dqG2[1])-creal(odqG[1]))<eps*(fabs(creal(dqG2[1]))+fabs(creal(odqG[1])))
				&& 2.0*fabs(cimag(dqG2[1])-cimag(odqG[1]))<eps*(fabs(cimag(dqG2[1]))+fabs(cimag(odqG[1]))) )
			&& ( 2.0*fabs(creal(dqG2[2])-creal(odqG[2]))<eps*(fabs(creal(dqG2[2]))+fabs(creal(odqG[2])))
				&& 2.0*fabs(cimag(dqG2[2])-cimag(odqG[2]))<eps*(fabs(cimag(dqG2[2]))+fabs(cimag(odqG[2]))) )
		) break;
	}
	if(n1==nc-EW_LIMIT){
		return -2;
	}

	*qG2*=0.25*I*fabs(idd);
	dqG2[0]*=-0.25*fabs(idd);
	dqG2[1]*=-0.25*fabs(idd);
	dqG2[2]*= 0.25*fabs(idd);

	return cc;
}

int ew_qG2_sum_n2(double complex *tqG,double complex *tdqG,int n1,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd,int *err)
{
	// function prototype
	int ew_qG2_cFn(double complex *Fn,double complex *dFn,int n1,int n2,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd);

	double complex Fn,dFn[3];
	int cc,n2,nc;

	nc=(int)round(-((-qd->vk[0]+(double)n1*vg1[0])*vg2[0]+(-qd->vk[1]+(double)n1*vg1[1])*vg2[1])/(vg2[0]*vg2[0]+vg2[1]*vg2[1]));
	//printf("nc=%d\n",nc);

	cc=0;
	n2=nc;
	cc+=1;
	*err=ew_qG2_cFn(&Fn,dFn,n1,n2,r,vg1,vg2,eps,veps,qd);	if(*err<0) return 0;
	*tqG+=Fn;
	tdqG[0]+=dFn[0];
	tdqG[1]+=dFn[1];
	tdqG[2]+=dFn[2];

	// n2 > nc
	for(n2=nc+1;n2<nc+EW_LIMIT;n2++){
		cc+=1;
		*err=ew_qG2_cFn(&Fn,dFn,n1,n2,r,vg1,vg2,eps,veps,qd); if(*err<0) return 0;
		*tqG+=Fn;
		tdqG[0]+=dFn[0];
		tdqG[1]+=dFn[1];
		tdqG[2]+=dFn[2];

		if(*err>0){
			if(  (fabs(creal(Fn))<fabs(creal(*tqG))*eps && fabs(cimag(Fn))<fabs(cimag(*tqG))*eps)
				&& (fabs(creal(dFn[0]))<fabs(creal(tdqG[0]))*eps && fabs(cimag(dFn[0]))<fabs(cimag(tdqG[0]))*eps)
				&& (fabs(creal(dFn[1]))<fabs(creal(tdqG[1]))*eps && fabs(cimag(dFn[1]))<fabs(cimag(tdqG[1]))*eps)
				&& (fabs(creal(dFn[2]))<fabs(creal(tdqG[2]))*eps && fabs(cimag(dFn[2]))<fabs(cimag(tdqG[2]))*eps)	)break;
		}
	}
	if(n2==nc+EW_LIMIT){
		*err=-2;
		return 0; // summation limit
	}

	// n2 < nc
	for(n2=nc-1;n2>nc-EW_LIMIT;n2--){
		cc+=1;
		*err=ew_qG2_cFn(&Fn,dFn,n1,n2,r,vg1,vg2,eps,veps,qd); if(*err<0) return 0;
		*tqG+=Fn;
		tdqG[0]+=dFn[0];
		tdqG[1]+=dFn[1];
		tdqG[2]+=dFn[2];

		if(*err>0){
			if(  (fabs(creal(Fn))<fabs(creal(*tqG))*eps && fabs(cimag(Fn))<fabs(cimag(*tqG))*eps)
				&& (fabs(creal(dFn[0]))<fabs(creal(tdqG[0]))*eps && fabs(cimag(dFn[0]))<fabs(cimag(tdqG[0]))*eps)
				&& (fabs(creal(dFn[1]))<fabs(creal(tdqG[1]))*eps && fabs(cimag(dFn[1]))<fabs(cimag(tdqG[1]))*eps)
				&& (fabs(creal(dFn[2]))<fabs(creal(tdqG[2]))*eps && fabs(cimag(dFn[2]))<fabs(cimag(tdqG[2]))*eps)	)break;
		}
	}
	if(n2==nc-EW_LIMIT){
		*err=-2;
		return 0; // summation limit
	}

	return cc;
}

int ew_qG2_cFn(double complex *Fn,double complex *dFn,int n1,int n2,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd)
{
	double complex ce,w,ff;
	double an[2],bn2,bn,arg,azve,bnaz,ebz,ef,f1,f2;

	an[0]=-qd->vk[0]+(double)n1*vg1[0]+(double)n2*vg2[0];
	an[1]=-qd->vk[1]+(double)n1*vg1[1]+(double)n2*vg2[1];
	bn2=qd->k*qd->k-an[0]*an[0]-an[1]*an[1];
	arg=an[0]*r[0]+an[1]*r[1];
	ce=cos(arg)+I*sin(arg);

	if(bn2>0.0){
		bn=sqrt(bn2);
		azve=r[2]/(2.0*veps);
		arg=bn*r[2];
		ef=exp(bn2*veps*veps-azve*azve);
		ff=cos(arg)+I*sin(arg);
		w=Faddeeva_w(bn*veps+I*azve,eps);

		*Fn=2.0*ce/bn*(I*ef*cimag(w)+ff);
		dFn[0]=*Fn*an[0];
		dFn[1]=*Fn*an[1];
		dFn[2]=2.0*ce*(ef*creal(w)-ff);

		return 0;
	}
	else if(bn2<0.0){
		bn=sqrt(-bn2); // |bn|
		azve=r[2]/(2.0*veps);
		bnaz=bn*r[2];
		ebz=exp(bnaz);
		f1=ebz*erfc(bn*veps+azve);
		f2=1.0/ebz*erfc(bn*veps-azve);

		*Fn=ce/(I*bn)*(f1+f2);
		dFn[0]=*Fn*an[0];
		dFn[1]=*Fn*an[1];
		dFn[2]=ce*(f1-f2);

		return 1;
	}
	else { // bn=0.0;
		return -3;
	}

}

