/*
 * d3_qpgf_d2.h
 *
 *  Created on: Jun 11, 2019
 *      Author: ohta
 */

#ifndef SRC_D3_QPGF_D2_H_
#define SRC_D3_QPGF_D2_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "Faddeeva.h"

// limit number
#define S_LIMIT 15

// Fourier domain summation
#define FD_LIMIT 100

// Ewald method
#define EW_LIMIT 100
#define EW_EXPM 3.0


// struct
typedef struct qp_data_v2{
  double k;             // wave number
  double vk[3];         // wave vector
  double vd1[2],vd2[2]; // 2d-lattice vector ( z-component is zero )
}QPD2;


int d3hm_qpgf_d2_fd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd);
/* Fourier domain summation
      qG : value of quasi-periodic green function
    *dqG : value of differential quasi-periodic green function
           dqG[0]=dqG/dx, dqG[1]=dqG/dy, dqG[2]=dqG/dz
       r : (x,y,z)
           r[0]=x, r[1]=y, r[2]=z
     eps : requested relative error
     *qd : pointer of quasi-periodic data (struct QPD2)

   return code > 0 : normal termination. total summation number
               =-1 : abnormal termination. reached to summation limit
                -2 : abnormal termination. \beta_n == 0.0 ( divergence )
                -3 : abnormal termination. determinant d == 0.0
                -4 : abnormal termination. z == 0.0
*/

int d3hm_qpgf_d2_ew(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd);
/* Ewald method
      qG : value of quasi-periodic green function
    *dqG : value of differential quasi-periodic green function
           dqG[0]=dqG/dx, dqG[1]=dqG/dy, dqG[2]=dqG/dz
       r : (x,y,z)
           r[0]=x, r[1]=y, r[2]=z
     eps : requested relative error
     *qd : pointer of quasi-periodic data (struct QPD2)

   return code > 0 : normal termination. total summation number
               =-1 : abnormal termination. qG1 reached to summation limit
                -2 : abnormal termination. qG2 reached to summation limit
                -3 : abnormal termination. \beta_n == 0.0
                -4 : abnormal termination. determinant d == 0.0
*/

int d3hm_qpgf_d2(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd);
/* auto select
   return code >  0 : total summation number
               =-11 : Fourier domain summation, abnormal termination. reached to summation limit
                -12 : Fourier domain summation, abnormal termination. \beta_n == 0.0
                -13 : Fourier domain summation, abnormal termination. determinant d == 0.0
                -14 : Fourier domain summation, abnormal termination. z == 0.0
                -21 : Ewald method, abnormal termination. qG1 reached to summation limit
                -22 : Ewald method, abnormal termination. qG2 reached to summation limit
                -23 : Ewald method, abnormal termination. \beta_n == 0.0
                -24 : Ewald method, abnormal termination. determinant d == 0.0
*/


#endif /* SRC_D3_QPGF_D2_H_ */
