/* Copyright (C) 1998
   Berwin A Turlach <bturlach@stats.adelaide.edu.au> */

/* This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public License
   as published by the Free Software Foundation; either version 2 of
   the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details. */

/* You should have received a copy of the GNU Library General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston,
   MA 02111-1307, USA. */

#ifndef BT_LASSO_H
#define BT_LASSO_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifdef __cplusplus
extern "C" {
#endif

#define Calloc(n,t)       malloc((n)*sizeof(t))
#define Realloc(p,n,t)    realloc(p,(n)*sizeof(t))
#define Free(p)           free(p), p=NULL
#define Memcpy(p,q,n)     memcpy(p,q,(n)*sizeof(*(p)))


void lasso(double *x, int *pn, int *pm, double *pt,
		      double *beta, double *y, double *yhat1, double *r,
		      double *lagrangian, int *psuc,  int *pverb, int *pas_sub);
void mult_lasso(double *x, int *pn, int *pm, double *pt, int *pl,
			   double *beta, double *y, double *yhat1, double *r,
			   double *lagrangian, int *psuc, int *pverb);

#define TRUE 1
#define FALSE 0
#define RMAT(i,j) (*(rmat+(j)*((j)+1)/2+(i)))
#define RCOL(j)   (rmat+(j)*((j)+1)/2)
#define QCOL(j)   (qmat+(j)*q_nrow)
#define QR_CHUNK 10

#ifdef __cplusplus
}
#endif

#endif
