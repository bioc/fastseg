#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <stack>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

using namespace std;


extern "C" SEXP segment(SEXP xS, SEXP epsS, SEXP deltaS, SEXP maxIntS,
		SEXP minSegS, SEXP squashingS, SEXP cyberWeightS) {
	int j, d, ll;
	long n=LENGTH(xS);
	long i;

	double eps = REAL(epsS)[0];
	//Rprintf("eps: %lf\n", eps);

	double cyberWeight = REAL(cyberWeightS)[0];
	//Rprintf("cyberWeight: %lf\n", cyberWeight);
	if (cyberWeight < 0.01){
		cyberWeight = 0.01;
	}

	int delta= INTEGER(deltaS)[0];
	int maxInt= INTEGER(maxIntS)[0];
	int minSeg= INTEGER(minSegS)[0];
	int squashing= INTEGER(squashingS)[0];

	double* x = REAL(xS);
	double *partialSumValues=(double *) R_alloc(n, sizeof(double));
	double *partialSumSquares=(double *) R_alloc(n, sizeof(double));
	double *pValue=(double *) R_alloc(n, sizeof(double));
	//long *leftBorders=(long *) R_alloc(n, sizeof(long));
	//long *rightBorders=(long *) R_alloc(n, sizeof(long));

	SEXP x_RET;
	PROTECT(x_RET = Rf_allocVector(REALSXP, n));
	double *xx=REAL(x_RET);

	SEXP savedStatistic_RET;
	PROTECT(savedStatistic_RET = Rf_allocVector(REALSXP, n));
	double *savedStatistic=REAL(savedStatistic_RET);

	//SEXP hist_RET;
	//PROTECT(hist_RET = Rf_allocVector(REALSXP, n));
	//double *hist=REAL(hist_RET);

	SEXP leftright_RET;
	PROTECT(leftright_RET = Rf_allocVector(INTSXP, n));
	int *leftright=INTEGER(leftright_RET);



	double globalMean,globalSd,diff,M2,globalVariance;
	double oldStatistic, meanLeft,meanRight,varLeft,varRight;
	double newStatistic,meanDiff,maxStatistic,DOF,a,b,eps1;
	double newPValue, maxPValue,oldPValue,maxIdx;
	//double newStatisticBptLeft,newStatisticBptRight,
	double beta,nn;

	eps1 = 1e-15;
	globalMean=0;
	globalSd=0;
	partialSumValues[0]=x[0];
	partialSumSquares[0]=x[0]*x[0];
	M2 = 0;
	for (i=0;i<n;i++){
		diff = x[i] - globalMean;
		globalMean = globalMean + diff/(i+1);
		//Rprintf("Mean: %lf\n", globalMean);
		M2 = M2 + diff*(x[i] - globalMean);
		//Rprintf("M2: %lf\n", M2);
		//Rprintf("x: %lf\n", x[i]);
		if (i>0){
			partialSumValues[i]=partialSumValues[i-1]+x[i];
			partialSumSquares[i]=partialSumSquares[i-1]+x[i]*x[i];
		}
		//Rprintf("PartialSumValues: %lf\n", partialSumValues[i]);
		//Rprintf("PartialSumSquares: %lf\n", partialSumSquares[i]);

		//hist[i]=2;
		xx[i]=x[i];


	}
	globalVariance = M2/(n-1);

	if (squashing > 0){
		// Experimental - will be completely removed in the next version.

		/*
		//beta = -log(2.0/1.8-1)/((double) squashing * sqrt(globalVariance));
		//Rprintf("Beta: %lf\n", beta);
		beta = (double) squashing;

		for (i=0;i<n;i++){
			x[i]=(2/(1+exp(-1/beta*((x[i]-globalMean)/sqrt(globalVariance))))-1);
		}
		globalMean=0;
		globalSd=0;
		partialSumValues[0]=x[0];
		partialSumSquares[0]=x[0]*x[0];
		M2 = 0;
		for (i=0;i<n;i++){
			diff = x[i] - globalMean;
			globalMean = globalMean + diff/(i+1);
			//Rprintf("Mean: %lf\n", globalMean);
			M2 = M2 + diff*(x[i] - globalMean);
			//Rprintf("M2: %lf\n", M2);
			//Rprintf("x: %lf\n", x[i]);
			if (i>0){
				partialSumValues[i]=partialSumValues[i-1]+x[i];
				partialSumSquares[i]=partialSumSquares[i-1]+x[i]*x[i];
			}
			//Rprintf("PartialSumValues: %lf\n", partialSumValues[i]);
			//Rprintf("PartialSumSquares: %lf\n", partialSumSquares[i]);

		}
		globalVariance = M2/(n-1);
		//Rprintf("Squashing values.\n");
		*/
	}

	//Rprintf("Using original values.\n");



	if (globalVariance < eps1){
		//Rprintf("Global Variance is zero!\n");
		globalVariance = eps1;
	}


	i = 0;
	while (i < (n-1)){
		//if (i > minSeg && i < n-minSeg-1){
		if ( fabs(x[i+1]-x[i])> eps && i > minSeg && i < n-minSeg-1){
			j = minSeg;
			d = 0;
			oldStatistic=0.0;
			maxStatistic=0.0;
			newPValue=0.0;
			maxPValue=0.0;
			oldPValue=0.0;
			maxIdx=0;

			while (d<=delta && j<=maxInt && (i+j+1) < n && (i-j-1)>=0){
				nn = ((double) j)+cyberWeight-1.0;

				meanLeft=(partialSumValues[i]-partialSumValues[i-j-1])/(j+1);
				varLeft=((partialSumSquares[i]-partialSumSquares[i-j-1])-(j+1)*meanLeft*meanLeft);
				varLeft=((varLeft)+cyberWeight*globalVariance)/(nn);

				//Rprintf("PartialSumSquaresLeft: %lf\n", (partialSumSquares[i]-partialSumSquares[i-j-1]));
				//Rprintf("MeanLeft: %lf\n", meanLeft);
				//Rprintf("VarLeft: %lf\n", varLeft);


				meanRight=(partialSumValues[i+j+1]-partialSumValues[i])/(j+1);
				varRight=((partialSumSquares[i+j+1]-partialSumSquares[i])-(j+1)*meanRight*meanRight);
				varRight=((varRight)+cyberWeight*globalVariance)/(nn);

				//Rprintf("PartialSumSquaresRight: %lf\n", (partialSumSquares[i+j+1]-partialSumSquares[i]));
				//Rprintf("MeanRight: %lf\n", meanRight);
				//Rprintf("VarRight: %lf\n", varRight);

				meanDiff=(meanLeft-meanRight);
				newStatistic=fabs(meanDiff)/sqrt(varLeft/(nn+1.0)+varRight/(nn+1.0)+eps1);

				leftright[i]=1;


				/*
				// bptLeft
				meanLeft=(partialSumValues[i]-partialSumValues[i-j-1])/(j+1);
				varLeft=((partialSumSquares[i]-partialSumSquares[i-j-1])-(j+1)*meanLeft*meanLeft);
				varLeft=(varLeft+cyberWeight*globalVariance)/(nn);


				meanRight=(partialSumValues[i+j]-partialSumValues[i])/(j);
				varRight=((partialSumSquares[i+j]-partialSumSquares[i])-(j)*meanRight*meanRight);
				varRight=(varRight+cyberWeight*globalVariance)/(nn-1.0);


				meanDiff=(meanLeft-meanRight);
				newStatisticBptLeft=fabs(meanDiff)/sqrt(varLeft/(nn+1.0)+varRight/(nn)+eps1);

				// bptRight
				meanLeft=(partialSumValues[i-1]-partialSumValues[i-j-1])/(j);
				varLeft=((partialSumSquares[i-1]-partialSumSquares[i-j-1])-(j)*meanLeft*meanLeft);
				varLeft=(varLeft+cyberWeight*globalVariance)/(nn-1.0);


				meanRight=(partialSumValues[i+j]-partialSumValues[i-1])/(j+1);
				varRight=((partialSumSquares[i+j]-partialSumSquares[i-1])-(j+1)*meanRight*meanRight);
				varRight=(varRight+cyberWeight*globalVariance)/(nn);


				meanDiff=(meanLeft-meanRight);
				newStatisticBptRight=fabs(meanDiff)/sqrt(varLeft/nn+varRight/(nn+1)+eps1);

				if (newStatisticBptLeft > newStatisticBptRight){
					newStatistic = newStatisticBptLeft;
					leftright[i]=0;
				} else{
					newStatistic = newStatisticBptRight;
					leftright[i]=1;
				}
				 */



				a=varLeft/(nn+1.0);
				b=varRight/(nn+1.0);
				DOF=(a+b)*(a+b)/((a*a)/nn +(b*b)/nn );

				newPValue= -pt(newStatistic,DOF,0,1);


				if (newPValue>maxPValue){
					maxStatistic=newStatistic;
					maxPValue=newPValue;
					maxIdx=((double) j);
				}

			
				//Rprintf("NewStatistic: %lf\n", newStatistic);

				if (newPValue>oldPValue){
					d = 0;
				} else {
					d = d+1;
				}
				oldPValue = newPValue;
				j = j+1;
			}

			//starts[i] = i;
			pValue[i]=maxPValue;
			//leftBorders[i]=i-((long) maxIdx)-1;
			//rightBorders[i]=i+((long) maxIdx)+1;
			//hist[leftBorders[i]]++;
			//hist[rightBorders[i]]++;

			i=i+1;
		} else{
			//starts[i] = i;
			pValue[i]=0; //actually log p value
			//savedDOF[i]=i;
			leftright[i]=-1;
			i=i+1;
		}

	}
	pValue[n-1]=0;
	leftright[n-1]=-1;
	// Determine local maxima
	i=0;
	if (minSeg > 2){
		while(i<n){
			savedStatistic[i]=pValue[i];
			ll = ((int) floor(((double) minSeg)/ 2.0));
			//ll = minSeg;

			if ((i -ll > 0) &&  (i+ll < n)){
				for (j=1;j<=ll;j++){
					if (pValue[i-j] > savedStatistic[i] || pValue[i+j] > savedStatistic[i] ){
						savedStatistic[i]=0;
					}
				}
			}

			//hist[i]=(savedStatistic[i]*log2(hist[i]))/2;
			i=i+1;
		}
	} else{
		for (i=0;i<n;i++){
			savedStatistic[i]=pValue[i];
			//hist[i]=(savedStatistic[i]*log2(hist[i]))/2;

		}
	}

	SEXP namesRET;
	PROTECT(namesRET = Rf_allocVector(STRSXP, 3));
	SET_STRING_ELT(namesRET, 0, Rf_mkChar("x"));
	SET_STRING_ELT(namesRET, 1, Rf_mkChar("stat"));
	//SET_STRING_ELT(namesRET, 2, Rf_mkChar("stat2"));
	SET_STRING_ELT(namesRET, 2, Rf_mkChar("leftright"));

	SEXP RET;
	PROTECT(RET = Rf_allocVector(VECSXP, 3));
	SET_VECTOR_ELT(RET, 0, x_RET);
	SET_VECTOR_ELT(RET, 1, savedStatistic_RET);
	//SET_VECTOR_ELT(RET, 2, hist_RET);
	SET_VECTOR_ELT(RET, 2, leftright_RET);
	Rf_setAttrib(RET, R_NamesSymbol, namesRET);
	UNPROTECT(5);
	return(RET);

}
