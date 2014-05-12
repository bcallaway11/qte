#include <RcppGSL.h>
#include <iostream>
#include <stdio.h>
#include <R.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_pow_int.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>

using namespace Rcpp;

/**
 * @file qtecpp.cpp
 *
 * This is the C++ code that extends (mainly improves speed)
 * for R QTE package.
 *
 * @author Brantly Callaway
 * 
 * @version 1.0
 *
 */

// [[Rcpp::depends(RcppGSL)]]
/*
// [[Rcpp::export]]
NumericVector foo(RcppGSL::vector<double> x, RcppGSL::vector<double> y,
 Function f) {
    RcppGSL::vector<double> v(1);
    printf("%lu", x.size());
    v[0] = 1.23;
    gsl_interp_accel *acc 
      = gsl_interp_accel_alloc ();
    gsl_spline *spline 
      = gsl_spline_alloc (gsl_interp_cspline, 50);
    double *xarray = gsl_vector_ptr(x,0);
    double *yarray = gsl_vector_ptr(y,0);
    int j;
    for(j=0; j<10;j++) {
      printf("%g\n",xarray[j]);
    }
    double ret = v[0];
    //double *xarray = &x.front();
    //double *yarray = &y.front();
    gsl_spline_init(spline, xarray, yarray, 50);
    double retty = gsl_spline_eval(spline, 12, acc);
    v.free();
    x.free();
    y.free();
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
//    return xarray[4];
    return 0;
}
*/



/**
 * partial1CopualCPP is the C++ implementation of the partial
 * derivative of the copula function with respect to its first 
 * argument.  We need to do this to simulate random draws
 * from the joint distribution following the procedure 
 * outlined in Nelson (2005).
 *
 * @param[in] u NumericVector of values for the first copula function argument.
 * Here u should be a scalar (=> probably should change the parameter type to 
 * double), but need to test that this doesn't break anything.
 * @param[in] v NumericVector of values for the second copula function argument
 * @param[in] h double the step size used.  It is the same as used throughout 
 * the method call
 * @param[in] copfun R function that is generated in the panelDiD function
 * @return NumericVector for the derivative of the copula function in the 
 * direction of the first argument
 */
// [[Rcpp::export]]
NumericVector partial1CopulaCPP(NumericVector u,
  NumericVector v, double h, Function copfun) {
    
    //RNGScope scope;
    //Environment base("package:base");
      
    NumericVector out(v.size());
    NumericVector temp1;
    NumericVector temp2;
    int i;
    temp1 = copfun(u+h,v);
    temp2 = copfun(u-h,v);
    out = (temp1-temp2)/(2*h);//(copfun(u+h,v)-copfun(u-h,v))/(2*h);
    //out = 1*(out>1) + out*(out<=1);
    return out;
}

/**
 *getListPartialQuantCPP is a wrapper for partial1.copula.  It takes in a 
 *NumericVector of u's, and then calls partial1.copulaCPP for each
 *of those u's individually.
 *
 *@param u NumericVector of length = #probevals
 *@param v NumericVector of grid points to evaluate copula function at
 *@param h double step size so don't get too close to boundaries
 *@param copfun points to copula function generated previously
 *
 *@return List of partial derivatives of copula function for each of the 
 *values of u passed in.
 */
// [[Rcpp::export]]
List getListPartialQuantCPP(NumericVector u, NumericVector v, double h,
 	Function copfun) {

  //RNGScope scope;
  //Environment base("package:base");
  //Function print = base["print"];

  int n = u.size();
  List retList(n);
  for (int i=0; i<n; i++) {
    retList(i) = partial1CopulaCPP(NumericVector::create(u(i)),v,h,copfun);
  }
  return retList;
}

/**
 * Not implemented yet, but potential performance gains by 
 * implementing the joint distribution in C++ rather than R.
 */
// [[Rcpp::export]]
NumericVector getJointCPP(NumericVector x,
  NumericVector y) {

  //int xsize;
  //int ysize;

  //xsize = x.size();
  //ysize = y.size();

  //int i, j;

  //for (i=0; i<xsize; i++) {
  //  for (j=0; j<ysize; j++) {
      //do nothing
  //  }
  //}
  //return 0;
}

/**
 * getJointUVCPP takes in two points and returns the value
 *
 * of the empirical joint CDF (of x1 and y1) for those two values
 * @param u double a point in x
 * @param v double a point in y
 * @param x1 NumericVector a set of datapoints of x
 * @param y1 NumericVector a set of datapoints of y
 *
 * @return NumericVector scalar value in joint distribution
 */
// [[Rcpp::export]]
NumericVector getJointUVCPP(double u, NumericVector v,
 NumericVector x1, NumericVector y1) {

  //this part does empirical joint distribution
  NumericVector out = NumericVector(v.size());
  int vsize = v.size();
  int ysize = y1.size();
  double thisout;
  for (int i=0; i<vsize; i++) {
    thisout = 0;
    for (int j=0; j<ysize; j++) {
      thisout += (x1(j)<u && y1(j)<v(i));
    }
    out(i) = thisout/ysize;
  }
  //add code for computing kernel density estimate.
  return out;
}

/**
 * Wrapper for using GSL package to compute quantiles
 *
 *@param x RcppGSL vector of function values
 *@param probs NumericVector of quantiles to compute
 *
 *@return NumericVector of quantiles.  Return will be of same length as probs.
 *
 */
// [[Rcpp::export]]
NumericVector quantileCPP(RcppGSL::vector<double> x,
 NumericVector probs) {

    //double *xarray = gsl_vector_ptr(x,0);
    //double *yarray = gsl_vector_ptr(y,0);

  gsl_sort_vector(x);
  //gsl_sort_vector(y);
  double *xarray = gsl_vector_ptr(x,0);
  int n = x.size();
  int probsize = probs.size();
  NumericVector out = NumericVector(probsize);

  for (int i=0; i<probsize; i++) {
    out(i) = gsl_stats_quantile_from_sorted_data(xarray,1,n,probs(i));
  }
  
  //x.free();

  return(out);
}


/**
 * Overloads quantileCPP so that it can be called
 * for a NumericVector rather than GSL vector.
 *
 *@param x NumericVector
 *@param probs NumericVector to compute quantiles for
 *
 *@return NumericVector of same length as probs.
 *
 */
NumericVector quantileCPP(NumericVector x,
 NumericVector probs) {
  
  int xsize = x.size();
  gsl_vector *v = gsl_vector_alloc(xsize);

  for (int i=0; i<xsize; i++) {
    gsl_vector_set(v,i,x(i));
  }

  NumericVector out = quantileCPP(v,probs);

  gsl_vector_free(v);

  return out;

}

/**
 * getListQuantilesCPP takes the List of 1st partial derivatives
 * created in the previous steps, and gets (#probevals = length(t))
 * random quantiles from them.  The purpose is for simulating from
 * joint distribution of the change and the intial outcomes.
 * 
 *@param partialvals List of vectors containing the 1st partial 
 * derivative of the copula function for the random uniforms 
 * generated previously
 *@param t NumericVector for random quantiles between 0 & 1
 */
// [[Rcpp::export]]
List getListQuantilesCPP(List partialvals,
 NumericVector t) {
  
  //RNGScope scope;
  //Environment base("package:base");
  //Function print = base["print"];

  int n = t.size();
  List out = List(n);
  //RcppGSL::vector<double> temp;

  //use these just to hold variables for next fun
  NumericVector temp1;
  NumericVector temp2;
  
  for (int i=0; i<n; i++) {
    temp1 = partialvals(i);
    temp2 = t(i);
    //print(temp);
    out(i) = quantileCPP(temp1, temp2);
  }
  return out;
}

int retInt(int inty) {
 return 5;
}

int (*retFun(int i))(int j) {
  return &retInt;
}


// [[Rcpp::export]]
int testy() {
  int (*out)(int) = retFun(20);
  return out(9);
}

// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i] / n;
  }
  return total;
}

/*
// [[Rcpp::export]]
NumericVector getFuncVal(NumericVector uSeq, NumericVector vSeq, double h) {
	int uLength = uSeq.size()
	int vLength = vSeq.size()
	for (int i=0; i<uLength, i++) {
		for (int j=0; j<vLength, j++) {
		
		}
	}
}
*/


int main() {

	double x = 5.0;
	return 0;
}


