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


// [[Rcpp::export]]
List getFuncValCPP(NumericVector u, NumericVector v, double h,
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

// [[Rcpp::export]]
NumericVector getJointUVCPP(double u, NumericVector v,
 NumericVector x1, NumericVector y1) {

  NumericVector out = NumericVector(v.size());
  int vsize = v.size();
  int ysize = y1.size();
  double thisout;
  for (int i=0; i<vsize; i++) {
    thisout = 0;
    //out(i) = sum(1*(x1<u & y1<v[i]))/length(y1)
    for (int j=0; j<ysize; j++) {
      thisout += (x1(j)<u && y1(j)<v(i));
      //printf("%g",thisout);
    }
    //printf("\n\n%g",thisout);
    out(i) = thisout/ysize;
  }
  return out;
}

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


/*Overloads quantileCPP so that is can be 
called for a NumericVector*/


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

//funcvals and t should be of same length 
// [[Rcpp::export]]
List getListQuantilesCPP(List funcvals,
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
    temp1 = funcvals(i);
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


