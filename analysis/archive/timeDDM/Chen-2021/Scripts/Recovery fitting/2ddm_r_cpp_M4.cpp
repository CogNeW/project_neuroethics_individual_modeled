#include <Rcpp.h>
#include <RcppParallel.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <time.h>

using namespace Rcpp;
using namespace RcppParallel;

typedef boost::mt19937                     ENG;    // Mersenne Twister
typedef boost::normal_distribution<double> DIST;   // Normal Distribution
typedef boost::variate_generator<ENG,DIST> GEN;    // Variate generator



// [[Rcpp::depends(RcppParallel)]]
struct ddm2w : public Worker {

  // input vectors/matrices to read from
  const double d_v;
  const double d_h;
  const double d_a;
  const double thres;
  const double nDT;
  const double tIn;
  const double bias;
  const double vd;
  const double hd;
  const double ad;
  const double sd_n;
  GEN gen;

  // output vector to write to
  RVector<double> vecOut;

  // initialize from Rcpp input and output matrices/vectors (the RMatrix/RVector class
  // can be automatically converted to from the Rcpp matrix/vector type)
  ddm2w(const double d_v, const double d_h, const double d_a, const double thres, const double nDT, const double tIn, const double bias, const double vd, const double hd, const double ad, const double sd_n, NumericVector vecOut , GEN gen)
    : d_v(d_v), d_h(d_h), d_a(d_a), thres(thres), nDT(nDT), tIn(tIn), bias(bias), vd(vd), hd(hd), ad(ad), sd_n(sd_n), gen(gen), vecOut(vecOut) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    double T = 5.2, dt = 0.001, lt;
    lt = (int)(T/dt);

    std::vector<double> ven_tMP(lt,1);
    std::vector<double> ven_tA(lt,1);
    int aux = std::abs((int)(tIn/dt));

    if (tIn > 0) {
      for (int t=0; t<aux; t++) {
        ven_tMP[t] = 0;
      }
    }
    else if (tIn < 0) {
      for (int t=0; t<aux; t++) {
        ven_tA[t] = 0;
      }
    }


    for (std::size_t i = begin; i < end; i++) {
      vecOut[i] = T;
      double X = bias*thres;
      int flag = 0;
      double cont = 0;
      double noise = 0;

      while (flag==0 && cont<lt) {

        noise=gen()*sqrt(dt);
        X = X + ((d_v*vd + d_h*hd)*ven_tMP[cont] + d_a*ad*ven_tA[cont])*dt + noise;

        if (X > thres) {
          flag=1;
          vecOut[i] = nDT + cont*dt;
        }
        else if (X < 0) {
          flag=1;
          vecOut[i] = -nDT -cont*dt;
        }
        cont++;

      }
    }
  }
};


// [[Rcpp::export]]
NumericVector ddm2_parallel(double d_v, double d_h, double d_a, double thres, double nDT, double tIn, double bias, double vd, double hd, double ad, double sd_n, unsigned int N) {

  //const double sd_n = 1.4;
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  ENG  eng;
  eng.seed(time.tv_nsec);
  DIST dist(0,sd_n);
  GEN  gen(eng,dist);

  //output vector
  NumericVector vecOut(N);

  // create the worker
  ddm2w ddm2w(d_v, d_h, d_a, thres, nDT, tIn, bias, vd, hd, ad, sd_n, vecOut, gen);

  // call the worker
  parallelFor(0, N, ddm2w);

  return vecOut;
}
