#include <weighted_interp.hpp>
#include <iostream>
#include <cmath>

namespace MathTools {

double vec_distance(int n, double v1[], double v2[])
{
  double sum = 0.0;
  for(int i=0; i<n; ++i)
    sum += pow(v1[i]-v2[i],2);

  return sqrt(sum);
}

void weighted_interp(int dim, int nd, double xd[], double r0,
                     void phi(int n, double r[], double r0, double v[]), double fd[],
                     int nq, double xq[], 
                     double fq[])
{
  for(int i=0; i<nq; ++i) {
    double sd = 0.0;
    double sn = 0.0;
    for(int j=0; j<nd; ++j) {
      double r = vec_distance(dim, xq+(i*dim), xd+(j*dim));
      double w = 0.0;
      phi(1, &r, r0, &w);
      sd += w;
      sn += w*fd[j];
    }
    fq[i] = sn/sd;
  }
}

}
