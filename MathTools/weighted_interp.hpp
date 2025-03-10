using namespace std;

namespace MathTools {

// helper functions
double vec_distance(int n, double v1[], double v2[]);

//interpolation functions
void weighted_interp(int dim, int nd, double xd[], double r0, 
  void phi(int n, double r[], double r0, double v[]), double fd[],
  int nq, double xq[], double fq[]);

}
