#include<affine_transformation.h>

// Scales the data in the incoming array to [0, 1].
void MathTools::affine_transformation(int n, int dim, double in[], double out[])
{

  Eigen::MatrixXd input(n, dim);
  Eigen::MatrixXd output(n, dim);

  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mat;

  input  = Eigen::Map<Mat>(in, n, dim);
  //output = Eigen::Map<Mat>(out, n, dim);

  // get the bounding box for each dimension
  Eigen::VectorXd min = input.colwise().minCoeff();
  Eigen::VectorXd max = input.colwise().maxCoeff();

  // scale input matrix to [0, 1]
  for(int i=0; i<input.rows(); ++i) {
  
    // scale each dimension
    for(int j=0; j<input.cols(); ++j) {

      if(std::abs(max(j) - min(j)) < 1e-8) {
        output(i, j) = 0.0;
      	continue;
      }
      output(i, j) = (input(i, j) - min(j)) / (max(j) - min(j));

    }

  }

  // update output
  for(int i=0; i<n; ++i)
    for(int j=0; j<dim; ++j)
      out[i*dim + j] = output(i, j);

}
