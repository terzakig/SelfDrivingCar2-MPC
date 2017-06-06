#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"

using namespace std;
using namespace Eigen;

// Evaluate a polynomial.
inline double polyeval(const Eigen::VectorXd &coeffs, double x) 
{
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) 
  {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
inline Eigen::VectorXd polyfit(const Eigen::VectorXd &xvals, 
			       const Eigen::VectorXd &yvals,
			       int order) 
{
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) 
  {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) 
  {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}


inline double polyderiveval(const Eigen::VectorXd &coeffs, double x) 
{
  assert(coeffs.size() > 0 && "zero-length coefficients!");
  double p = 1;
  double d = 0;
  for (int i = 1; i < coeffs.size(); i++)
  {
    d += i * coeffs[i] * p;
    p *= x;
  }
  return d;
}

// transform a bunch of points to the car's local coordinate frame
// and return them as a pair of Eigen vectors
inline std::pair<VectorXd, VectorXd > TransformPoints2Local(const std::vector<double> &xvals, 
							    const std::vector<double> &yvals,
							    double px,
							    double py,
							    double phi
  )
{
    assert(xvals.size() == yvals.size() && " Size of coordinate vectors should be the same!");

    VectorXd tx(xvals.size());
    VectorXd ty(yvals.size());
    
    for (int i = 0; i < xvals.size(); i++) 
    {
	double x =  cos(phi) * (xvals[i] - px) + sin(phi) * (yvals[i] - py);
	double y = -sin(phi) * (xvals[i] - px) + cos(phi) * (yvals[i] - py);
	tx[i] = x;
	ty[i] = y;	
	
    }
    
    return std::pair<VectorXd, VectorXd>(tx, ty);
}

class MPC {
 public:
  
  
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(const Eigen::VectorXd &state, const Eigen::VectorXd &coeffs);
};

#endif /* MPC_H */
