#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;




// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;
 
const double T = 5; // seconds look-ahead time horizon
const size_t N = 18;
const double dt = T / N;
const size_t x_start = 0;
const size_t y_start = x_start + N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;

const double v_ref = 30 * ( 1609.34 / 3600.0); 


class FG_eval 
{
 public:
   
   
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(const Eigen::VectorXd &coeffs ) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) 
  {
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;
  
    // The part of the cost based on the reference state.
    for (int i = 0; i < N; i++) 
    {
      fg[0] += 100*CppAD::pow(vars[cte_start + i] , 2);
      fg[0] += 100 * CppAD::pow(vars[epsi_start + i] , 2);
      fg[0] += CppAD::pow(vars[v_start + i] - v_ref, 2);
    }

    // Minimize the use of actuators.
    for (int i = 0; i < N - 1; i++) 
    {
      fg[0] +=  CppAD::pow(vars[delta_start + i], 2);
      fg[0] +=  20*CppAD::pow(vars[a_start + i], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int i = 0; i < N - 2; i++) 
    {
      fg[0] += 1000000*CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
      fg[0] += 100 * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
    }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

   
    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // The rest of the constraints
    for (int i = 0; i < N - 1; i++) 
    {
      // The state at time t+1 .
      AD<double> x_n = vars[x_start + i + 1];
      AD<double> y_n = vars[y_start + i + 1];
      AD<double> psi_n = vars[psi_start + i + 1];
      AD<double> v_n = vars[v_start + i + 1];
      AD<double> cte_n = vars[cte_start + i + 1];
      AD<double> epsi_n = vars[epsi_start + i + 1];

      // The state at time t.
      AD<double> x_n_1 = vars[x_start + i];
      AD<double> y_n_1 = vars[y_start + i];
      AD<double> psi_n_1 = vars[psi_start + i];
      AD<double> v_n_1 = vars[v_start + i];
      AD<double> cte_n_1 = vars[cte_start + i];
      AD<double> epsi_n_1 = vars[epsi_start + i];

      // Only consider the actuation at time t.
      AD<double> delta_n_1 = vars[delta_start + i];
      AD<double> a_n_1 = vars[a_start + i];

      AD<double> f_n_1 = coeffs[0]; 
      for (int k = 1; k < coeffs.size(); k++)
	f_n_1 += coeffs[k] * CppAD::pow(x_n_1, k);
      
      AD<double> df_n_1 = coeffs[1];
      for (int k = 2; k < coeffs.size(); k++)
	df_n_1 += k * coeffs[k] * CppAD::pow(x_n_1, k-1);
      
      AD<double> psides_n_1 = CppAD::atan(df_n_1);

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      fg[2 + x_start + i] = x_n - ( x_n_1 + v_n_1 * CppAD::cos(psi_n_1) * dt );
      fg[2 + y_start + i] = y_n - ( y_n_1 + v_n_1 * CppAD::sin(psi_n_1) * dt );
      fg[2 + psi_start + i] = psi_n - ( psi_n_1 + v_n_1 * delta_n_1 / Lf * dt );
      fg[2 + v_start + i] = v_n - ( v_n_1 + a_n_1 * dt );
      fg[2 + cte_start + i] =
          cte_n - ( (f_n_1 - y_n_1) + ( v_n_1 * CppAD::sin(epsi_n_1) * dt ) );
      fg[2 + epsi_start + i] =
          epsi_n - ( ( psi_n_1 - psides_n_1 ) + v_n_1 * delta_n_1 / Lf * dt );
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
			   
MPC::~MPC() {}

vector<double> MPC::Solve(const Eigen::VectorXd &state, const Eigen::VectorXd &coeffs) {
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  
  double x_0 = state[0];
  double y_0 = state[1];
  double psi_0 = state[2];
  double v_0 = state[3];
  double cte_0 = state[4];
  double epsi_0 = state[5];

  
  
  // number of independent variables
  // N timesteps == N - 1 actuations
  size_t n_vars = N * 6 + (N - 1) * 2;
  // Number of constraints
  size_t n_constraints = N * 6;

  // Zeroing all (solver) variables
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) 
    vars[i] = 0.0;
 
    
  
  // Set the initial variable values
  vars[x_start] = x_0;
  vars[y_start] = y_0;
  vars[psi_start] = psi_0;
  vars[v_start] = v_0;
  vars[cte_start] = cte_0;
  vars[epsi_start] = epsi_0;

   
  // Lower and upper limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; i++) 
  {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; i++) 
  {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start; i < n_vars; i++) 
  {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for constraints
  // All of these should be 0 except the initial
  // state indices.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) 
  {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  
  
  // A. Assigning lower and upper bounds for the INITIAL state variables in the solver.
  //    They are both set to the same (initial) values.
  constraints_lowerbound[x_start] = x_0;
  constraints_lowerbound[y_start] = y_0;
  constraints_lowerbound[psi_start] = psi_0;
  constraints_lowerbound[v_start] = v_0;
  constraints_lowerbound[cte_start] = cte_0;
  constraints_lowerbound[epsi_start] = epsi_0;

  constraints_upperbound[x_start] = x_0;
  constraints_upperbound[y_start] = y_0;
  constraints_upperbound[psi_start] = psi_0;
  constraints_upperbound[v_start] = v_0;
  constraints_upperbound[cte_start] = cte_0;
  constraints_upperbound[epsi_start] = epsi_0;

  // B. Configuring objective and constraints
  FG_eval fg_eval(coeffs);

  // options
  std::string options;
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  options += "Numeric max_cpu_time          0.05\n";

  // Solution object
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>( options, 
					 vars, 
					 vars_lowerbound, 
					 vars_upperbound, 
					 constraints_lowerbound,
					 constraints_upperbound, 
					 fg_eval, 
					 solution
				       );


  // Check if the solver was successful
  bool ok = true;
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;
  // building a vector contaoining some values to display...
  
  /*return {solution.x[x_start + 1],   solution.x[y_start + 1],
          solution.x[psi_start + 1], solution.x[v_start + 1],
          solution.x[cte_start + 1], solution.x[epsi_start + 1],
          solution.x[delta_start],   solution.x[a_start]};
	  */
  size_t step = 1;
  vector<double> ret;
  for (int i = 0; i < N-2; i += step) 
  {
    ret.push_back( solution.x[x_start + i + 1] );
    ret.push_back( solution.x[y_start + i + 1] );
    
  }
  
  ret.push_back( solution.x[delta_start] );
  ret.push_back( solution.x[a_start] );
  return ret;
}

