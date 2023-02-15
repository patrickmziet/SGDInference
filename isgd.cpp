#include <RcppArmadillo.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
using namespace Rcpp;

  
arma::rowvec elementwise_pow(arma::rowvec base, const double p) {
  arma::rowvec result;
  result.copy_size(base);
  for (std::size_t i = 0; i < result.n_elem; ++i) {
    result[i] = std::pow(base[i], p);
  }
  return result;
}


/* defining key functions */
arma::rowvec implicit_sgd(arma::mat dataX, arma::vec dataY, arma::rowvec C, 
													std::string model_name, double gamma, int n_pass, uint32_t seed, bool use_permutation);

/* testing out Cpp multivariate normal generating 
 * Currently not used
 * */
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

/* tr(var) asymptotic ========================================================== */
/* to calculate trace of asymptotic covariance 
 * Currently not used
 * */
// [[Rcpp::export]]
double tr_asymp_var(arma::mat dataX, arma::vec dataY, std::string model_name,
                    double gamma, int n_pass, arma::vec seed_ls) {
  /* dimension parameters */
  int N = dataX.n_rows;
  int p = dataX.n_cols;
	int R = seed_ls.n_elem;
  
  /* updated variables */
  uint32_t seed;
  arma::rowvec theta = arma::zeros<arma::rowvec>(p);
  arma::rowvec U_t   = arma::zeros<arma::rowvec>(p);
  arma::rowvec M_t   = arma::zeros<arma::rowvec>(p);
  arma::rowvec V_t   = arma::zeros<arma::rowvec>(p);
	int idx;

	arma::rowvec C = arma::ones<arma::rowvec>(p);
  /* main loop */
  for(double r=1; r <= R; r++) {
    idx = r-1;
		seed  = seed_ls[idx]; 
    theta = implicit_sgd(dataX, dataY, C, model_name, gamma, n_pass, seed, true);
		U_t   = M_t % M_t; /* element-wise multiplication */
		M_t   = (1/r) * ( (r-1) * M_t + theta);
		V_t   = (1/r) * ( (r-1) * V_t + theta % theta + (r-1) * U_t ) - M_t % M_t;
		/*
		Rcout << "theta: " << theta << std::endl;
		Rcout << "U_t: " << U_t << std::endl;
		Rcout << "M_t: " << M_t << std::endl;
		Rcout << "V_t: " << V_t << std::endl;
		*/
  }

	/* normalize variance */
	V_t = N * V_t;
	Rcout << "tr should be: " << ( (gamma/2.0) * (p * 1.0) ) << std::endl;
	Rcout << "tr is: " << arma::accu(V_t) << std::endl;
	Rcout << "variance is: " << V_t << std::endl;
	return(arma::accu(V_t));
}

// [[Rcpp::export]]
void test_stuff(arma::vec dataY) {
	Rcout << "element wise mult: " << dataY % dataY << std::endl;
	Rcout << "scalar times vec (2.0 * dataY): " << 2.0 * dataY << std::endl;
}



/* ISGD ======================================================================= */
/* define zeroin function */ 
double zeroin(double ax, double bx, double (*f)(double x), double tol, int max_iter);

/* for gaussian glm */
double identity(double predictor) {
	return predictor;
}

/* for binomial glm */
double expit(double predictor) {
	if (predictor > 60) {
		return 1;
	} else if (predictor < -60) {
		return 0;
	} else {
		return std::exp(predictor) / (1 + std::exp(predictor));
	}
}

 /* Lmin_bounds cpp helper function */
// [[Rcpp::export]]
arma::rowvec Lmin_bounds_cpp(arma::mat dataX, arma::vec dataY,
                             std::string model_name,
                             arma::rowvec theta) {
  // dimension parameters
  int N = dataX.n_rows;
  int p = dataX.n_cols;
  
  // chose glm_link
  double (*glm_link_loc) (double); // local version in function to not mess with isgd function
  if (model_name.compare("gaussian") == 0) {
    glm_link_loc = identity;
  } else if (model_name.compare("binomial") == 0) {
    glm_link_loc = expit;
  } else if (model_name.compare("poisson") == 0) {
    glm_link_loc = std::exp; 
  } else {
    throw std::invalid_argument("model_name not gaussian or binomial or poisson");
  }
  
  // Estimate fisher information
  arma::mat J = arma::mat(p, p, arma::fill::zeros);
  double yj;
  arma::rowvec xj;
  arma::rowvec grad_j;
  for (int i=0; i<N; i++) {
    yj = dataY(i);
    xj = dataX.row(i);
    grad_j = (yj - glm_link_loc(arma::sum(xj % theta))) * xj;
    J = J + grad_j.t() * grad_j / N;
  }
  arma::mat J1 = J.i(); // inverse
  double trJ1 = arma::trace(J1); // trace J
  double trJ2 = arma::trace(J1 * J1); // trace J^-2
  // Laguerre bound
  double a = p * trJ2 / std::pow(trJ1,2) - 1;
  double b = std::sqrt((p-1) * a);
  double lbound = (1 / trJ1) * p / (1 + b);
  double trA = arma::trace(J);
  double ubound = trA / p;
  arma::rowvec output;
  output << lbound << ubound << arma::endr;
  return(output);
}


/* global variables for implicit update */
double gamma_n;
double y_n;
double pred_n;
//double norm_xn;
double C_norm_xn;
double (*glm_link) (double);

// returns scalar value of yn - h(theta_{n-1}' xn + xn^2 ξ)  -- for a GLM
double implicitBound(double ksi) {
	return (gamma_n * (y_n - glm_link(pred_n + C_norm_xn * ksi)));
}

/* for implicit update calculation */
double implicitFn(double u) {
	return (u - implicitBound(u));
}

// try to pass in seed, may make things faster
// [[Rcpp::export]]
arma::rowvec implicit_sgd(arma::mat dataX, arma::vec dataY, 
                          arma::rowvec C,
													std::string model_name, 
													double gamma, int n_pass, 
													uint32_t seed,
													bool use_permutation) {
  // dimension parameters
  int N = dataX.n_rows;
  int p = dataX.n_cols;
  int niters = N * n_pass;

	if (C.has_nan()) {
		C = arma::ones<arma::rowvec>(p); 
	}	
		
  // chose glm_link
  if (model_name.compare("gaussian") == 0) {
    glm_link = identity;
  } else if (model_name.compare("binomial") == 0) {
    glm_link = expit;
  } else if (model_name.compare("poisson") == 0) {
		glm_link = std::exp; 
	} else {
    throw std::invalid_argument("model_name not gaussian or binomial or poisson");
  }
  
  // Random generator
  boost::random::mt19937 gen(seed);
  boost::random::uniform_int_distribution<> dist(0, N-1);
  int i;
  
  // udpate variables
  arma::rowvec x_n(p);
  arma::rowvec theta = arma::zeros<arma::rowvec>(p);
  double bound_left;
  double bound_right;
  double ksi_star;
  for(int n=1; n < niters; n++) {
    i = n-1; 
    if(use_permutation) {
      i = dist(gen);
    }
    x_n = dataX.row(i);
    y_n = dataY[i];
    pred_n = arma::dot(x_n,theta);
    // <------  TODO: IMPORTANT
    //norm_xn = arma::dot(x_n,x_n);
		C_norm_xn = arma::sum(x_n % x_n % C);
		gamma_n = gamma / n;
    
    /* search interval */
    bound_left  = std::min(implicitBound(0), 0.0);
    bound_right = std::max(implicitBound(0), 0.0);
    
    /* solve implicit equation */
    ksi_star = 0;
    if (bound_left != bound_right) {
      ksi_star = zeroin(bound_left,bound_right,implicitFn,0.0,1000);
    }	
    
    if(!arma::is_finite(ksi_star)) {
      printf("->LB+%.1f - RB=%.1f\n %.1f\t%.1f\n", 
             bound_left, bound_right, 
             pred_n, C_norm_xn);
      
      n = niters;
    }
    /* Main update */
    theta = theta + ksi_star * (x_n % C);
  }
  
  return theta;
}	

// [[Rcpp::export]]
arma::rowvec implicit_sgd_const(arma::mat dataX, arma::vec dataY, 
                          arma::rowvec C,
                          std::string model_name, 
                          double gamma, int n_pass, 
                          uint32_t seed,
                          bool use_permutation) {
  // dimension parameters
  int N = dataX.n_rows;
  int p = dataX.n_cols;
  int niters = N * n_pass;
  
  if (C.has_nan()) {
    C = arma::ones<arma::rowvec>(p); 
  }	
  
  // chose glm_link
  if (model_name.compare("gaussian") == 0) {
    glm_link = identity;
  } else if (model_name.compare("binomial") == 0) {
    glm_link = expit;
  } else if (model_name.compare("poisson") == 0) {
    glm_link = std::exp; 
  } else {
    throw std::invalid_argument("model_name not gaussian or binomial or poisson");
  }
  
  // Random generator
  boost::random::mt19937 gen(seed);
  boost::random::uniform_int_distribution<> dist(0, N-1);
  int i;
  
  // udpate variables
  arma::rowvec x_n(p);
  arma::rowvec theta = arma::zeros<arma::rowvec>(p);
  double bound_left;
  double bound_right;
  double ksi_star;
  for(int n=1; n < niters; n++) {
    i = n-1; 
    if(use_permutation) {
      i = dist(gen);
    }
    x_n = dataX.row(i);
    y_n = dataY[i];
    pred_n = arma::dot(x_n,theta);
    //norm_xn = arma::dot(x_n,x_n);
    C_norm_xn = arma::sum(x_n % x_n % C);
    gamma_n = gamma;
    
    /* search interval */
    bound_left  = std::min(implicitBound(0), 0.0);
    bound_right = std::max(implicitBound(0), 0.0);
    
    /* solve implicit equation */
    ksi_star = 0;
    if (bound_left != bound_right) {
      ksi_star = zeroin(bound_left,bound_right,implicitFn,0.0,1000);
    }	
    
    if(!arma::is_finite(ksi_star)) {
      printf("->LB+%.1f - RB=%.1f\n %.1f\t%.1f\n", 
             bound_left, bound_right, 
             pred_n, C_norm_xn);
      
      n = niters;
    }
    /* Main update */
    theta = theta + ksi_star * (x_n % C);
  }
  
  return theta;
}	

// [[Rcpp::export]]
arma::rowvec implicit_sgd_0(arma::rowvec theta0,
                            arma::mat dataX, arma::vec dataY, 
                            arma::rowvec C,
                            std::string model_name, 
                            double gamma, int n_pass, 
                            uint32_t seed,
                            bool use_permutation) {
  // dimension parameters
  int N = dataX.n_rows;
  int p = dataX.n_cols;
  int niters = N * n_pass;
  
  if (C.has_nan()) {
    C = arma::ones<arma::rowvec>(p); 
  }	
  
  // chose glm_link
  if (model_name.compare("gaussian") == 0) {
    glm_link = identity;
  } else if (model_name.compare("binomial") == 0) {
    glm_link = expit;
  } else if (model_name.compare("poisson") == 0) {
    glm_link = std::exp; 
  } else {
    throw std::invalid_argument("model_name not gaussian or binomial or poisson");
  }
  
  // Random generator
  boost::random::mt19937 gen(seed);
  boost::random::uniform_int_distribution<> dist(0, N-1);
  int i;
  
  // udpate variables
  arma::rowvec x_n(p);
  arma::rowvec theta = theta0;
  double bound_left;
  double bound_right;
  double ksi_star;
  for(int n=1; n < niters; n++) {
    i = n-1;
    if(use_permutation) {
      i = dist(gen);
    }
    
    x_n = dataX.row(i);
    y_n = dataY[i];
    pred_n = arma::dot(x_n,theta);
    //norm_xn = arma::dot(x_n,x_n);
    C_norm_xn = arma::sum(x_n % x_n % C);
    gamma_n = gamma / n;
    
    /* search interval */
    bound_left  = std::min(implicitBound(0), 0.0);
    bound_right = std::max(implicitBound(0), 0.0);
    
    /* solve implicit equation */
    ksi_star = 0;
    if (bound_left != bound_right) {
      ksi_star = zeroin(bound_left,bound_right,implicitFn,0.0,1000);
    }	
    
    if(!arma::is_finite(ksi_star)) {
      printf("->LB+%.1f - RB=%.1f\n %.1f\t%.1f\n", 
             bound_left, bound_right, 
             pred_n, C_norm_xn);
      
      n = niters;
    }
    /* Main update */
    theta = theta + ksi_star * (x_n % C);
  }
  
  return theta;
}	

// [[Rcpp::export]]
arma::rowvec implicit_sgd_avg(arma::mat dataX, arma::vec dataY, 
                          arma::rowvec C,
                          std::string model_name, 
                          double gamma, int n_pass, 
                          uint32_t seed,
                          bool use_permutation) {
  // dimension parameters
  int N = dataX.n_rows;
  int p = dataX.n_cols;
  int niters = N * n_pass;
  
  if (C.has_nan()) {
    C = arma::ones<arma::rowvec>(p); 
  }	
  
  // chose glm_link
  if (model_name.compare("gaussian") == 0) {
    glm_link = identity;
  } else if (model_name.compare("binomial") == 0) {
    glm_link = expit;
  } else if (model_name.compare("poisson") == 0) {
    glm_link = std::exp; 
  } else {
    throw std::invalid_argument("model_name not gaussian or binomial or poisson");
  }
  
  // Random generator
  boost::random::mt19937 gen(seed);
  boost::random::uniform_int_distribution<> dist(0, N-1);
  int i;
  
  // udpate variables
  arma::rowvec x_n(p);
  arma::rowvec theta = arma::zeros<arma::rowvec>(p);
  arma::rowvec theta_avg = arma::zeros<arma::rowvec>(p);
  double bound_left;
  double bound_right;
  double ksi_star;
  for(double n=1; n < niters; n++) {
    i = n-1;
    if(use_permutation) {
      i = dist(gen);
    }
    x_n = dataX.row(i);
    y_n = dataY[i];
    pred_n = arma::dot(x_n,theta);
    //norm_xn = arma::dot(x_n,x_n);
    C_norm_xn = arma::sum(x_n % x_n % C);
    gamma_n = gamma / sqrt(n);
    
    /* search interval */
    bound_left  = std::min(implicitBound(0), 0.0);
    bound_right = std::max(implicitBound(0), 0.0);
    
    /* solve implicit equation */
    ksi_star = 0;
    if (bound_left != bound_right) {
      ksi_star = zeroin(bound_left,bound_right,implicitFn,0.0,1000);
    }	
    
    if(!arma::is_finite(ksi_star)) {
      printf("->LB+%.1f - RB=%.1f\n %.1f\t%.1f\n", 
             bound_left, bound_right, 
             pred_n, C_norm_xn);
      
      n = niters;
    }
    /* Main update */
    theta = theta + ksi_star * (x_n % C);
    theta_avg = (1/n) * ((n-1) * theta_avg + theta);
  }
  
  return theta_avg;
}	

// [[Rcpp::export]]
arma::rowvec implicit_sgd_avg_0(arma::rowvec theta0,
                              arma::mat dataX, arma::vec dataY, 
                              arma::rowvec C,
                              std::string model_name, 
                              double gamma, int n_pass, 
                              uint32_t seed,
                              bool use_permutation) {
  // dimension parameters
  int N = dataX.n_rows;
  int p = dataX.n_cols;
  int niters = N * n_pass;
  
  if (C.has_nan()) {
    C = arma::ones<arma::rowvec>(p); 
  }	
  
  // chose glm_link
  if (model_name.compare("gaussian") == 0) {
    glm_link = identity;
  } else if (model_name.compare("binomial") == 0) {
    glm_link = expit;
  } else if (model_name.compare("poisson") == 0) {
    glm_link = std::exp; 
  } else {
    throw std::invalid_argument("model_name not gaussian or binomial or poisson");
  }
  
  // Random generator
  boost::random::mt19937 gen(seed);
  boost::random::uniform_int_distribution<> dist(0, N-1);
  int i;
  
  // udpate variables
  arma::rowvec x_n(p);
  arma::rowvec theta = theta0;
  arma::rowvec theta_avg = theta0;
  double bound_left;
  double bound_right;
  double ksi_star;
  for(double n=1; n < niters; n++) {
    i = n-1;
    if(use_permutation) {
      i = dist(gen);
    }
    x_n = dataX.row(i);
    y_n = dataY[i];
    pred_n = arma::dot(x_n,theta);
    //norm_xn = arma::dot(x_n,x_n);
    C_norm_xn = arma::sum(x_n % x_n % C);
    gamma_n = gamma / sqrt(n);
    
    /* search interval */
    bound_left  = std::min(implicitBound(0), 0.0);
    bound_right = std::max(implicitBound(0), 0.0);
    
    /* solve implicit equation */
    ksi_star = 0;
    if (bound_left != bound_right) {
      ksi_star = zeroin(bound_left,bound_right,implicitFn,0.0,1000);
    }	
    
    if(!arma::is_finite(ksi_star)) {
      printf("->LB+%.1f - RB=%.1f\n %.1f\t%.1f\n", 
             bound_left, bound_right, 
             pred_n, C_norm_xn);
      
      n = niters;
    }
    /* Main update */
    theta = theta + ksi_star * (x_n % C);
    theta_avg = (1/n) * ((n-1) * theta_avg + theta);
  }
  
  return theta_avg;
}	

/* ==========================================
 * Unit tests to ensure that zeroin() is working
 */
int counter;
/* test functions */
/* f1 with global variables */
double y = 2.0;
double z = 5.0;
double f1(double x) {
  counter++;
  return (pow(x,y)-y)*x - z;
}
double f2(double x) {
  counter++;
  return cos(x) - x;
}
double f3(double x) {
  counter++;
  return sin(x) - x;
}
double f4(double x) {
  counter++;
  return std::log(x);
}
double f5(double x) {
  counter++;
  return x * std::sin(x);
}
double f6(double x) {
  counter++;
  return std::pow(x,3.0) * std::cos(x) * std::sin(x);
}
void test_zeroin(
    double a, double b, /* range root seeked for */
    double (*f) (double x), /* function */
    const char * msg, /* symbolic function */
    double truth /* wolfram calculated root */
    ) {
  counter = 0;
  double root = zeroin(a,b,f,0.0,1000);
  double diff = std::abs(root - truth);

  /* check that root calc good enough */
  if (diff > 1e-5) {
		printf("Function\t%s",msg);
    throw std::runtime_error("Difference between uniroot() and  wolfram");
  }
  /* make sure that function at root is zero */
  if (std::abs((*f)(root) - 0.0) > 1e-5) {
		printf("Function\t%s",msg);
    throw std::runtime_error("f(root) not at zero");
  }
}

/* ensures that zeroin() working */
// [[Rcpp::export]]
void test_uniroot(){
  Rcout << "machine precision<double>: "
    << std::numeric_limits<double>::epsilon() << std::endl;
  Rcout << "sqrt(machine precision<double>): "
    << std::sqrt( std::numeric_limits<double>::epsilon() ) << std::endl;

  test_zeroin(2.0,3.0,f1,"x^3 - 2*x - 5",2.0945514815);
  test_zeroin(-1.0,3.0,f2,"cos(x)-x",0.739085);
  test_zeroin(-1.0,3.0,f3,"sin(x)-x",0.0);
  test_zeroin(0.5,1.5,f4,"log(x)",1.0);
  test_zeroin(2.0,4.0,f5,"xsin(x)",3.14159265359);
  test_zeroin(-2.0,-1.0,f6,"x^3sin(x)cos(x)",-1.57079632679);

  printf("\nAll tests GOOD for zeroin()!\n");
}

/* ============================================
 * Unit tests to ensure implicit_sgd() working
 */
// [[Rcpp::export]]
void test_isgd(int N, int p) {
	arma::vec theta_star(p); 
	for (int i=0; i<p; i++) {
		theta_star[i] = std::exp(-0.75 * i);
	}
	arma::mat X;
	X.load("X_test.csv",arma::csv_ascii);
	arma::vec eps(N,arma::fill::randn); 
	arma::vec Y = X * theta_star + 0.5 * eps; 
	arma::rowvec C = arma::ones<arma::rowvec>(p);
	arma::rowvec theta_est = implicit_sgd(X,Y,C,
                                       "gaussian",
                                       2.0,1,
                                       std::time(0),
                                       false);
	arma::vec diff = theta_est.t() - theta_star;
	if ( std::abs( arma::dot(diff,diff) / p ) > 1e-2 ) {
		Rcout << arma::dot(diff,diff) << std::endl; 
		throw std::runtime_error("implicit_sgd() does not give est close enough");
	}
	printf("\nTest GOOD for implicit_sgd!\n");
}


/*** R
*/

// R code run with source(echo = TRUE) so don't need explicitly print output

/* JERRY CHEE made changes to max_iter, define epsilon
 ************************************************************************
 *	    		    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,tol)
 *	double ax; 			Root will be seeked for within
 *	double bx;  			a range [ax,bx]
 *	double (*f)(double x);		Name of the function whose zero
 *					will be seeked for
 *	double tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bissection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *		   the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bissection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is 
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bissection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */

double zeroin(		/* An estimate to the root	*/
	double ax,				/* Left border | of the range	*/
	double bx,  				/* Right border| the root is seeked*/
	double (*f)(double x),			/* Function under investigation	*/
	double tol,				/* Acceptable tolerance		*/
	int max_iter)			/* maximum number of iterations */
{
	double a,b,c;				/* Abscissae, descr. see above	*/
	double fa;				/* f(a)				*/
	double fb;				/* f(b)				*/
	double fc;				/* f(c)				*/

	a = ax;  b = bx;  fa = (*f)(a);  fb = (*f)(b);
	c = a;   fc = fa;
	double EPSILON = std::sqrt( std::numeric_limits<double>::epsilon() );

	for(int i=0; i<max_iter; i++)		/* Main iteration loop	*/
	{
		double prev_step = b-a;		/* Distance from the last but one*/
		/* to the last approximation	*/
		double tol_act;			/* Actual tolerance		*/
		double p;      			/* Interpolation step is calcu- */
		double q;      			/* lated in the form p/q; divi- */
		/* sion operations is delayed   */
		/* until the last moment	*/
		double new_step;      		/* Step at this iteration       */

		if( std::abs(fc) < std::abs(fb) )
		{                         		/* Swap data for b to be the 	*/
			a = b;  b = c;  c = a;          /* best approximation		*/
			fa=fb;  fb=fc;  fc=fa;
		}
		tol_act = 2*EPSILON*std::abs(b) + tol/2;

		new_step = (c-b)/2;

		if( std::abs(new_step) <= tol_act || fb == (double)0 )
		{
			return b;				/* Acceptable approx. is found	*/
		}

		/* Decide if the interpolation can be tried	*/
		if( std::abs(prev_step) >= tol_act	/* If prev_step was large enough*/
				&& std::abs(fa) > std::abs(fb) )	/* and was in true direction,	*/
		{					/* Interpolatiom may be tried	*/
			register double t1,cb,t2;
			cb = c-b;
			if( a==c )			/* If we have only two distinct	*/
			{				/* points linear interpolation 	*/
				t1 = fb/fa;			/* can only be applied		*/
				p = cb*t1;
				q = 1.0 - t1;
			}
			else				/* Quadric inverse interpolation*/
			{
				q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
				p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);
			}
			if( p>(double)0 )		/* p was calculated with the op-*/
			{
				q = -q;			/* posite sign; make p positive	*/
			}
			else				/* and assign possible minus to	*/
			{
				p = -p;			/* q				*/
			}
			if( p < (0.75*cb*q-std::abs(tol_act*q)/2)	/* If b+p/q falls in [b,c]*/
					&& p < std::abs(prev_step*q/2) )	/* and isn't too large	*/
			{
				new_step = p/q;			/* it is accepted	*/
			}
			/* If p/q is too large then the	*/
			/* bissection procedure can 	*/
			/* reduce [b,c] range to more	*/
			/* extent			*/
		}

		if( std::abs(new_step) < tol_act )	/* Adjust the step to be not less*/
		{
			if( new_step > (double)0 )	/* than tolerance		*/
			{
				new_step = tol_act;
			}
			else
			{
				new_step = -tol_act;
			}
		}

		a = b;  fa = fb;			/* Save the previous approx.	*/
		b += new_step;  fb = (*f)(b);	/* Do step to a new approxim.	*/
		if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
		{                 			/* Adjust c for it to have a sign*/
			c = a;  fc = fa;                  /* opposite to that of b	*/
		}
	}
	return b;

}
