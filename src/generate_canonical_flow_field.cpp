#include "common.h"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include <boost/timer/progress_display.hpp>

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include "boost/asio/thread_pool.hpp"
#include "boost/asio.hpp"

#include <Eigen/Dense>
using namespace Eigen;

#include "arb.h"
#include "arb_hypgeom.h"

#include <gsl/gsl_integration.h>

#include <thread>


mutex progress_bar_lock;
boost::timer::progress_display* progress_bar;


// Integral of Gamma(a, b) from t to infty is Upper_Incomplete_Gamma_Function(a, b * t) / Gamma(a)
// gamma_q_inv(a, q) = x ---> q = gamma_q(a, x) = Upper_Incomplete_Gamma_Function(a, x) / Gamma(a)

// Therefore, the inverse survival function at x is:
// q = Upper_Incomplete_Gamma_Function(a, b * t) / Gamma(a) ---> 
// gamma_isf = gamma_q_inv(a, q) / b 
//
double gamma_isf(
    double alpha,
    double beta,
    double q) {
    return boost::math::gamma_q_inv(alpha, q) / beta; 
}

// Exp CDF is 1-exp(-x). SF is 1-CDF = exp(-x). 
// isf is therefore: q = exp(-x) ---> -log(q) = x
double exp_isf(double q) {
    return -log(q);
}

void get_times_grid(
    double alpha, 
    double beta, 
    double integration_min_time,
    double integration_sf_at_max_time_gamma, 
    int integration_n_steps,
    double integration_sf_at_max_time_exp,
    int integration_n_steps_exp,
    vector<double>& times_grid) 
{

    double x = integration_min_time;
    double first_step_size = (gamma_isf(alpha, beta, integration_sf_at_max_time_gamma) - integration_min_time) / integration_n_steps;
    for (int i = 0; i < integration_n_steps + 1; i++) {
        times_grid.push_back(x);
        x += first_step_size;
    }

    double second_step_size = (exp_isf(integration_sf_at_max_time_exp) / (integration_n_steps_exp - 1));
    for (int i = 0; i < integration_n_steps_exp; i++) {
        x += second_step_size;
        times_grid.push_back(x);   // Order is reversed to not repeat the last element of the first grid
    }
}


double distribution_difference_pdf(
    double t,
    double alpha, 
    double beta) {

    slong prec;
    arb_t x, a, b, z, lx;
    arb_init(x);
    arb_init(a);
    arb_init(b);
    arb_init(z);
    arb_init(lx);

    /*
    first_part = np.exp(
        (-t) + alpha * np.log(t * beta) - scipy.special.loggamma(alpha+1) + \
            float(flint.good(lambda: flint.arb(-(beta-1)*t).hypgeom_1f1(alpha, alpha+1, regularized=False), maxprec=1000000).log())
    )
    */
    arb_set_d(a, alpha);
    arb_set_d(b, alpha+1);
    arb_set_d(z, -(beta-1)*t);
 
    prec = 20;
    while (prec <= 1000000) {        
        arb_hypgeom_1f1(x, a, b, z, false, prec * 1.01 + 2 * 10);
        if (arb_rel_accuracy_bits(x) > prec) {
            break;
        }
        prec *= 2;
    }
    arb_log(lx, x, prec);

    double first_part = std::exp(
        (-t) + alpha * std::log(t * beta) - boost::math::lgamma(alpha+1) + arf_get_d(arb_midref(lx), ARF_RND_NEAR)
    );
    // arb_printn(x, 20, 0); cout << endl;
    // arb_printn(lx, 20, 0); cout << endl;
    // cout << boost::format("first_part = %f\n") % arf_get_d(arb_midref(lx), ARF_RND_NEAR);

    /*
    second_part = np.exp(
        (-t) + alpha * np.log(t * beta) - scipy.special.loggamma(alpha+1) + \
            float(flint.good(lambda: flint.arb(-(beta+1)*t).hypgeom_1f1(alpha, alpha+1, regularized=False), maxprec=1000000).log())
    )
    */
    arb_set_d(a, alpha);
    arb_set_d(b, alpha+1);
    arb_set_d(z, -(beta+1)*t);
 
    prec = 20;
    while (prec <= 1000000) {        
        arb_hypgeom_1f1(x, a, b, z, false, prec * 1.01 + 2 * 10);
        if (arb_rel_accuracy_bits(x) > prec) {
            break;
        }
        prec *= 2;
    }
    arb_log(lx, x, prec);

    double second_part = std::exp(
        (-t) + alpha * std::log(t * beta) - boost::math::lgamma(alpha+1) + arf_get_d(arb_midref(lx), ARF_RND_NEAR)
    );
    // cout << boost::format("second_part = %f\n") % second_part;

    double res = first_part - second_part;

    // res += (np.exp(-2*t)/2 - 0.5 - t) * np.exp(log_gamma_pdf)
    double log_gamma_pdf = alpha * std::log(beta) - boost::math::lgamma(alpha) + (alpha - 1) * std::log(t) - beta * t;
    res += (std::exp(-2*t)/2 - 0.5 - t) * std::exp(log_gamma_pdf);

    // cout << boost::format("res2 = %f\n") % res;


    // res += (1 - np.exp(-2*t)) * scipy.special.gammaincc(alpha, beta*t) # includes factor of scipy.special.gamma(alpha)            
    res += (1 - std::exp(-2*t)) * boost::math::gamma_q(alpha, beta*t);

    // cout << boost::format("res3 = %f\n") % res;

    // exit(-1);



    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
    arb_clear(z);
    arb_clear(lx);

    return res;

}



double distribution_difference_pdf_arb_integration(
    double t,
    double alpha, 
    double beta) {

    slong prec;
    arb_t x, a, b, z, lx;
    arb_init(x);
    arb_init(a);
    arb_init(b);
    arb_init(z);
    arb_init(lx);

    /*
    first_part = np.exp(
        (-t) + alpha * np.log(t * beta) - scipy.special.loggamma(alpha+1) + \
            float(flint.good(lambda: flint.arb(-(beta-1)*t).hypgeom_1f1(alpha, alpha+1, regularized=False), maxprec=1000000).log())
    )
    */
    arb_set_d(a, alpha);
    arb_set_d(b, alpha+1);
    arb_set_d(z, -(beta-1)*t);
 
    prec = 20;
    while (prec <= 1000000) {        
        arb_hypgeom_1f1_integration(x, a, b, z, false, prec * 1.01 + 2 * 10);
        if (arb_rel_accuracy_bits(x) > prec) {
            break;
        }
        prec *= 2;
    }
    arb_log(lx, x, prec);

    double first_part = std::exp(
        (-t) + alpha * std::log(t * beta) - boost::math::lgamma(alpha+1) + arf_get_d(arb_midref(lx), ARF_RND_NEAR)
    );
    // arb_printn(x, 20, 0); cout << endl;
    // arb_printn(lx, 20, 0); cout << endl;
    // cout << boost::format("first_part = %f\n") % arf_get_d(arb_midref(lx), ARF_RND_NEAR);

    /*
    second_part = np.exp(
        (-t) + alpha * np.log(t * beta) - scipy.special.loggamma(alpha+1) + \
            float(flint.good(lambda: flint.arb(-(beta+1)*t).hypgeom_1f1(alpha, alpha+1, regularized=False), maxprec=1000000).log())
    )
    */
    arb_set_d(a, alpha);
    arb_set_d(b, alpha+1);
    arb_set_d(z, -(beta+1)*t);
 
    prec = 20;
    while (prec <= 1000000) {        
        arb_hypgeom_1f1_integration(x, a, b, z, false, prec * 1.01 + 2 * 10);
        if (arb_rel_accuracy_bits(x) > prec) {
            break;
        }
        prec *= 2;
    }
    arb_log(lx, x, prec);

    double second_part = std::exp(
        (-t) + alpha * std::log(t * beta) - boost::math::lgamma(alpha+1) + arf_get_d(arb_midref(lx), ARF_RND_NEAR)
    );
    // cout << boost::format("second_part = %f\n") % second_part;

    double res = first_part - second_part;

    // res += (np.exp(-2*t)/2 - 0.5 - t) * np.exp(log_gamma_pdf)
    double log_gamma_pdf = alpha * std::log(beta) - boost::math::lgamma(alpha) + (alpha - 1) * std::log(t) - beta * t;
    res += (std::exp(-2*t)/2 - 0.5 - t) * std::exp(log_gamma_pdf);

    // cout << boost::format("res2 = %f\n") % res;


    // res += (1 - np.exp(-2*t)) * scipy.special.gammaincc(alpha, beta*t) # includes factor of scipy.special.gamma(alpha)            
    res += (1 - std::exp(-2*t)) * boost::math::gamma_q(alpha, beta*t);

    // cout << boost::format("res3 = %f\n") % res;

    // exit(-1);



    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
    arb_clear(z);
    arb_clear(lx);

    return res;

}

struct integral_f_params { double alpha; double beta; double t; };

double integral_f(double s, void *p) {
    struct integral_f_params* params = (struct integral_f_params *)p;
    double alpha = (params->alpha);
    double beta = (params->beta);
    double t = (params->t);

    double ret = (std::exp(s-t) - std::exp(-s-t)) * 
        std::exp(alpha * std::log(beta) - boost::math::lgamma(alpha) + (alpha - 1) * std::log(s) - beta * s);

    // if ((alpha >= 10000) && (beta >= 100000)) {
    //     cout << boost::format("alpha, beta, t -> ret: %f, %f, %f -> %1.20f\n")
    //         % alpha % beta % t % ret;
    // }

    return ret;
}

double distribution_difference_pdf_gsl(
    double t,
    double alpha, 
    double beta) {

    double res = 0.0;

    // Calculate the integral by numerical integration
    gsl_function F;
    struct integral_f_params params = { alpha, beta, t };

    F.function = &integral_f;
    F.params = &params;

    int max_subintervals = 1000;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(max_subintervals);
    double result, error;

    // double breakpoints[2];
    // breakpoints[0] = 0.0;
    // breakpoints[1] = t;

    // gsl_integration_qagp(
    //     &F,                 // Function to integrate
    //     breakpoints,        // Known singularity points
    //     2,                  // Number of singularity points     
    //     0,                  // Absolute error. Set to 0 to only activate relative error
    //     1e-7,               // Relative error
    //     max_subintervals,   // Maximum subintervals
    //     w,                  // Workspace 
    //     &result, &error     // Return values
    // );

    gsl_integration_qag(
        &F,                 // Function to integrate
        1e-20,              // Lower bound
        t,                  // Upper bound
        0,                  // Absolute error. Set to 0 to only activate relative error
        1e-7,               // Relative error
        max_subintervals,   // Maximum subintervals
        GSL_INTEG_GAUSS61,  // Integration rule
        w,                  // Workspace 
        &result, &error     // Return values
    );

    res = result;

    double log_gamma_pdf = alpha * std::log(beta) - boost::math::lgamma(alpha) + (alpha - 1) * std::log(t) - beta * t;
    res += (std::exp(-2*t)/2 - 0.5 - t) * std::exp(log_gamma_pdf);

    res += (1 - std::exp(-2*t)) * boost::math::gamma_q(alpha, beta*t);

    gsl_integration_workspace_free(w);


    return res;

}


void calculate_updated_gamma_parameters(
    double mean,
    double cv,
    double integration_min_time,
    double integration_sf_at_max_time,
    int integration_n_steps,
    bool log_coords,
    bool use_numeric_integration,
    double* u,
    double* v
    ) {
    
    double std = cv * mean;
    double alpha = (mean * mean) / (std * std);
    double beta = mean / (std * std);

    vector<double> times_grid;
    get_times_grid(
        alpha, 
        beta, 
        integration_min_time,
        integration_sf_at_max_time, 
        integration_n_steps,
        1e-3,
        integration_n_steps,
        times_grid);

    uint n = times_grid.size();

    VectorXd pdf_grid(n);
    MatrixXd partials_grid(n, 2);

    for (uint i = 0; i < n; i++) {
        double t = times_grid[i];
        double time_diff = 0.0;
        if (i < n-1) {
            time_diff = std::sqrt(times_grid[i+1] - times_grid[i]);
        }

        if (use_numeric_integration) {
            // pdf_grid(i) = distribution_difference_pdf_gsl(t, alpha, beta) * time_diff;
            pdf_grid(i) = distribution_difference_pdf_arb_integration(t, alpha, beta) * time_diff;            
        } else {
            pdf_grid(i) = distribution_difference_pdf(t, alpha, beta) * time_diff;
        }

        double gamma_pdf = std::exp(alpha * std::log(beta) - boost::math::lgamma(alpha) + (alpha - 1) * std::log(t) - beta * t);
        partials_grid(i, 0) = gamma_pdf * (-boost::math::digamma(alpha) + std::log(beta) + std::log(t)) * time_diff;
        partials_grid(i, 1) = gamma_pdf * (alpha/beta - t) * time_diff;

        if (log_coords) {
            partials_grid(i, 0) *= alpha / std::log10(std::exp(1));
            partials_grid(i, 1) *= beta / std::log10(std::exp(1));
        }
    }

    VectorXd uv(partials_grid.bdcSvd(ComputeThinU | ComputeThinV).solve(pdf_grid));

    
    // *u = uv(0);
    // *v = uv(1);   
    *u = uv(0) - uv(1);
    *v = -0.5*uv(0);


    progress_bar_lock.lock();
    (*progress_bar) += 1;
    progress_bar_lock.unlock();
}

void main_processing(const po::variables_map& vm, int n_threads) {
    double log10mean, mean, log10cv, cv; 
    
    int mean_n_steps = vm["mean_n_steps"].as<int>();
    double mean_step = (std::log10(vm["mean_max"].as<double>()) -  std::log10(vm["mean_min"].as<double>())) / (mean_n_steps - 1);    
    int cv_n_steps = vm["cv_n_steps"].as<int>();
    double cv_step = (std::log10(vm["cv_max"].as<double>()) -  std::log10(vm["cv_min"].as<double>())) / (cv_n_steps - 1);    

    vector<double> results_u(mean_n_steps * cv_n_steps);
    vector<double> results_v(mean_n_steps * cv_n_steps);

    boost::asio::thread_pool thread_pool(n_threads);

    progress_bar = new boost::timer::progress_display(mean_n_steps * cv_n_steps);
    
    double* ptr_u = results_u.data();
    double* ptr_v = results_v.data();
    log10mean = std::log10(vm["mean_min"].as<double>());
    for (int i = 0; i < mean_n_steps; i++) {
        mean = std::pow(10, log10mean);
        log10cv = std::log10(vm["cv_min"].as<double>());
        for (int j = 0; j < cv_n_steps; j++) {
            cv = std::pow(10, log10cv);
            
            //cout << boost::format("%d %d - %f %f\n") % i % j % mean % cv;

            
            boost::asio::post(thread_pool, boost::bind(
                calculate_updated_gamma_parameters, 
                mean,
                cv,
                vm["integration_min_time"].as<double>(),
                vm["integration_sf_at_max_time"].as<double>(),
                vm["integration_n_steps"].as<int>(),
                vm["log_coords"].as<bool>(),
                vm["use_numeric_integration"].as<bool>(),
                ptr_u,
                ptr_v
            ));
            

            
            
            // calculate_updated_gamma_parameters(
            //     mean,
            //     cv,
            //     vm["integration_min_time"].as<double>(),
            //     vm["integration_sf_at_max_time"].as<double>(),
            //     vm["integration_n_steps"].as<int>(),
            //     vm["log_coords"].as<bool>(),
            //     vm["use_numeric_integration"].as<bool>(),
            //     ptr_u,
            //     ptr_v
            // );
            
            

            log10cv += cv_step;
            ptr_u++;
            ptr_v++;
        }
        log10mean += mean_step;
    }

    thread_pool.join();   

    ofstream output_file;
    output_file.open(vm["output"].as<string>()); 
    output_file << boost::format("%f %f %d\n") % vm["mean_min"].as<double>() % vm["mean_max"].as<double>() % vm["mean_n_steps"].as<int>();
    output_file << boost::format("%f %f %d\n") % vm["cv_min"].as<double>() % vm["cv_max"].as<double>() % vm["cv_n_steps"].as<int>();

    int ptr = 0;
    for (int i = 0; i < mean_n_steps; i++) {
        for (int j = 0; j < cv_n_steps; j++) {
            output_file << results_u[ptr] << " ";
            ptr++;
        }
        output_file << "\n";
    }

    ptr = 0;
    for (int i = 0; i < mean_n_steps; i++) {
        for (int j = 0; j < cv_n_steps; j++) {
            output_file << results_v[ptr] << " ";
            ptr++;
        }
        output_file << "\n";
    }

    output_file.close();
}

int main(int argc, char** argv) {
    // Print command line
    cout << "Command line:" << endl;
    for (int i = 0; i < argc; ++i) {
        cout << argv[i] << ' ';
    }
    cout << endl << endl;

    //
    // Parse flags
    //
    po::options_description desc("Allowed options");
    desc.add_options()
        ("mean_min", po::value<double>()->default_value(1e-2), "The minimum mean to evaluate")
        ("mean_max", po::value<double>()->default_value(1e2), "The maximum mean to evaluate")
        ("mean_n_steps", po::value<int>()->default_value(5), "The number of grid steps for mean")
        ("cv_min", po::value<double>()->default_value(1e-2), "The minimum coefficient of variation to evaluate")
        ("cv_max", po::value<double>()->default_value(1), "The maximum coefficient of variation to evaluate")
        ("cv_n_steps", po::value<int>()->default_value(3), "The number of grid steps for coefficient of variation")
        ("integration_min_time", po::value<double>()->default_value(1e-10), "Defines the lower bound for integration")
        ("integration_sf_at_max_time", po::value<double>()->default_value(1e-3), "Defines the upper bound of integration when this survival function (1-cdf) is obtained")
        ("integration_n_steps", po::value<int>()->default_value(1000), "Number of points to use for integration")
        ("n_threads", po::value<int>()->default_value(-1), "Number of threads")
        ("output,o", po::value<string>(), "Output filename path")
        ("log_coords,l", po::value<bool>()->default_value(true), "Use log coordinates in output?")
        ("use_numeric_integration,n", po::value<bool>()->default_value(false), "Use numeric integration")        
    ;


    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
          options(desc).run(), vm);
    po::notify(vm);

    // TODO: Validate flags thoroughly
    // Figure out number of threads
	int n_threads = 1;
	int concurrentThreadsSupported = std::thread::hardware_concurrency();
	if (concurrentThreadsSupported > 0) {
		n_threads = concurrentThreadsSupported ;
	}
	if (vm["n_threads"].as<int>() > 0) {
		n_threads = vm["n_threads"].as<int>();
	}	
    cout << boost::format("Using %d threads...\n") % n_threads;

    auto t1 = std::chrono::high_resolution_clock::now();

    main_processing(vm, n_threads);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    cout << boost::format("Done (%f sec)") % (ms_double.count()/1000) << endl;

    return 0;
}
