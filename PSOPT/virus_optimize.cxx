#include "psopt.h"

const double epsilon_0 = 8.854e-12;


#define u0_state 0
#define u1_state 1

#define x0_v_state 2
#define x1_v_state 3

#define x0_h_state 4
#define x1_h_state 5

double T0 = 1e-8;
double U0 = 1.0;
double X0 = 1e-6;

Workspace * workspace;

std::vector<double> u0_store;
std::vector<double> u1_store;
std::vector<double> u2_store;
std::vector<double> hit_store;

struct Cell{
    double extracellular_conductivity; // S/m
    double extracellular_permittivity; // relative
    double intracellular_conductivity; // S/m
    double intracellular_permittivity; // relative
    double membrane_conductivity; // S/m
    double membrane_permittivity; // relative
    double cell_diameter;  // meters
    double membrane_thickness;


    double a1 = 0;
    double a2 = 0;
    double a3 = 0;
    double b1 = 0;
    double b2 = 0;
    double b3 = 0;
    double R = 0;

    double alpha = 0;
    double beta = 0;
    double gamma = 0;
    double phi = 0;
    double xi = 0;

    void init();

};

void Cell::init(){
    double e_o = extracellular_permittivity * epsilon_0;
    double e_i = intracellular_permittivity * epsilon_0;
    double e_m = membrane_permittivity * epsilon_0;
    R = cell_diameter / 2.0;
    // double selfR = R;

    double l_o = extracellular_conductivity; // S/m
    double l_i = intracellular_conductivity; //S/m
    double l_m = membrane_conductivity; //S/m

    double d = membrane_thickness;
    // epsilon_0

    double sub1 = (3.0 * (R*R) - 3.0 * d * R + (d*d));
    double sub2 = (3.0 * d * R - d*d);

    a1 = 3.0 * d * l_o * ((l_i * (sub1)) + l_m*(sub2)); //eq.9a
    a2 = 3.0 * d * ((l_i * e_o + l_o * e_i) * sub1 + (l_m * e_o + l_o * e_m) * sub2);
    a3 = 3.0 * d * e_o * (e_i * (sub1) + e_m * sub2);

    b1 = 2.0 * (R*R*R) * (l_m + 2.0*l_o) * (l_m + 0.5 * l_i) + 2.0 * ((R-d)*(R-d)*(R-d)) * (l_m - l_o) * (l_i - l_m);

    b2 = 2.0 * (R*R*R) * (l_i * (0.5 * e_m + e_o) + l_m * (0.5*e_i + 2.0*e_m + 2*e_o) + l_o * (e_i + 2.0 * e_m)) + (2.0 * ((R-d)*(R-d)*(R-d))
    * (l_i * (e_m - e_o) + l_m * (e_i - 2.0*e_m + e_o) - l_o * (e_i - e_m))); // is this truly a multiply, or a cross?

    b3 = 2.0 * (R*R*R) * (e_m + 2.0*e_o) * (e_m + 0.5 * e_i) + 2.0 * ((R-d)*(R-d)*(R-d)) * (e_m - e_o) * (e_i - e_m);

    alpha = (R*a3/b3);
    beta = (R*a2/b3);
    gamma = (R*a1/b3);
    phi = (b2/b3);
    xi = (b1/b3);

}



Cell* virus;
Cell* host;

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace){
  // No linkages as this is a single phase problem
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase,Workspace* workspace){

   return -(final_states[ x0_v_state ]*final_states[ x0_v_state ]) + (final_states[ x0_h_state ]*final_states[ x0_h_state ]);
   // return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace){

    //is there a fancy sqrt?
    //take out derivative terms
    // return -states[ x0_v_state ]*100 + states[ x0_h_state ] + sqrt(states[u0_state]*states[u0_state]);

    // return -smooth_fabs(states[ x0_v_state ], 1e-9) + smooth_fabs(states[ x0_h_state ], 1e-9); //also seems to work okay
    //
    //return states[ x0_v_state ]- + (states[ x0_h_state ]*states[ x0_h_state ]) + (states[ x1_h_state ]);

    // return sqrt((states[ x0_v_state ] - (1e-8/X0))*(states[ x0_v_state ] - (1e-8/X0))) + sqrt(states[ x1_v_state ]*states[ x1_v_state ]) + sqrt(states[ x1_h_state ]*states[ x1_h_state ]) + sqrt(states[x0_h_state]*states[x0_h_state]);

    // return ((states[ x0_v_state ] - (1e-4/X0))*(states[ x0_v_state ] - (1e-4/X0))) + (states[ x1_v_state ]*states[ x1_v_state ]) + (states[ x1_h_state ]*states[ x1_h_state ]) + (states[x0_h_state]*states[x0_h_state]);


    // return -(states[ x0_v_state ]*states[ x0_v_state ]) ;

    // double rho = virus->R / host->R;
    adouble input = (states[ u0_state ]*states[ u0_state ]) + 1e-6;
    // double rho = 1;
    // return -((rho*rho)*(states[ x0_v_state ]*states[ x0_v_state ])/input) + ((states[ x0_h_state ]*states[ x0_h_state ])/input);

    double rho = 1;
    return -((rho*rho)*(states[ x0_v_state ]*states[ x0_v_state ])) + (((states[ x0_h_state ]*states[ x0_h_state ])));


    // return (states[ x0_h_state ]*states[ x0_h_state ]) / (states[ x0_v_state ]*states[ x0_v_state ]);



    // return (1.0-(states[ x0_v_state ]*states[ x0_v_state ])) + states[ x0_h_state ]*states[ x0_h_state ];
    // return sqrt((states[ x0_v_state ] - (1e-4/X0))*(states[ x0_v_state ] - (1e-4/X0)));


    // return 0;
}



////////////////////////////////////////////////////////////////////////////
/////////////////// Define the integrand of the integral constraint ///////
////////////////////////////////////////////////////////////////////////////
adouble integrand( adouble* states, adouble* controls, adouble* parameters,
                    adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    //integrand from isoperimetric.cxx
    adouble g;
    adouble u0   = states[ u0_state ];
    adouble u1   = states[ u1_state ];
    adouble u2 = controls[ 0 ];
    // g = (u0*u0);
    // return (((U0 / (T0*T0))*host->alpha*u2 + (U0 / T0)*host->beta*u1 + host->gamma*U0*u0));


    return 0;

    // return (U0 / (T0*T0))*host->alpha * u2 + (U0 / T0)* host->beta * u1 + host->gamma * u0;
    //why doesn't this work?
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace){

    // adouble u0 = initial_states[ u0_state ];
    // e[ 0 ] = u0;
    // adouble u1 = initial_states[ u1_state ];
    // e[ 1 ] = u1;

    adouble controls[1];
    get_initial_controls(controls,xad,iphase,workspace);
    adouble q4 = controls[ 0 ];
    e[ 0 ] = q4;

    adouble x0_v = initial_states[ x0_v_state ];
    e[ 1 ] = x0_v;

    adouble x0_h = initial_states[ x0_h_state ];
    e[ 2 ] = x0_h;

    adouble x1_v = initial_states[ x1_v_state ];
    e[ 3 ] = x1_v;

    adouble x1_h = initial_states[ x1_h_state ];
    e[ 4 ] = x1_h;

    // adouble u0_f = final_states[ 3 ];
    // e[ 3 ] = u0_f;
    //Compute the integral to be constrained
    adouble Q;
    Q = integrate(integrand, xad, iphase, workspace);
    e[ 5 ] = Q;
}



// std::cout << "{\n";
// for(int p = 0; p < 100; p++){
//     std::cout << cj[p] << "\n";
// }
// std::cout << "{\n";

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace){


    //https://mathoverflow.net/a/87902/176668
    // adouble u0   = states[ u0_state ];
    // adouble u1   = states[ u1_state ];
    // adouble u2 = controls[ 0 ];
    // derivatives[ u0_state ] = u1;
    // derivatives[ u1_state ] = u2; //multiply by dt?
    //
    adouble t0, tf;
    get_times(&t0, &tf, xad, iphase, workspace);

    adouble* cj = workspace->single_trajectory_tmp;
    // MatrixXd& control_scaling = workspace->problem->phase[0].scale.controls;
    // cj[0] = xad[0]/control_scaling(0); //control trajectory array - only valid on order 1
    get_individual_control_trajectory(cj, 0, 1, xad, workspace); // needs simplifying I think


    adouble dt = (tf-t0) / workspace->problem->phases(1).nodes(0);

    //this only works with local collocation (trapezoidal, H-S)
    int cur_idx = time.value() / (dt.value());
    // std::cout << time.value() << "\n";
    // std::cout << dt.value() << "\n";
    // std::cout <<  time.value() / (dt.value()) << "\n";
    // std::cout <<  cj[cur_idx] << "\n";
    // std::cout << cur_idx << "\n";


    int h = 1;
    if ( time == tf ) {
        h = -1;
    }

    adouble u0 = controls[0];

    adouble u1 = (-h*cj[cur_idx] + h*cj[cur_idx+h]) / (dt);
    // std::cout << u1 << "\n";

    // https://en.wikipedia.org/wiki/Finite_difference # Higher-order differences
    // second-order forward/backward
    adouble u2 = (cj[cur_idx+2*h] - 2*cj[cur_idx+h] + cj[cur_idx]) / (dt*dt);

    // states[u0_state] = u1; //offset by 1 since
    // states[u1_state] = u2;

    // derivatives[ u0_state ] = u1;
    // derivatives[ u1_state ] = u2;

    //from the PSOPT implementation of get_state_derivative,


    adouble x0_v = states[ x0_v_state ];
    adouble x1_v = states[ x1_v_state ];

    derivatives[ x0_v_state ] = x1_v; // m.Equation(x1_v_v==x0_v_v.dt())
    derivatives[ x1_v_state ] = ((U0 / (T0*T0))*virus->alpha*u2 + (U0 / T0)*virus->beta*u1 + virus->gamma*U0*u0 - virus->phi*(X0 / T0)*x1_v - virus->xi*X0*x0_v)/(X0 / (T0*T0));

    //
    adouble x0_h    = states[ x0_h_state ];
    adouble x1_h    = states[ x1_h_state ];

    derivatives[ x0_h_state ] = x1_h; // m.Equation(x1_v==x0_v.dt())
    derivatives[ x1_h_state ] = ((U0 / (T0*T0))*host->alpha*u2 + (U0 / T0)*host->beta*u1 + host->gamma*U0*u0 - host->phi*(X0 / T0)*x1_h - host->xi*X0*x0_h)/(X0 / (T0*T0));

    if(!hit_store[cur_idx]){
        u0_store[cur_idx] = u0.value();
        u1_store[cur_idx] = u1.value();
        u2_store[cur_idx] = u2.value();
        hit_store[cur_idx] = 1;
    }

    // path[ 0 ] = //path constraint unused here
    // path[ 0 ] = (U0 / (T0*T0))*host->alpha * u2 + (U0 / T0)* host->beta * u1 + host->gamma * u0;
}

static std::string eigentoString(const Eigen::MatrixXd& mat){
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

double normalized_gaussian_pulse(double t){
    t = t - ((1e-8/T0)/2);//-8
    double fwhm = 1e-9/T0; //-9
    double sigma = fwhm/2.355;
    return exp(-((t*t)/(2.0*(sigma*sigma))));
}

double ssin(double t){
    return sin((t/1e-6) * 2*3.14);
}

#define TEST 0

//SPRAL is apparently the best, needs
//https://github.com/jump-dev/Ipopt.jl/issues/223

int main(void)
{

    Alg  algorithm;
    Sol  solution;
    Prob problem;
    problem.name        		          = "Virus_optimization";
    problem.outfilename                 = "virus_optimize.txt";

    problem.nphases   			          = 1;
    problem.nlinkages                   = 0;
    psopt_level1_setup(problem);


    problem.phases(1).nstates   		= 6;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 6;
    problem.phases(1).npath         = 0;
    int nnodes    			             = 300;

    problem.phases(1).nodes         << nnodes;

    psopt_level2_setup(problem, algorithm);



    //in the twoburn.cxx example, rk4_propagate takes a NULL.
    workspace = new Workspace{problem, algorithm,solution};
    hit_store.resize(nnodes);

    u0_store.resize(nnodes);
    u1_store.resize(nnodes);
    u2_store.resize(nnodes);


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////


    virus = new Cell{0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9};
    virus->init();

    host = new Cell{0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9};
    host->init();

    double end_time = 1e-8 / T0;

    double control_bounds = 100;

    double output_bounds = 4;
    double derivative_scaling = 40;
    double second_derivative_scaling = 40;

    //bounds are questionable.
    problem.phases(1).bounds.lower.states << -control_bounds, -derivative_scaling, -output_bounds, -derivative_scaling, -output_bounds, -derivative_scaling;
    problem.phases(1).bounds.upper.states << control_bounds, derivative_scaling, output_bounds, derivative_scaling, output_bounds, derivative_scaling;

    problem.phases(1).bounds.lower.controls << -0.5;
    problem.phases(1).bounds.upper.controls << 0.5;

    // double x0_initial_value = 0.0;
    // double u0_initial_value = 0.0;
    // double u0_integral_constraint = 0;
    // double u0_integral_constraint = 0;

    problem.phases(1).bounds.lower.events << 0,0,0,0,0,0; //2
    problem.phases(1).bounds.upper.events << 0,0,0,0,0,0;

    // problem.phases(1).bounds.lower.events << 0,0,0; //2
    // problem.phases(1).bounds.upper.events << 0,0,0;


    // problem.phases(1).bounds.lower.path << 0.0;
    // problem.phases(1).bounds.upper.path << 0.0;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = end_time;
    problem.phases(1).bounds.upper.EndTime      = end_time;



    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem functions  ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	   = &endpoint_cost;
    problem.dae             	= &dae;
    problem.events 		      = &events;
    problem.linkages		      = &linkages;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////


    //just define control bounds not an energy integral!



    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd u_guess    =  zeros(ncontrols,nnodes);

    // u_guess(0,nnodes/2) = 1;
    MatrixXd x_guess    =  zeros(nstates,nnodes);
    // MatrixXd u_guess =  MatrixXd::Random(ncontrols,nnodes);
    // MatrixXd x_guess = MatrixXd::Random(ncontrols,nnodes);

    // u_guess.
    // MatrixXf::
    //
    MatrixXd time_guess =  linspace(0.0,end_time,nnodes);
    // MatrixXd u_guess    = RandomGaussian(ncontrols, nnodes);


    MatrixXd guess_controls = time_guess.unaryExpr(&normalized_gaussian_pulse)*0.5;
    u_guess = guess_controls;
    // MatrixXd guess_controls_d1 = zeros(problem.phases(1).ncontrols,nnodes);
    // MatrixXd guess_controls_d2 = zeros(problem.phases(1).ncontrols,nnodes);
    //
    // for (int i=1;i<nnodes-1;i++){ //take first gaussian derivative, central differences
    //     guess_controls_d1(i)=(guess_controls(i+1)-guess_controls(i-1))/(2.0*(end_time/nnodes));
    // }
    // for (int i=1;i<nnodes-1;i++){ //take second gaussian derivative (since u2 is our control, not u0!)
    //     guess_controls_d2(i)=(guess_controls_d1(i+1)-guess_controls_d1(i-1))/(2.0*(end_time/nnodes));
    // }

    problem.phases(1).guess.controls       = u_guess; //guess_controls_d2
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = time_guess;


    //try just _h and _v?

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter algorithm options  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    algorithm.nlp_iter_max                = 100;
    // algorithm.nlp_tolerance               = 1e-15;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    // algorithm.collocation_method          = "Legendre";
    algorithm.collocation_method          = "trapezoidal";
    // algorithm.collocation_method          ="Hermite-Simpson";
    algorithm.mesh_refinement             = "automatic";

    // RowVectorXi my_vector(1);
    // my_vector  << nnodes;
    // problem.phases(1).nodes = my_vector;

    algorithm.mr_max_iterations = 1;
    // algorithm.mr_M1 = 30;

    // algorithm.ode_tolerance             = 1.e-13;//increases mesh refinement depth - relative

    // algorithm.nsteps_error_integration  = 20;
    // // algorithm.mr_kappa = 0.4;
    // algorithm.mr_max_increment_factor = 0.05;
    // // algorithm.mr_M1 = 40;
    // algorithm.mr_initial_increment = 50;

    //spral yields an "Emeergency mode" error if
    // export OMP_CANCELLATION=TRUE
    // export OMP_NESTED=TRUE
    // export OMP_PROC_BIND=TRUE
    // is not set

    //see devel/doc/options.dox
    algorithm.ipopt_linear_solver = "ma57";

    // algorithm.ipopt_solver_GPU = 0;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////       Do a test run       //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    if(TEST){
        int test_nnodes = 200;
        double test_end_time = 1e-8;
        MatrixXd initial_test_state    =  zeros(problem.phases(1).nstates,1);
        MatrixXd test_parameters    =  ones(0,1);
        MatrixXd test_state_trajectory    =  zeros(problem.phases(1).nstates,test_nnodes);
        MatrixXd test_time_vector =  linspace(0.0,test_end_time,test_nnodes);
        MatrixXd test_controls = test_time_vector.unaryExpr(&normalized_gaussian_pulse);

        MatrixXd test_controls_derivative_1 = zeros(problem.phases(1).ncontrols,test_nnodes);
        MatrixXd test_controls_derivative_2 = zeros(problem.phases(1).ncontrols,test_nnodes);

        for (int i=1;i<test_nnodes-1;i++){ //take first gaussian derivative, central differences
            test_controls_derivative_1(i)=(test_controls(i+1)-test_controls(i-1))/(2.0*(test_end_time/test_nnodes));
        }
        for (int i=1;i<test_nnodes-1;i++){ //take second gaussian derivative (since u2 is our control, not u0!)
            test_controls_derivative_2(i)=(test_controls_derivative_1(i+1)-test_controls_derivative_1(i-1))/(2.0*(test_end_time/test_nnodes));
        }


        rk4_propagate( dae,
            test_controls_derivative_2,
            test_time_vector,
            initial_test_state,
            test_parameters,
            problem,
            1,
            test_state_trajectory,
            workspace);

        plot(test_time_vector,test_state_trajectory.row(x0_v_state),problem.name, "time (s)", "states", "x0_v");
        plot(test_time_vector,test_state_trajectory.row(x0_h_state),problem.name, "time (s)", "states", "x0_h");
        plot(test_time_vector,test_controls_derivative_2,problem.name, "time (s)", "states", "x");
        return 0;
    }

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Now call PSOPT to solve the problem   /////////////////
    ////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////


    MatrixXd states    = solution.get_states_in_phase(1);
    MatrixXd controls  = solution.get_controls_in_phase(1);
    MatrixXd t         = solution.get_time_in_phase(1);


    MatrixXd x0_v = states.row(x0_v_state);
    MatrixXd x0_h = states.row(x0_h_state);
    // MatrixXd u0 = states.row(u0_state);
    // MatrixXd u1 = states.row(u1_state);
    // MatrixXd u2 = controls.row(0);

    MatrixXd u0 = controls.row(0);
    MatrixXd u1 = states.row(u0_state);
    MatrixXd u2 = states.row(u1_state);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Save solution data to files if desired ////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Save(x0_v,"x0_v.dat");
    Save(x0_h,"x0_h.dat");
    Save(u0, "u0.dat");
    Save(u1, "u1.dat");
    Save(u2, "u2.dat");
    Save(t, "t.dat");



    ////////////////////////////////////////////////////////////////////////////
    ///////////  Plot some results if desired (requires gnuplot) ///////////////
    ////////////////////////////////////////////////////////////////////////////

    MatrixXd sum = (U0 / (T0*T0))*host->alpha * u2 + (U0 / T0)* host->beta * u1 + host->gamma * u0;

    plot(t * T0,x0_v*X0,problem.name, "time (s)", "states", "x0_v");
    plot(t * T0,x0_h*X0,problem.name, "time (s)", "states", "x0_h");
    plot(t * T0,u0,problem.name, "time (s)", "states", "u0");
    Eigen::MatrixXd u0_ = Eigen::Map<Eigen::MatrixXd>(u0_store.data(), 1, nnodes);
    Eigen::MatrixXd u1_ = Eigen::Map<Eigen::MatrixXd>(u1_store.data(), 1, nnodes);
    Eigen::MatrixXd u2_ = Eigen::Map<Eigen::MatrixXd>(u2_store.data(), 1, nnodes);
    plot(t * T0,u0_,problem.name, "time (s)", "states", "u0_2");

    plot(t * T0,u1_,problem.name, "time (s)", "states", "u1");
    plot(t * T0,u2_,problem.name, "time (s)", "states", "u2");
    plot(t * T0,sum,problem.name, "time (s)", "states", "sum");

    Eigen::MatrixXd hit_store_ = Eigen::Map<Eigen::MatrixXd>(hit_store.data(), 1, nnodes);
    plot(t * T0,hit_store_,problem.name, "time (s)", "states", "sum");

}
