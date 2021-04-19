#include "psopt.h"

const double epsilon_0 = 8.854e-12;


#define u0_state 0
#define u1_state 1

#define x0_v_state 2
#define x1_v_state 3

#define x0_h_state 4
#define x1_h_state 5


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


    // double tau_1;
    // double tau_2;


    void init();

};

void Cell::init(){
    double e_o = extracellular_permittivity * epsilon_0; // S/m
    double e_i = intracellular_permittivity * epsilon_0; //S/m
    double e_m = membrane_permittivity * epsilon_0; //S/m
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

    b1 = 2.0 * (R*R*R) * (l_m +     2.0*l_o) * (l_m + 0.5 * l_i) + 2.0 * ((R-d)*(R-d)*(R-d)) * (l_m - l_o) * (l_i - l_m);

    b2 = 2.0 * (R*R*R) * (l_i * (0.5 * e_m + e_o) + l_m * (0.5*e_i + 2.0*e_m + 2*e_o) + l_o * (e_i + 2.0 * e_m)) + (2.0 * ((R-d)*(R-d)*(R-d))
    * (l_i * (e_m - e_o) + l_m * (e_i - 2.0*e_m + e_o) - l_o * (e_i - e_m))); // is this truly a multiply, or a cross?

    b3 = 2.0 * (R*R*R) * (e_m + 2.0*e_o) * (e_m + 0.5 * e_i) + 2.0 * ((R-d)*(R-d)*(R-d)) * (e_m - e_o) * (e_i - e_m);

    // tau_1 = tau_1_f(b1, b2, b3); //yes, this is correct; only b's involved.
    // tau_2 = tau_2_f(b1, b2, b3);
}



Cell* virus;
Cell* host;





///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase,Workspace* workspace)
{
   // return -(final_states[ x0_v_state ]);
   return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{

    // return states[ x0_h_state ]/states[ x0_v_state ]; // does not converge
    // return -smooth_fabs(states[ x0_v_state ], 1e-7) + smooth_fabs(states[ x0_h_state ], 1e-7); //also seems to work okay
    return -(states[ x0_v_state ]*states[ x0_v_state ]) + (states[ x0_h_state ]*states[ x0_h_state ]);
    return 0;
}



////////////////////////////////////////////////////////////////////////////
/////////////////// Define the integrand of the integral constraint ///////
////////////////////////////////////////////////////////////////////////////
adouble integrand( adouble* states, adouble* controls, adouble* parameters,
                    adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    //integrand from isoperimetric.cxx
    adouble g;
    adouble u0 = states[ u0_state ];
    g = (u0*u0);
    return g;
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace){

    adouble u0 = initial_states[ u0_state ];
    e[ 0 ] = u0;
    adouble u1 = initial_states[ u1_state ];
    e[ 1 ] = u1;

    adouble x0 = initial_states[ x0_v_state ];
    e[ 2 ] = x0;

    adouble x0_h = initial_states[ x0_h_state ];
    e[ 3 ] = x0_h;

    adouble x1 = initial_states[ x1_v_state ];
    e[ 4 ] = x1;

    adouble x1_h = initial_states[ x1_h_state ];
    e[ 5 ] = x1_h;

    // adouble u0_f = final_states[ 3 ];
    // e[ 3 ] = u0_f;
    // Compute the integral to be constrained
    adouble Q;
    Q = integrate(integrand, xad, iphase, workspace);
    e[ 6 ] = Q;
}



void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace){

    // vector<Cell * > cells = (vector<Cell * >) workspace->user_data;
    // Cell virus = (Cell) cells[0];

    //https://mathoverflow.net/a/87902/176668
    adouble u0   = states[ u0_state ];
    adouble u1   = states[ u1_state ];
    adouble u2 = controls[ 0 ];
    derivatives[ u0_state ] = u1;
    derivatives[ u1_state ] = u2;

    adouble x0_v    = states[ x0_v_state ];
    adouble x1_v    = states[ x1_v_state ];

    derivatives[ x0_v_state ] = x1_v; // m.Equation(x1_v_v==x0_v_v.dt())
    derivatives[ x1_v_state ] = ((virus->R*virus->a1*u2 + virus->R*virus->a2*u1 + virus->R*virus->a3*u0 - virus->b2*x1_v - virus->b3*x0_v)/virus->b1);


    adouble x0_h    = states[ x0_h_state ];
    adouble x1_h    = states[ x1_h_state ];

    derivatives[ x0_h_state ] = x1_h; // m.Equation(x1_v==x0_v.dt())
    derivatives[ x1_h_state ] = ((host->R*host->a1*u2 + host->R*host->a2*u1 + host->R*host->a3*u0 - host->b2*x1_h - host->b3*x0_h)/host->b1);

    // path[ 0 ] = //path constraint unused here
}

static std::string eigentoString(const Eigen::MatrixXd& mat){
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

double normalized_gaussian_pulse(double t){
    t = t - ((1e-8)/2);//-8
    double fwhm = 1e-9; //-9
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
    problem.phases(1).nevents   		= 7;
    problem.phases(1).npath         = 0;
    int nnodes    			             = 150;

    problem.phases(1).nodes         << nnodes;

    psopt_level2_setup(problem, algorithm);


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////


    virus = new Cell{0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9};
    virus->init();

    host = new Cell{0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9};
    host->init();

    // cells.push_back(virus);
    // cells.push_back(host);

    // problem.user_data = (void *) cells;

    double end_time = 1e-8;

    double control_bounds = 2;

    double output_bounds = 10;
    double derivative_scaling = 1.0/(1e-12); //highest permissible derivative value - gets very high!
    double second_derivative_scaling = derivative_scaling*derivative_scaling;

    //bounds are questionable.
    problem.phases(1).bounds.lower.states << -control_bounds, -derivative_scaling, -output_bounds, -derivative_scaling,  -output_bounds, -derivative_scaling;
    problem.phases(1).bounds.upper.states << control_bounds, derivative_scaling, output_bounds, derivative_scaling, output_bounds, derivative_scaling;


    problem.phases(1).bounds.lower.controls << -second_derivative_scaling;
    problem.phases(1).bounds.upper.controls << second_derivative_scaling;

    double x0_initial_value = 0.0;
    double u0_initial_value = 0.0;
    // double u0_integral_constraint = end_time/2.0;
    double u0_integral_constraint = end_time/2.0;

    problem.phases(1).bounds.lower.events << 0,0,0, u0_initial_value, x0_initial_value, x0_initial_value , u0_integral_constraint; //2
    problem.phases(1).bounds.upper.events << 0,0,0, u0_initial_value, x0_initial_value, x0_initial_value, u0_integral_constraint;

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


    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd u_guess    =  ones(ncontrols,nnodes) * 0.001;
    MatrixXd x_guess    =  zeros(nstates,nnodes);
    MatrixXd time_guess =  linspace(0.0,end_time,nnodes);


    // MatrixXd guess_controls = time_guess.unaryExpr(&ssin);
    //
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
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6; //is this relative? I don't think so
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    // algorithm.collocation_method          = "Legendre";
    // algorithm.collocation_method          = "trapezoidal";
    // algorithm.collocation_method          ="Hermite-Simpson";
    algorithm.mesh_refinement             = "automatic";

    // RowVectorXi my_vector(5);
    // my_vector  << 20, 40, 100, 500, 1500;
    // problem.phases(1).nodes = my_vector;

    algorithm.mr_max_iterations = 10;
    // algorithm.mr_M1 = 30;
    algorithm.ode_tolerance             = 1.e-8;//increases mesh refinement depth - relative
    // algorithm.nsteps_error_integration  = 20;
    // // algorithm.mr_kappa = 0.4;
    // algorithm.mr_max_increment_factor = 0.05;
    // // algorithm.mr_M1 = 40;
    // algorithm.mr_initial_increment = 50;

    // algorithm.ipopt_linear_solver = "ma";
    // algorithm.ipopt_linear_solver = "paradiso";
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


        //in the twoburn.cxx example, rk4_propagate takes a NULL.
        unique_ptr<Workspace> workspace_up{ new Workspace{problem, algorithm,solution} };

        rk4_propagate( dae,
            test_controls_derivative_2,
            test_time_vector,
            initial_test_state,
            test_parameters,
            problem,
            1,
            test_state_trajectory,
            workspace_up.get());

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
    MatrixXd u0 = states.row(u0_state);
    MatrixXd u1 = states.row(u1_state);
    MatrixXd u2 = controls.row(0);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Save solution data to files if desired ////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Save(x0_v,"x0_v.dat");
    Save(x0_h,"x0_h.dat");
    Save(u0, "u0.dat");
    Save(u2, "u2.dat");
    Save(t, "t.dat");



    ////////////////////////////////////////////////////////////////////////////
    ///////////  Plot some results if desired (requires gnuplot) ///////////////
    ////////////////////////////////////////////////////////////////////////////

    plot(t,x0_v,problem.name, "time (s)", "states", "x0_v");
    plot(t,x0_h,problem.name, "time (s)", "states", "x0_h");
    plot(t,u0,problem.name, "time (s)", "states", "u0");
    plot(t,u1,problem.name, "time (s)", "states", "u1");
    plot(t,u2,problem.name, "time (s)", "states", "u2");

}
