#include "psopt.h"

const double epsilon_0 = 8.854e-12;



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


    double tau_1;
    double tau_2;


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


// typedef struct {
//     double a1_v;
//     double a2_v;
//     double a3_v;
//     double b1_v;
//     double b2_v;
//     double b3_v;
//     double a1_h;
//     double a2_h;
//     double a3_h;
//     double b1_h;
//     double b2_h;
//     double b3_h;
// } Constants;




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
   return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{

    // Constants* C = (Constants*) workspace->user_data;


    return  0;
}



////////////////////////////////////////////////////////////////////////////
/////////////////// Define the integrand of the integral constraint ///////
////////////////////////////////////////////////////////////////////////////
adouble integrand( adouble* states, adouble* controls, adouble* parameters,
                    adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    //integrand from isoperimetric.cxx
    adouble g;
    adouble u0 = states[ 3 ];
    g = u0*u0;
    return g;
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace){

    adouble x0 = initial_states[ 0 ];
    e[ 0 ] = x0;
    adouble u0 = initial_states[ 3 ];
    e[ 1 ] = u0;

    // adouble xf = final_states[ 0 ];

    // Compute the integral to be constrained
    adouble Q;
    Q = integrate(integrand, xad, iphase, workspace);
    e[ 2 ] = Q;
}



void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace){

   // vector<Cell * > cells = (vector<Cell * >) workspace->user_data;
   // Cell virus = (Cell) cells[0];


  //https://mathoverflow.net/a/87902/176668
  adouble u0   = states[ 3 ];
  adouble u1   = states[ 4 ];
  adouble u2 = controls[ 0 ];
  derivatives[ 3 ] = u1;
  derivatives[ 4 ] = u2;

   adouble x0    = states[ 0 ];
   adouble x1    = states[ 1 ];
   //
   // // adouble x2_v    = states[ 2 ];
   derivatives[ 0 ] = x1; // m.Equation(x1_v==x0_v.dt())
   derivatives[ 1 ] = ((virus->R*virus->a1*u2 + virus->R*virus->a2*u1 + virus->R*virus->a3*u0 - virus->b2*x1 - virus->b3*x0)/virus->b1);
   // m.Equation(x2_v==x1_v.dt()) // x2_v = (R_v*a1_v*u2 + R_v*a2_v*u1 + R_v*a3_v*u0 - b2_v*x1_v - b3_v*x0_v)/b1_v);

   // path[ 0 ] = //path constraint unused here
}

static std::string eigentoString(const Eigen::MatrixXd& mat){
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

double normalized_gaussian_pulse(double t){
    t = t - ((1e-8)/2);
    double fwhm = 1e-9;
    double sigma = fwhm/2.355;
    return exp(-((t*t)/(2.0*(sigma*sigma))));
}

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


    problem.phases(1).nstates   		= 5;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 3;
    problem.phases(1).npath         = 0;
    int nnodes    			             = 50;
    problem.phases(1).nodes         << nnodes;

    psopt_level2_setup(problem, algorithm);


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////


    virus = new Cell{0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9};
    virus->init();

    host = new Cell{0.3, 80, 0.3, 80, 1e-7, 5, 50e-9, 5e-9};
    host->init();

    // cells.push_back(virus);
    // cells.push_back(host);

    // problem.user_data = (void *) cells;

    double control_bounds = 1e9;

    double output_bounds = 1e9;

    problem.phases(1).bounds.lower.states << -output_bounds, -output_bounds*100, -output_bounds*100, -control_bounds, -control_bounds; //fix bounds!
    problem.phases(1).bounds.upper.states << output_bounds, output_bounds*100, output_bounds*100, control_bounds, control_bounds; //fix bounds!


    problem.phases(1).bounds.lower.controls << -control_bounds;
    problem.phases(1).bounds.upper.controls << control_bounds;

    double x0_initial_value = 0.0;
    double u0_initial_value = 0.0;
    double u0_integral_constraint = 1.0;

    problem.phases(1).bounds.lower.events << x0_initial_value, u0_initial_value, u0_integral_constraint; //2
    problem.phases(1).bounds.upper.events << x0_initial_value, u0_initial_value, u0_integral_constraint;

    // problem.phases(1).bounds.lower.path << 0.0;
    // problem.phases(1).bounds.upper.path << 0.0;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1e-10;
    problem.phases(1).bounds.upper.EndTime      = 3.0;



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

    MatrixXd u_guess    =  ones(ncontrols,nnodes);
    MatrixXd x_guess    =  ones(nstates,nnodes);
    MatrixXd time_guess =  linspace(0.0,3.0,nnodes);

    //
    // u_guess << linspace(theta0,thetaf,nnodes),
    //            linspace(phi0,phif,nnodes);
    //
    // for (int i = 0;i< nnodes;i++) {
    //
    //   x_guess(0,i) = a*sin(u_guess(0,i))*cos(u_guess(1,i));
    //   x_guess(1,i) = a*sin(u_guess(0,i))*sin(u_guess(1,i));
    //   x_guess(2,i) = b*cos(u_guess(0,i));
    //
    // }


    problem.phases(1).guess.controls       = u_guess;
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = time_guess;


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter algorithm options  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-9;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////       Do a test run       //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    int test_nnodes = 200;
    MatrixXd initial_test_state    =  zeros(problem.phases(1).nstates,1);
    MatrixXd test_parameters    =  ones(0,1);
    MatrixXd test_state_trajectory    =  zeros(problem.phases(1).nstates,test_nnodes);
    MatrixXd test_time_vector =  linspace(0.0,1e-8,test_nnodes);
    MatrixXd test_controls = test_time_vector.unaryExpr(&normalized_gaussian_pulse);
    MatrixXd test_controls_derivative_1 = zeros(problem.phases(1).ncontrols,test_nnodes);
    MatrixXd test_controls_derivative_2 = zeros(problem.phases(1).ncontrols,test_nnodes);

    for (int i=1;i<test_nnodes-1;i++){ //take first gaussian derivative, central differences
        test_controls_derivative_1(i)=(test_controls(i+1)-test_controls(i-1))/2;
    }
    for (int i=1;i<test_nnodes-1;i++){ //take second gaussian derivative (since u2 is our control, not u0!)
        test_controls_derivative_2(i)=(test_controls_derivative_1(i+1)-test_controls_derivative_1(i-1))/2;
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

    // std::cout << eigentoString(test_state_trajectory) + "\n";
    plot(test_time_vector,test_state_trajectory.row(0).normalized(),problem.name, "time (s)", "states", "x y z s b");
    plot(test_time_vector,test_controls_derivative_2,problem.name, "time (s)", "states", "x");

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


    MatrixXd x0 = states.row(0);
    MatrixXd u0 = states.row(3);

    MatrixXd u2 = controls.row(0);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Save solution data to files if desired ////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    // Save(x,"x.dat");
    // Save(y,"y.dat");
    // Save(z, "z.dat");
    // Save(theta,"theta.dat");
    // Save(phi, "phi.dat");
    // Save(t,"t.dat");


    ////////////////////////////////////////////////////////////////////////////
    ///////////  Plot some results if desired (requires gnuplot) ///////////////
    ////////////////////////////////////////////////////////////////////////////

    plot(t,x0.normalized(),problem.name, "time (s)", "states", "x0");
    plot(t,u0.normalized(),problem.name, "time (s)", "states", "u0");

}
