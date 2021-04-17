#include "psopt.h"



typedef struct {
    double a1_v;
    double a2_v;
    double a3_v;
    double b1_v;
    double b2_v;
    double b3_v;
    double a1_h;
    double a2_h;
    double a3_h;
    double b1_h;
    double b2_h;
    double b3_h;
} Constants;





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

    Constants* C = (Constants*) workspace->user_data;

    double V = C->V;

    adouble theta = controls[ 0 ];
    adouble phi   = controls[ 1 ];



    return  L;
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
    // adouble xf = final_states[ 0 ];

    adouble Q;

    // Compute the integral to be constrained
    Q = integrate(integrand, xad, iphase, workspace);

    e[ 0 ] = x0;
    e[ 1 ] = Q;
}



void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

   Constants* C = (Constants*) workspace->user_data;


   adouble x0_v    = states[ 0 ];
   adouble x1_v    = states[ 1 ];
   adouble x2_v    = states[ 2 ];
   derivatives[ 0 ] = x1_v; // m.Equation(x1_v==x0_v.dt())
   derivatives[ 1 ] = x2_v; // m.Equation(x2_v==x1_v.dt())


   //https://mathoverflow.net/a/87902/176668
   adouble u0    = states[ 3 ];
   adouble u1    = states[ 4 ];
   adouble u2 = controls[ 0 ];
   derivatives[ 3 ] = u1;
   derivatives[ 4 ] = u2;


   double a1_v = C->a1_v;
   double a2_v = C->a2_v;
   double a3_v = C->a3_v;
   double b1_v = C->b1_v;
   double b2_v = C->b2_v;
   double b3_v = C->b3_v;

   // double a1_h = C->a1_h;
   // double a2_h = C->a2_h;
   // double a3_h = C->a3_h;
   // double b1_h = C->b1_h;
   // double b2_h = C->b2_h;
   // double b3_h = C->b3_h;


   x2_v == (R_v*a1_v*u2 + R_v*a2_v*u1 + R_v*a3_v*u0 - b2_v*x1_v - b3_v*x0_v)/b1_v)



   // path[ 0 ] = //path constraint unused here
}

double normalized_gaussian_pulse(double t,double fwhm):
    double sigma = fwhm/2.355;
    return exp(-((t**2.0)/(2.0*(sigma**2.0))));


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

    problem.phases(1).nstates   		= 3;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 2;
    problem.phases(1).npath         = 1;
    problem.phases(1).nodes         << 20;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

   Constants* C = new Constants;

   problem.user_data = (void*) C;


   C->V = 900.00; // Speed in km/h
   C->a = 6384.0; // Earth's semi-major axis in km
   C->b = 6353.0; // Earth's semi-minor axis in km


    problem.phases(1).bounds.lower.states << xL, yL, zL;
    problem.phases(1).bounds.upper.states << xU, yU, zU;


    problem.phases(1).bounds.lower.controls << -1;
    problem.phases(1).bounds.upper.controls << 1;

    double x0_initial_value = 0.0;
    double u0_integral_constraint = 1.0;

    problem.phases(1).bounds.lower.events << x0_initial_value, u0_integral_constraint;
    problem.phases(1).bounds.upper.events << x0_initial_value, u0_integral_constraint;

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

    int nnodes    			             = 30;
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd u_guess    =  zeros(ncontrols,nnodes);
    MatrixXd x_guess    =  zeros(nstates,nnodes);
    MatrixXd time_guess =  linspace(0.0,3.0,nnodes);


    u_guess << linspace(theta0,thetaf,nnodes),
               linspace(phi0,phif,nnodes);

    for (int i = 0;i< nnodes;i++) {

      x_guess(0,i) = a*sin(u_guess(0,i))*cos(u_guess(1,i));
      x_guess(1,i) = a*sin(u_guess(0,i))*sin(u_guess(1,i));
      x_guess(2,i) = b*cos(u_guess(0,i));

    }


    problem.phases(1).guess.controls       = u_guess;
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = time_guess;


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-4;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";


    MatrixXd test_time_vector =  linspace(0.0,1e-8,nnodes);
    MatrixXd test_controls = time.binaryExpr(1e-9,&normalized_gaussian_pulse);

    rk4_propagate( &dae,
        MatrixXd& test_controls,
        MatrixXd& test_time_vector,
        MatrixXd& initial_state,
        MatrixXd& parameters,
        Prob & problem,
        int iphase,
        MatrixXd& state_trajectory);





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


    MatrixXd x = states.row(0);
    MatrixXd y = states.row(1);
    MatrixXd z = states.row(2);

    MatrixXd theta = controls.row(0);
    MatrixXd phi   = controls.row(1);

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

	 plot(t,states,problem.name, "time (s)", "states", "x y z");

    plot(t,controls,problem.name, "time (s)", "controls", "theta phi");


	 plot(t,states,problem.name, "time (s)", "states", "x y z",
                           "pdf", "geodesic_states.pdf");

    plot(t,controls,problem.name, "time (s)", "controls", "theta phi",
                           "pdf", "geodesic_controls.pdf");

}








///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}
