#include <acado_code_generation.hpp>

int main( )
{
  USING_NAMESPACE_ACADO

  DifferentialState   ey, ephi, vx, vy, w, d, delta, t, s;//, dist; //, dummy;
  IntermediateState   dsdt, Frx, Fry, Ffy, alphaF, alphaR;
  Control             dd, ddelta; //, slack;

  OnlineData          ks;//, vx_obs, ey_obs;

  DifferentialEquation  f;

  const double m = 0.041;
  const double Iz = 27.8e-6;
  const double lf = 0.029;
  const double lr = 0.033;

  const double Cm1=0.287;
  const double Cm2=0.0545;
  const double Cr0=0.0518;
  const double Cr2=0.00035;

  const double Br = 3.3852;
  const double Cr = 1.2691;
  const double Dr = 0.1737;

  const double Bf = 2.579;
  const double Cf = 1.2;
  const double Df = 0.192;

  const double e_long = 0.9;
  const double e_eps = 0.95;

  const double min_alphaF = -0.6;
  const double max_alphaF = 0.6;

  dsdt = ( vx*cos(ephi) - vy*sin(ephi) ) / ( 1 - ey*ks );

  alphaF = -atan((w*lf+vy)/vx) + delta;
  alphaR = atan((w*lr-vy)/vx);

  f << dot(ey)      == 1/dsdt * ( vx*sin(ephi) + vy*cos(ephi) );
  f << dot(ephi)    == 1/dsdt * w - ks;

  /*
  IntermediateState a = (Cm1-Cm2*vx)*d/m;
  f << dot(vx)      == 1/dsdt * a;
  f << dot(vy)      == 1/dsdt * ( lr/(lr+lf) * ( ddelta*vx + delta*a) );
  f << dot(w)       == 1/dsdt * (  1/(lr+lf) * ( ddelta*vx + delta*a) );
  */

  Frx = (Cm1-Cm2*vx)*d - Cr0 - Cr2*vx*vx;
  //Frx = (Cm1-Cm2*vx)*d;
  Fry = Dr*sin(Cr*atan(Br*alphaR));
  Ffy = Df*sin(Cf*atan(Bf*alphaF));
  f << dot(vx)      == 1/dsdt * (w*vy + 1/m*(Frx-Ffy*sin(delta)));
  f << dot(vy)      == 1/dsdt * (-w*vx + 1/m*(Fry+Ffy*cos(delta)));
  f << dot(w)       == 1/dsdt * (1/Iz * (lf*Ffy*cos(delta)-lr*Fry));

  f << dot(d)       == 1/dsdt * dd;
  f << dot(delta)   == 1/dsdt * ddelta;
  f << dot(t)       == 1/dsdt;
  f << dot(s)       == 1;
  //f << dot(dist)    == 1/dsdt * (-vx*cos(ephi)+vy*sin(ephi)+vx_obs);

  const double s_start =  0.0;
  const double s_end   =  0.9;
  OCP ocp( s_start, s_end, 15 );
  ocp.minimizeMayerTerm( t );
  //ocp.minimizeLagrangeTerm( exp(1e2(-dist*dist-(ey-ey_obs)*(ey-ey_obs) + 0.06*0.06)) );
  ocp.setModel( f );

  //double dey = 0.045;
  double dey = 0.035;
  double max_v = 1.6;
  ocp.subjectTo( -0.17 <= ey <= +0.17 );
  ocp.subjectTo( -1.5 <= ephi <= +1.5 );
  ocp.subjectTo( +0.05 <= vx <= +max_v );
  ocp.subjectTo( -1.0 <= vy <= +1.0 );
  ocp.subjectTo( -8.0 <= w <= +8.0 );
  ocp.subjectTo( -1.0 <= d <= +1.0 );
  ocp.subjectTo( -0.6 <= delta <= +0.6 );
  ocp.subjectTo( +0.0 <= t <= 100.0 );

  ocp.subjectTo( -10.0 <= dd <= +10.0 );
  ocp.subjectTo( -10.0 <= ddelta <= +10.0 );

  ocp.subjectTo( vx*vx+vy*vy <= max_v*max_v );

  ocp.subjectTo( (e_long*Frx)*(e_long*Frx) + Fry*Fry - (e_eps*Dr)*(e_eps*Dr) <= 0.0 );
  //ocp.subjectTo( (e_long*Frx)*(e_long*Frx) - (e_eps*Dr)*(e_eps*Dr) <= 0.0 );
  ocp.subjectTo( min_alphaF <= alphaF <= max_alphaF );
/*
  ocp.subjectTo( 0.0 <= slack <= 0.055 );
  ocp.subjectTo( ey + dey - slack <= +0.17 );
  ocp.subjectTo( -0.17 <= ey - dey + slack );
*/
  ocp.subjectTo( ey + dey <= +0.17 );
  ocp.subjectTo( -0.17 <= ey - dey );

  //ocp.setNOD(3);
  ocp.setNOD(1);

  OCPexport mpc( ocp );

  mpc.set( HESSIAN_APPROXIMATION,       EXACT_HESSIAN     );
  mpc.set( DISCRETIZATION_TYPE,         MULTIPLE_SHOOTING );
  mpc.set( INTEGRATOR_TYPE,             INT_RK4           );
  mpc.set( NUM_INTEGRATOR_STEPS,        100               );
  mpc.set( QP_SOLVER,                   QP_QPOASES        );
  mpc.set( HOTSTART_QP,                 NO               	);
  mpc.set( GENERATE_TEST_FILE,          NO                );
  mpc.set( GENERATE_MAKE_FILE,          NO                );
  mpc.set( GENERATE_MATLAB_INTERFACE,   NO                );
  mpc.set( SPARSE_QP_SOLUTION, 		      FULL_CONDENSING_N2);
  //mpc.set( DYNAMIC_SENSITIVITY, 		  	FORWARD_OVER_BACKWARD				  );
  mpc.set( DYNAMIC_SENSITIVITY, 		  	SYMMETRIC				  );

  mpc.set( GENERATE_SIMULINK_INTERFACE, NO                );
  mpc.set( CG_HARDCODE_CONSTRAINT_VALUES, NO 					    );
  mpc.set( CG_USE_VARIABLE_WEIGHTING_MATRIX, YES 				  );

  mpc.set( LEVENBERG_MARQUARDT, 1e-6 );

  if (mpc.exportCode( "../../src/acado_generated_code_N15" ) != SUCCESSFUL_RETURN)
    exit( EXIT_FAILURE );

  mpc.printDimensionsQP( );

  return EXIT_SUCCESS;
}
