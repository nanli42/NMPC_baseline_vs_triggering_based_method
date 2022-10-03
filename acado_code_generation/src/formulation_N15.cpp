#include <acado_code_generation.hpp>

typedef struct CarParams {
    double wheelbase;
    double friction_coeff;
    double h_cg; // height of car's CG
    double l_f; // length from CG to front axle
    double l_r; // length from CG to rear axle
    double cs_f; // cornering stiffness coeff for front wheels
    double cs_r; // cornering stiffness coeff for rear wheels
    double mass;
    double I_z; // moment of inertia about z axis from CG
} Param;

int main( )
{
  USING_NAMESPACE_ACADO

  DifferentialState   ey, ephi, v, beta, omega, delta, t, s;
  IntermediateState   dsdt;
  Control             a, vd;

  OnlineData          ks;

  DifferentialEquation  f;

  double g = 9.81; // m/s^2

  Param p = {
        .wheelbase = 0.3302,
        .friction_coeff = 0.523,
        .h_cg = 0.074,
        .l_f = 0.3302 - 0.17145,
        .l_r = 0.17145,
        .cs_f = 4.718,
        .cs_r = 5.4562,
        .mass = 3.47,
        .I_z = 0.04712
  };

  IntermediateState rear_val = g * p.l_r - a * p.h_cg;
  IntermediateState front_val = g * p.l_f + a * p.h_cg;

  IntermediateState beta_dot = (p.friction_coeff / (v * (p.l_r + p.l_f))) *
          (p.cs_f * delta * (rear_val) -
           beta * (p.cs_r * (front_val) + p.cs_f * (rear_val)) +
           (omega/v) * (p.cs_r * p.l_r * (front_val) - p.cs_f * p.l_f * (rear_val))) -
          omega;

  IntermediateState omega_dot = (p.friction_coeff * p.mass / (p.I_z * p.wheelbase)) *
          (p.l_f * p.cs_f * delta * (rear_val) +
           beta * (p.l_r * p.cs_r * (front_val) - p.l_f * p.cs_f * (rear_val)) -
           (omega/v) * (std::pow(p.l_f, 2) * p.cs_f * (rear_val) + std::pow(p.l_r, 2) * p.cs_r * (front_val)));

  dsdt = v*cos(beta + ephi) / ( 1 - ey*ks );

  f << dot(ey)      == 1/dsdt * ( v*sin(beta + ephi) );
  f << dot(ephi)    == 1/dsdt * omega - ks;
  f << dot(v)       == 1/dsdt * a;
  f << dot(beta)    == 1/dsdt * beta_dot;
  f << dot(omega)   == 1/dsdt * omega_dot;
  f << dot(delta)   == 1/dsdt * vd;
  f << dot(t)       == 1/dsdt;
  f << dot(s)       == 1;

  const double s_start =  0.0;
  const double s_end   =  6.0;
  OCP ocp( s_start, s_end, 20 );
  ocp.minimizeMayerTerm( t );
  ocp.setModel( f );

  double max_v = 5.0;
  ocp.subjectTo( -1.3 <= ey <= +1.3 );
  ocp.subjectTo( -0.75 <= ephi <= +0.75 );
  ocp.subjectTo( +0.05 <= v <= +max_v );
  ocp.subjectTo( -0.4 <= beta <= +0.4 ); // inferred value
  ocp.subjectTo( -4.0 <= omega <= +4.0 ); // inferred value
  ocp.subjectTo( -0.41 <= delta <= +0.41 );
  ocp.subjectTo( +0.0 <= t <= 1000.0 );
  ocp.subjectTo( +0.0 <= s <= 1000.0 );

  ocp.subjectTo( -13.26 <= a <= +9.51 );
  ocp.subjectTo( -3.2 <= vd <= +3.2 );

  IntermediateState v_para = v*cos(beta + ephi);
  IntermediateState a_vert = v_para*v_para/(1/ks-ey);
  IntermediateState a_para = a*cos(beta + ephi);
  IntermediateState a_total_2 = a_para*a_para + a_vert*a_vert;

  ocp.subjectTo( a_total_2 <= p.friction_coeff*g * p.friction_coeff*g );

  ocp.setNOD(1);

  OCPexport mpc( ocp );

  mpc.set( HESSIAN_APPROXIMATION,       EXACT_HESSIAN     );
  mpc.set( DISCRETIZATION_TYPE,         MULTIPLE_SHOOTING );
  mpc.set( INTEGRATOR_TYPE,             INT_RK4           );
  mpc.set( NUM_INTEGRATOR_STEPS,        500               );
  mpc.set( QP_SOLVER,                   QP_QPOASES        );
  mpc.set( HOTSTART_QP,                 NO               	);
  mpc.set( GENERATE_TEST_FILE,          NO                );
  mpc.set( GENERATE_MAKE_FILE,          NO                );
  mpc.set( GENERATE_MATLAB_INTERFACE,   NO                );
  mpc.set( SPARSE_QP_SOLUTION, 		      FULL_CONDENSING_N2);
  mpc.set( DYNAMIC_SENSITIVITY, 		  	SYMMETRIC				  );

  mpc.set( GENERATE_SIMULINK_INTERFACE, NO                );
  mpc.set( CG_HARDCODE_CONSTRAINT_VALUES, NO 					    );
  mpc.set( CG_USE_VARIABLE_WEIGHTING_MATRIX, YES 				  );

  //mpc.set( LEVENBERG_MARQUARDT, 1e-6 );

  if (mpc.exportCode( "../../src/acado_generated_code_N15" ) != SUCCESSFUL_RETURN)
    exit( EXIT_FAILURE );

  mpc.printDimensionsQP( );

  return EXIT_SUCCESS;
}
