/*
 *    This file was auto-generated using the ACADO Toolkit.
 *
 *    While ACADO Toolkit is free software released under the terms of
 *    the GNU Lesser General Public License (LGPL), the generated code
 *    as such remains the property of the user who used ACADO Toolkit
 *    to generate this code. In particular, user dependent data of the code
 *    do not inherit the GNU LGPL license. On the other hand, parts of the
 *    generated code that are a direct copy of source code from the
 *    ACADO Toolkit or the software tools it is based on, remain, as derived
 *    work, automatically covered by the LGPL license.
 *
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */



/*

IMPORTANT: This file should serve as a starting point to develop the user
code for the OCP solver. The code below is for illustration purposes. Most
likely you will not get good results if you execute this code without any
modification(s).

Please read the examples in order to understand how to write user code how
to run the OCP solver. You can find more info on the website:
www.acadotoolkit.org

*/

#include "acado_common.h"
#include "acado_auxiliary_functions.h"

#include "system_config.h"
#include "track.h"

#include <stdio.h>
#include <stdlib.h>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

/* Some convenient definitions. */
#define NX          ACADO_NX  /* Number of differential state variables.  */
#define NXA         ACADO_NXA /* Number of algebraic variables. */
#define NU          ACADO_NU  /* Number of control inputs. */
#define NOD         ACADO_NOD  /* Number of online data values. */

#define NY          ACADO_NY  /* Number of measurements/references on nodes 0..N - 1. */
#define NYN         ACADO_NYN /* Number of measurements/references on node N. */

#define N           ACADO_N   /* Number of intervals in the horizon. */

#define VERBOSE     1         /* Show iterations: 1, silent: 0.  */

/* Global variables used by the solver. */
ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

using namespace nmpc4rc;

#define NUM_STEPS 61    /* Number of real-time iterations. */
#define TRIGGER_MODE_1 1 /* 1 - curvature, 2 - uncertainty */
#define TRIGGER_MODE_2 1 /* 1 - curvature, 2 - uncertainty */
#define TRIGGER_CURVATURE 0.4
#define TRIGGER_PROGRESS_TIME 0.1

double get_curature(ArcLengthSpline sp, double s) {
    double dx = sp.getDerivative(s)(0);
    double dy = sp.getDerivative(s)(1);
    double ddx = sp.getSecondDerivative(s)(0);
    double ddy = sp.getSecondDerivative(s)(1);
    double k = (dx*ddy-dy*ddx)/std::pow((dx*dx+dy*dy), 3/2);

    if (abs(k)<1e-5) {
      return 1e-5;
    }

    if (abs(k)>5) {
      if (k < 0) k = -5;
      else k = 5;
    }

    return k;
};

//////////////////////////////////////////////////////////////////////////////
////////                        Boost integration                     ////////
//////////////////////////////////////////////////////////////////////////////
using namespace boost::numeric::odeint;
double m = 0.041;double Iz = 27.8e-6;double lf = 0.029;double lr = 0.033;
double Cm1=0.287;double Cm2=0.0545;double Cr0=0.0518;double Cr2=0.00035;
double Br = 3.3852;double Cr = 1.2691;double Dr = 0.1737;
double Bf = 2.579;double Cf = 1.2;double Df = 0.192;

typedef std::vector<double> state_type;

void obs_odes_in_term_of_s( const state_type &x , state_type &dxds , double s) {
  double dsdt = (x[2]*cos(x[1])-x[3]*sin(x[1])) / (1-x[0]*x[11]); //x[13]);
  double alphaF = -atan((x[4]*lf+x[3])/x[2]) + x[6];
  double alphaR = atan((x[4]*lr-x[3])/x[2]);

  double Frx = (Cm1-Cm2*x[2])*x[5] - Cr0 - Cr2*x[2]*x[2];
  double Fry = Dr*sin(Cr*atan(Br*alphaR));
  double Ffy = Df*sin(Cf*atan(Bf*alphaF));

  dxds[0] = 1/dsdt * (x[2]*sin(x[1]) + x[3]*cos(x[1]));
  dxds[1] = 1/dsdt * x[4] - x[11];
  dxds[2] = 1/dsdt * (x[4]*x[3] + 1/m*(Frx-Ffy*sin(x[6])));
  dxds[3] = 1/dsdt * (-x[4]*x[2] + 1/m*(Fry+Ffy*cos(x[6])));
  dxds[4] = 1/dsdt * (1/Iz * (lf*Ffy*cos(x[6])-lr*Fry));
  dxds[5] = 1/dsdt * x[9]; //x[11];
  dxds[6] = 1/dsdt * x[10]; //x[12];
  dxds[7] = 1/dsdt * 1;

  dxds[8] = 1;

  dxds[9] = 0;
  dxds[10] = 0;

  dxds[11] = 0;
}
State integrate_in_term_of_s(State state, double delta_s) {
  state_type x(12); // = { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, // 9 states
            //   1.0, 0.0, // 2 ctrls
            //   0.0 // 1 curvature
            // }; // initial conditions // 9 states + 2 ctrls + 1 curvature

  x[0] = state.ey;
  x[1] = state.epsi;
  x[2] = state.vx;
  x[3] = state.vy;
  x[4] = state.w;
  x[5] = state.d;
  x[6] = state.delta;
  x[7] = state.t;
  x[8] = state.s;

  x[9] = state.dd;
  x[10] = state.ddelta;

  x[11] = state.curvature;

  typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
  size_t steps =  integrate_adaptive( make_controlled< error_stepper_type >( 1.0e-10 , 1.0e-6 ) , // abs_err = 1.0e-10 , rel_err = 1.0e-6
                    obs_odes_in_term_of_s, x, 0.0, delta_s, 1e-6 );  //...start_time, end_time, default_dt...

  State result_state;
  result_state.ey = x[0];
  result_state.epsi = x[1];
  result_state.vx = x[2];
  result_state.vy = x[3];
  result_state.w = x[4];
  result_state.d = x[5];
  result_state.delta = x[6];
  result_state.t = x[7];
  result_state.s = x[8];

  return result_state;
};

/* A template for testing of the solver. */
int main( )
{
  Track track;
  TrackPos track_xy = track.loadTrack();
  track.fitTrackIntoSpline(track_xy);

	/* Some temporary variables. */
	int    i, j, step;
	acado_timer t;

  FILE * fp;
  FILE * fp_actual;
  fp = fopen ("../data/traj.txt", "w");
  fp_actual = fopen ("../data/traj_actual.txt", "w");

	int next_n = 1;

	/* Initialize the solver. */
	acado_initializeSolver();

  acadoVariables.x0[0] = 0.0; // ey
  acadoVariables.x0[1] = 0.0; // ehi
  acadoVariables.x0[2] = 1.0; // vx
  acadoVariables.x0[3] = 0.0; // vy
  acadoVariables.x0[4] = 0.0; // omega
  acadoVariables.x0[5] = 0.0; // d
  acadoVariables.x0[6] = 0.0; // delta
  acadoVariables.x0[7] = 0.0; // t
  acadoVariables.x0[8] = 0.0; // s

  for (i = 0; i < N; ++i) acadoVariables.od[i] = get_curature(track.sp, acadoVariables.x0[8] + i*0.06);

  for (i = 0; i < N + 1; ++i)  {
    for (j = 0; j < NX; ++j)
      acadoVariables.x[ i*NX + j ] = 0.0;
    acadoVariables.x[ i*NX + 2 ] = 1.0;
    acadoVariables.x[ i*NX + 5 ] = 1.0;
    acadoVariables.x[ i*NX + 6 ] = atan(acadoVariables.od[i]*0.06);
  }

  for (i = 0; i < NU * N; ++i)  acadoVariables.u[ i ] = 0.0;

	/* The "real-time iterations" loop. */
	for(step = 0; step < NUM_STEPS; ++step)
	{
		/* Get the time before start of the loop. */
		acado_tic( &t );

		double prec = 1e+6;
		int iter = 0;
		while(prec > 1e-4 && iter < 20 || prec == 0.0) {
			/* Prepare first step */
			acado_preparationStep();
			/* Perform the feedback step. */
			acado_feedbackStep();
			prec = acado_getKKT();
			iter++;
		}

if (TRIGGER_MODE_1) {
    for (i = 1; i < N; ++i)
      if ((abs(get_curature(track.sp, acadoVariables.x0[8] + (N+i)*0.06)-get_curature(track.sp, acadoVariables.x0[8] + N*0.06))
      > TRIGGER_CURVATURE) && (acadoVariables.x[NX*i+7]-acadoVariables.x0[7]>TRIGGER_PROGRESS_TIME)) {
        double a = get_curature(track.sp, acadoVariables.x0[8] + (N+i)*0.06);
        double b = get_curature(track.sp, acadoVariables.x0[8] + N*0.06);
        //printf("%f, %f - delta %f", a, b, a-b);
        break;
    }
    next_n = i;
}

		/* Read the elapsed time. */
		real_t te = acado_toc( &t );

		/* Apply the new control immediately to the process, first NU components. */
		if( VERBOSE ) {
			printf("\tStep %d: Total iter = %d, KKT Tolerance = %.3e\n", step, iter, acado_getKKT() );
      printf("\tprogress:   %.3f\n", acadoVariables.x0[8]);
      printf("\ttime:   %.3g ms\n\n", 1e3 * te);
		}

		for (i = 0; i < N + 1; ++i)
			fprintf(fp, "%f %f %f %f %f %f %f %f %f %f\n",
				acadoVariables.x[i*NX + 8], // s
				acadoVariables.x[i*NX + 0], // ey
				acadoVariables.x[i*NX + 1], // ephi
				sqrt(acadoVariables.x[i*NX + 2]*acadoVariables.x[i*NX + 2] + acadoVariables.x[i*NX + 3]*acadoVariables.x[i*NX + 3]), // v
        acadoVariables.x[i*NX + 5], // d
        acadoVariables.x[i*NX + 6], // delta
        acadoVariables.x[i*NX + 7], // t
				acadoVariables.od[i], // kappa
        i == N ? 0.0 : acadoVariables.u[i*NU + 0], // dd
        i == N ? 0.0 : acadoVariables.u[i*NU + 1] // ddelta
			);

    printf("next_n: %d\n", next_n);

    State state = {
      acadoVariables.x[0*NX+0],
      acadoVariables.x[0*NX+1],
      acadoVariables.x[0*NX+2],
      acadoVariables.x[0*NX+3],
      acadoVariables.x[0*NX+4],
      acadoVariables.x[0*NX+5],
      acadoVariables.x[0*NX+6],
      acadoVariables.x[0*NX+7],
      acadoVariables.x[0*NX+8],
      0.0,
      0.0,
      acadoVariables.u[0*NU+0],
      acadoVariables.u[0*NU+1],
      acadoVariables.od[0]
    };

if (TRIGGER_MODE_2) {

    for (i = 0; i < next_n; ++i) {
      fprintf(fp_actual, "%f %f %f %f %f %f %f\n",
        state.s,
        state.ey,
        state.epsi,
        sqrt(state.vx*state.vx+state.vy*state.vy),
        state.d,
        state.delta,
        state.t
      );

      /*
      state.ddelta += (rand() / (RAND_MAX + 1.)-0.5)*2 * 1.0;
      if (state.ddelta > 10.0) state.ddelta = 10.0;
      if (state.ddelta < -10.0) state.ddelta = -10.0;
      */
      state = integrate_in_term_of_s(state, 0.06);

      if (i>0) {
        if (abs(state.ey-acadoVariables.x[(i+1)*NX+0])>0.01 ||
            abs(state.vx-acadoVariables.x[(i+1)*NX+2])>0.20
          ) {
          next_n = (i+1);
          state.dd = acadoVariables.u[(i+1)*NU+0];
          state.ddelta = acadoVariables.u[(i+1)*NU+1];
          state.curvature = acadoVariables.od[i+1];
          state = integrate_in_term_of_s(state, 0.06);
          printf("next_n change to: %d\n", next_n);
          break;
        }
      }

      /*
      state.ey += (rand() / (RAND_MAX + 1.)-0.5)*2 * 0.02;
      if (state.ey > 0.135) state.ey = 0.135;
      if (state.ey < -0.135) state.ey = -0.135;
      */
      /*
      printf("%f\n", state.ey);
      printf("%f\n", state.epsi);
      printf("%f\n", state.vx);
      printf("%f\n", state.vy);
      printf("%f\n", state.w);
      printf("%f\n", state.d);
      printf("%f\n", state.delta);
      printf("%f\n", state.t);
      printf("%f\n", state.s);
      printf("---\n");
      */
      state.dd = acadoVariables.u[(i+1)*NU+0];
      state.ddelta = acadoVariables.u[(i+1)*NU+1];
      state.curvature = acadoVariables.od[i+1];
    }

    acadoVariables.x0[0] = state.ey;
    acadoVariables.x0[1] = state.epsi;
    acadoVariables.x0[2] = state.vx;
    acadoVariables.x0[3] = state.vy;
    acadoVariables.x0[4] = state.w;
    acadoVariables.x0[5] = state.d;
    acadoVariables.x0[6] = state.delta;
    acadoVariables.x0[7] = state.t;
    acadoVariables.x0[8] = state.s;
} else {
    fprintf(fp_actual, "%f %f %f %f %f %f %f\n",
      state.s,
      state.ey,
      state.epsi,
      sqrt(state.vx*state.vx+state.vy*state.vy),
      state.d,
      state.delta,
      state.t
    );

    for (i = 0; i < NX; ++i) {
      acadoVariables.x0[i] = acadoVariables.x[NX*next_n + i];
    }
}

    for (i = 0; i < N+1; ++i)
      acadoVariables.od[i] = get_curature(track.sp, acadoVariables.x0[8] + i*0.06);

    for (i = 0; i < next_n; i++) {
      acado_shiftStates(2, 0, 0);
      acado_shiftControls(0);
    }
    for (i = N + 1 - next_n; i < N + 1; ++i)  {
      for (j = 0; j < NX; ++j)
        acadoVariables.x[ i*NX + j ] = 0.0;
      acadoVariables.x[ i*NX + 2 ] = 1.0;
      acadoVariables.x[ i*NX + 5 ] = 1.0;
      acadoVariables.x[ i*NX + 6 ] = atan(acadoVariables.od[i]*0.06);
      if (acadoVariables.od[i]>1.0)
        acadoVariables.x[ i*NX + 0 ] = +0.1;
      else if (acadoVariables.od[i]<-1.0)
        acadoVariables.x[ i*NX + j ] = -0.1;

      for (j = 0; j < NU; ++j)
        acadoVariables.u[ (i-1)*NU + j ] = 0.0;
    }

    //for (i = 0; i < NU * N; ++i)  acadoVariables.u[ i ] = 0.0;

	}

  fclose(fp);
  fclose(fp_actual);

  return 0;
}
