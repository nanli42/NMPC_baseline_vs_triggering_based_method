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


#ifndef ACADO_COMMON_H
#define ACADO_COMMON_H

#include <math.h>
#include <string.h>

#ifndef __MATLAB__
#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */
#endif /* __MATLAB__ */

/** \defgroup ACADO ACADO CGT generated module. */
/** @{ */

/** qpOASES QP solver indicator. */
#define ACADO_QPOASES  0
#define ACADO_QPOASES3 1
/** FORCES QP solver indicator.*/
#define ACADO_FORCES   2
/** qpDUNES QP solver indicator.*/
#define ACADO_QPDUNES  3
/** HPMPC QP solver indicator. */
#define ACADO_HPMPC    4
#define ACADO_GENERIC    5

/** Indicator for determining the QP solver used by the ACADO solver code. */
#define ACADO_QP_SOLVER ACADO_QPOASES

#include "acado_qpoases_interface.hpp"


/*
 * Common definitions
 */
/** User defined block based condensing. */
#define ACADO_BLOCK_CONDENSING 0
/** Compute covariance matrix of the last state estimate. */
#define ACADO_COMPUTE_COVARIANCE_MATRIX 0
/** Flag indicating whether constraint values are hard-coded or not. */
#define ACADO_HARDCODED_CONSTRAINT_VALUES 0
/** Indicator for fixed initial state. */
#define ACADO_INITIAL_STATE_FIXED 1
/** Number of control/estimation intervals. */
#define ACADO_N 15
/** Number of online data values. */
#define ACADO_NOD 1
/** Number of path constraints. */
#define ACADO_NPAC 5
/** Number of control variables. */
#define ACADO_NU 2
/** Number of differential variables. */
#define ACADO_NX 9
/** Number of algebraic variables. */
#define ACADO_NXA 0
/** Number of differential derivative variables. */
#define ACADO_NXD 0
/** Number of references/measurements per node on the first N nodes. */
#define ACADO_NY 0
/** Number of references/measurements on the last (N + 1)st node. */
#define ACADO_NYN 0
/** Total number of QP optimization variables. */
#define ACADO_QP_NV 30
/** Number of integration steps per shooting interval. */
#define ACADO_RK_NIS 7
/** Number of Runge-Kutta stages per integration step. */
#define ACADO_RK_NSTAGES 4
/** Single versus double precision data type representation. */
#define ACADO_SINGLE_PRECISION 0
/** Providing interface for arrival cost. */
#define ACADO_USE_ARRIVAL_COST 0
/** Indicator for usage of non-hard-coded linear terms in the objective. */
#define ACADO_USE_LINEAR_TERMS 0
/** Indicator for type of fixed weighting matrices. */
#define ACADO_WEIGHTING_MATRICES_TYPE 2


/*
 * Globally used structure definitions
 */

/** The structure containing the user data.
 * 
 *  Via this structure the user "communicates" with the solver code.
 */
typedef struct ACADOvariables_
{
int dummy;
/** Matrix of size: 16 x 9 (row major format)
 * 
 *  Matrix containing 16 differential variable vectors.
 */
real_t x[ 144 ];

/** Matrix of size: 15 x 2 (row major format)
 * 
 *  Matrix containing 15 control variable vectors.
 */
real_t u[ 30 ];

/** Column vector of size: 16
 * 
 *  Matrix containing 16 online data vectors.
 */
real_t od[ 16 ];

/** Column vector of size: 9
 * 
 *  Current state feedback vector.
 */
real_t x0[ 9 ];

/** Column vector of size: 30
 * 
 *  Lower bounds values.
 */
real_t lbValues[ 30 ];

/** Column vector of size: 30
 * 
 *  Upper bounds values.
 */
real_t ubValues[ 30 ];

/** Column vector of size: 195
 * 
 *  Lower bounds values for affine constraints.
 */
real_t lbAValues[ 195 ];

/** Column vector of size: 195
 * 
 *  Upper bounds values for affine constraints.
 */
real_t ubAValues[ 195 ];

/** Matrix of size: 15 x 9 (row major format) */
real_t mu[ 135 ];


} ACADOvariables;

/** Private workspace used by the auto-generated code.
 * 
 *  Data members of this structure are private to the solver.
 *  In other words, the user code should not modify values of this 
 *  structure. 
 */
typedef struct ACADOworkspace_
{
/** Column vector of size: 2004 */
real_t rhs_aux[ 2004 ];

real_t rk_ttt;

/** Row vector of size: 186 */
real_t rk_xxx[ 186 ];

/** Matrix of size: 4 x 165 (row major format) */
real_t rk_kkk[ 660 ];

/** Row vector of size: 252 */
real_t rk_sweep1[ 252 ];

/** Row vector of size: 252 */
real_t rk_sweep2[ 252 ];

/** Row vector of size: 186 */
real_t state[ 186 ];

/** Column vector of size: 135 */
real_t d[ 135 ];

/** Matrix of size: 135 x 9 (row major format) */
real_t evGx[ 1215 ];

/** Matrix of size: 135 x 2 (row major format) */
real_t evGu[ 270 ];

/** Matrix of size: 165 x 11 (row major format) */
real_t EH[ 1815 ];

/** Matrix of size: 9 x 9 (row major format) */
real_t EH_N[ 81 ];

/** Column vector of size: 54 */
real_t objAuxVar[ 54 ];

/** Row vector of size: 12 */
real_t objValueIn[ 12 ];

/** Row vector of size: 91 */
real_t objValueOut[ 91 ];

/** Matrix of size: 135 x 9 (row major format) */
real_t Q1[ 1215 ];

/** Matrix of size: 30 x 2 (row major format) */
real_t R1[ 60 ];

/** Matrix of size: 135 x 2 (row major format) */
real_t S1[ 270 ];

/** Matrix of size: 9 x 9 (row major format) */
real_t QN1[ 81 ];

/** Column vector of size: 259 */
real_t conAuxVar[ 259 ];

/** Row vector of size: 17 */
real_t conValueIn[ 17 ];

/** Row vector of size: 181 */
real_t conValueOut[ 181 ];

/** Column vector of size: 75 */
real_t evH[ 75 ];

/** Matrix of size: 75 x 9 (row major format) */
real_t evHx[ 675 ];

/** Matrix of size: 75 x 2 (row major format) */
real_t evHu[ 150 ];

/** Matrix of size: 11 x 11 (row major format) */
real_t evDDH[ 121 ];

/** Column vector of size: 5 */
real_t evHxd[ 5 ];

/** Column vector of size: 144 */
real_t sbar[ 144 ];

/** Column vector of size: 9 */
real_t Dx0[ 9 ];

/** Matrix of size: 135 x 9 (row major format) */
real_t C[ 1215 ];

/** Matrix of size: 9 x 2 (row major format) */
real_t W1[ 18 ];

/** Matrix of size: 9 x 2 (row major format) */
real_t W2[ 18 ];

/** Matrix of size: 1080 x 2 (row major format) */
real_t E[ 2160 ];

/** Column vector of size: 144 */
real_t QDy[ 144 ];

/** Column vector of size: 9 */
real_t w1[ 9 ];

/** Column vector of size: 9 */
real_t w2[ 9 ];

/** Matrix of size: 30 x 30 (row major format) */
real_t H[ 900 ];

/** Matrix of size: 195 x 30 (row major format) */
real_t A[ 5850 ];

/** Column vector of size: 30 */
real_t g[ 30 ];

/** Column vector of size: 30 */
real_t lb[ 30 ];

/** Column vector of size: 30 */
real_t ub[ 30 ];

/** Column vector of size: 195 */
real_t lbA[ 195 ];

/** Column vector of size: 195 */
real_t ubA[ 195 ];

/** Column vector of size: 30 */
real_t x[ 30 ];

/** Column vector of size: 225 */
real_t y[ 225 ];


} ACADOworkspace;

/* 
 * Forward function declarations. 
 */


/** Performs the integration and sensitivity propagation for one shooting interval.
 *
 *  \param rk_eta Working array to pass the input values and return the results.
 *  \param resetIntegrator The internal memory of the integrator can be reset.
 *
 *  \return Status code of the integrator.
 */
int acado_integrate( real_t* const rk_eta, int resetIntegrator );

/** Export of an ACADO symbolic function.
 *
 *  \param in Input to the exported function.
 *  \param out Output of the exported function.
 */
void acado_acado_forward(const real_t* in, real_t* out);

/** Export of an ACADO symbolic function.
 *
 *  \param in Input to the exported function.
 *  \param out Output of the exported function.
 */
void acado_acado_backward(const real_t* in, real_t* out);

/** Preparation step of the RTI scheme.
 *
 *  \return Status of the integration module. =0: OK, otherwise the error code.
 */
int acado_preparationStep(  );

/** Feedback/estimation step of the RTI scheme.
 *
 *  \return Status code of the qpOASES QP solver.
 */
int acado_feedbackStep(  );

/** Solver initialization. Must be called once before any other function call.
 *
 *  \return =0: OK, otherwise an error code of a QP solver.
 */
int acado_initializeSolver(  );

/** Initialize shooting nodes by a forward simulation starting from the first node.
 */
void acado_initializeNodesByForwardSimulation(  );

/** Shift differential variables vector by one interval.
 *
 *  \param strategy Shifting strategy: 1. Initialize node 16 with xEnd. 2. Initialize node 16 by forward simulation.
 *  \param xEnd Value for the x vector on the last node. If =0 the old value is used.
 *  \param uEnd Value for the u vector on the second to last node. If =0 the old value is used.
 */
void acado_shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd );

/** Shift controls vector by one interval.
 *
 *  \param uEnd Value for the u vector on the second to last node. If =0 the old value is used.
 */
void acado_shiftControls( real_t* const uEnd );

/** Get the KKT tolerance of the current iterate.
 *
 *  \return The KKT tolerance value.
 */
real_t acado_getKKT(  );

/** Calculate the objective value.
 *
 *  \return Value of the objective function.
 */
real_t acado_getObjective(  );

/** EVD-based regularization of a Hessian block.
 *
 */
void acado_regularize( real_t* const hessian_block );


/* 
 * Extern declarations. 
 */

extern ACADOworkspace acadoWorkspace;
extern ACADOvariables acadoVariables;

/** @} */

#ifndef __MATLAB__
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */
#endif /* __MATLAB__ */

#endif /* ACADO_COMMON_H */
