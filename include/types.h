#ifndef NMPC4RacingCar_TYPES_H
#define NMPC4RacingCar_TYPES_H

namespace nmpc4rc {

  struct State {
    double ey;
    double epsi;
    double v;
    double beta;
    double omega;
    double delta;
    double t;
    double s;

    double a;
    double vd;

    double ks;

    double X;
    double Y;

    double vx;
    double vy;
    double w;
    double d;
    double dd;
    double ddelta;
    double curvature;
  };

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

}

#endif //NMPC4RacingCar_TYPES_H
