#ifndef NMPC4RacingCar_TYPES_H
#define NMPC4RacingCar_TYPES_H

namespace nmpc4rc {

  struct State {
    double ey;
    double epsi;
    double vx;
    double vy;
    double w;
    double d;
    double delta;
    double t;
    double s;

    double X; // used for Spline
    double Y; // used for Spline

    double dd;
    double ddelta;
    double curvature;
  };

}

#endif //NMPC4RacingCar_TYPES_H
