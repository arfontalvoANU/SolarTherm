#ifndef ST_RECEIVER_STARTUP_TIME_H
#define ST_RECEIVER_STARTUP_TIME_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

double energy(const double* dni, const double* eta_op, double A_h, double dt, int nsteps, double Q_min);

#endif
