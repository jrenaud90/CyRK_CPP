#pragma once

#include <cstring>

#include "common.hpp"

class CySolverDense
{
/* Attributes */
protected:


public:

    // Integrator info
    int integrator_int = -1;

    // y and t state info
    unsigned int num_y = 0;
    double t_old = 0.0;
    double t_now = 0.0;

    // Stored y values at interpolation step
    double y_stored[Y_LIMIT] = { };
    double* y_stored_ptr     = &y_stored[0];

/* Methods */
protected:

public:
    virtual ~CySolverDense() {};
    CySolverDense() {};
    CySolverDense(int integrator_int, double t_old, double t_now, double* y_in_ptr, unsigned int num_y);

    virtual void call(double t_interp, double* y_interped);
};
