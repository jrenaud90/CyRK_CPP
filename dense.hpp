#pragma once

#include <cstring>

#include "common.hpp"

class CySolverDense
{
/* Attributes */
protected:


public:
    double t_old = 0.0;
    double t_now = 0.0;
    unsigned int num_y = 0;

    // Stored y values at interpolation step
    double y_stored[Y_LIMIT] = {};
    double* y_stored_ptr     = &y_stored[0];

/* Methods */
protected:

    virtual void call_implementation(double t_interp, double* y_interped);

public:
    virtual ~CySolverDense() {};
    CySolverDense() {};
    CySolverDense(double t_old, double t_now, double* y_in_ptr, unsigned int num_y);

    void call(double t_interp, double* y_interped);
};