
#include "dense.hpp"

CySolverDense::CySolverDense(double t_old, double t_now, double* y_in_ptr, unsigned int num_y) :
        t_old(t_old),
        t_now(t_now),
        num_y(num_y)
{
    // We need to copy over the values of y_now because they will be used later on whenever this interpolator is called.
    std::memcpy(this->y_stored_ptr, y_in_ptr, sizeof(double) * this->num_y);
}

void CySolverDense::call(double t_interp, double* y_interped)
{
    // Perform any checks and call implementation function.
    this->call_implementation(t_interp, y_interped);
}

void  CySolverDense::call_implementation(double t_interp, double* y_interped)
{
    // Base class just returns current y
    std::memcpy(y_interped, this->y_stored_ptr, sizeof(double) * this->num_y);
}
