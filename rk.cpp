#include "rk.hpp"

// ########################################################################################################################
// RKSolver (Base)
// ########################################################################################################################
// Constructors
RKSolver::RKSolver() {}
RKSolver::RKSolver(
    // Base Class input arguments
    DiffeqFuncType diffeq_ptr,
    std::shared_ptr<CySolverResult> const storage_ptr,
    const double t_start,
    const double t_end,
    const double* y0_ptr,
    const unsigned int num_y,
    const unsigned int num_extra,
    const double* args_ptr,
    const size_t max_num_steps,
    const size_t max_ram_MB,
    const bool use_dense_output,
    const double* t_eval,
    const size_t len_t_eval,
    // RKSolver input arguments
    const double rtol,
    const double atol,
    const double* rtols_ptr,
    const double* atols_ptr,
    const double max_step_size,
    const double first_step_size) :
        CySolverBase(
            diffeq_ptr,
            storage_ptr,
            t_start,
            t_end,
            y0_ptr,
            num_y,
            num_extra,
            args_ptr,
            max_num_steps,
            max_ram_MB,
            use_dense_output,
            t_eval,
            len_t_eval),
        max_step_size(max_step_size),
        user_provided_first_step_size(first_step_size)
        
{
    // Check for errors
    if (this->user_provided_first_step_size != 0.0)
    {
        if (this->user_provided_first_step_size < 0.0)
        {
            this->storage_ptr->error_code = -1;
            this->storage_ptr->update_message("User-provided initial step size must be a positive number.");
        }
        else if (first_step_size > (this->t_delta_abs * 0.5))
        {
            this->storage_ptr->error_code = -1;
            this->storage_ptr->update_message("User-provided initial step size must be smaller than 50 % of the time span size.");
        }
    }

    // Setup tolerances
    // User can provide an array of relative tolerances, one for each y value.
    // The length of the pointer array must be the same as y0 (and <= 25).

    if (rtols_ptr)
    {
        // rtol for each y
        this->use_array_rtols = true;
        for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
        {
            double temp_double = rtols_ptr[y_i];
            if (temp_double < EPS_100)
            {
                temp_double = EPS_100;
            }
            this->rtols_ptr[y_i] = temp_double;
        }
    }
    else {
        // only one rtol
        double temp_double = rtol;
        if (temp_double < EPS_100)
        {
            temp_double = EPS_100;
        }
        this->rtols_ptr[0] = temp_double;
    }

    if (atols_ptr)
    {
        // atol for each y
        this->use_array_atols = true;
        std::memcpy(this->atols_ptr, atols_ptr, sizeof(double) * this->num_y);
    }
    else {
        // only one atol
        this->atols_ptr[0] = atol;
    }
}


// Destructors
RKSolver::~RKSolver()
{

}


// Protected Methods
void RKSolver::p_estimate_error()
{   
    // Initialize rtol and atol
    double rtol = this->rtols_ptr[0];
    double atol = this->atols_ptr[0];

    this->error_norm = 0.0;

    for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
    {
        if (this->use_array_rtols)
        {
            rtol = this->rtols_ptr[y_i];
        }

        if (this->use_array_atols)
        {
            atol = this->atols_ptr[y_i];
        }

        // Dot product between K and E
        double error_dot;
        const unsigned int stride_K = y_i * this->n_stages_p1;

        switch (this->n_stages)
        {
        // These loops go 1 more than `n_stages`.
        // Note: DOP853 is handled in an override by its subclass.
        case(3):
            // RK23
            error_dot =  this->E_ptr[0] * this->K_ptr[stride_K];
            error_dot += this->E_ptr[1] * this->K_ptr[stride_K + 1];
            error_dot += this->E_ptr[2] * this->K_ptr[stride_K + 2];
            error_dot += this->E_ptr[3] * this->K_ptr[stride_K + 3];

            break;
        case(6):
            // RK45
            error_dot =  this->E_ptr[0] * this->K_ptr[stride_K];
            error_dot += this->E_ptr[1] * this->K_ptr[stride_K + 1];
            error_dot += this->E_ptr[2] * this->K_ptr[stride_K + 2];
            error_dot += this->E_ptr[3] * this->K_ptr[stride_K + 3];
            error_dot += this->E_ptr[4] * this->K_ptr[stride_K + 4];
            error_dot += this->E_ptr[5] * this->K_ptr[stride_K + 5];
            error_dot += this->E_ptr[6] * this->K_ptr[stride_K + 6];

            break;
        default:
            // Resort to unrolled loops
            // Initialize
            error_dot = 0.0;
            // New or Non-optimized RK method. default to for loop.
            for (unsigned int j = 0; j < (this->n_stages + 1); j++)
            {
                error_dot += this->E_ptr[j] * this->K_ptr[stride_K + j];
            }
            break;
        }

        // Find scale of y for error calculations
        const double scale_inv = 1.0 / (atol + std::fmax(std::fabs(this->y_old_ptr[y_i]), std::fabs(this->y_now_ptr[y_i])) * rtol);

        // We need the absolute value but since we are taking the square, it is guaranteed to be positive.
        // TODO: This will need to change if CySolver ever accepts complex numbers
        // error_norm_abs = fabs(error_dot_1)
        error_dot *= scale_inv;

        this->error_norm += (error_dot * error_dot);
    }
    this->error_norm = this->step_size * std::sqrt(this->error_norm) / this->num_y_sqrt;
}


void RKSolver::p_step_implementation()
{
    // Run RK integration step
    
    // Create local variables instead of calling class attributes for pointer objects.
    double* l_K_ptr       = this->K_ptr;
    const double* l_A_ptr = this->A_ptr;
    const double* l_B_ptr = this->B_ptr;
    const double* l_C_ptr = this->C_ptr;
    double* l_y_now_ptr   = this->y_now_ptr;
    double* l_y_old_ptr   = this->y_old_ptr;
    double* l_dy_now_ptr  = this->dy_now_ptr;
    double* l_dy_old_ptr  = this->dy_old_ptr;

    // Other local variables
    double l_t_old     = this->t_old;
    double l_step      = this->step;
    double l_step_size = this->step_size;

    // Determine step size based on previous loop
    // Find minimum step size based on the value of t (less floating point numbers between numbers when t is large)
    const double min_step_size = 10. * std::fabs(std::nextafter(l_t_old, this->direction_inf) - l_t_old);
    // Look for over/undershoots in previous step size
    l_step_size = std::min<double>(l_step_size, this->max_step_size);
    l_step_size = std::max<double>(l_step_size, min_step_size);

    // Optimization variables
    // Define a very specific A (Row 1; Col 0) now since it is called consistently and does not change.
    const double A_at_10 = l_A_ptr[1 * this->len_Acols + 0];

    // Determine new step size
    bool step_accepted = false;
    bool step_rejected = false;
    bool step_error    = false;

    // !! Step Loop
    while (!step_accepted) {

        // Check if step size is too small
        // This will cause integration to fail: step size smaller than spacing between numbers
        if (l_step_size < min_step_size) {
            step_error = true;
            this->status = -1;
            break;
        }

        // Move time forward for this particular step size
        double t_delta_check;
        if (this->direction_flag) {
            l_step    = l_step_size;
            this->t_now_ptr[0]   = l_t_old + l_step;
            t_delta_check = this->t_now_ptr[0] - this->t_end;
        }
        else {
            l_step    = -l_step_size;
            this->t_now_ptr[0]   = l_t_old + l_step;
            t_delta_check = this->t_end - this->t_now_ptr[0];
        }

        // Check that we are not at the end of integration with that move
        if (t_delta_check > 0.0) {
            this->t_now_ptr[0] = this->t_end;

            // If we are, correct the step so that it just hits the end of integration.
            l_step = this->t_now_ptr[0] - l_t_old;
            if (this->direction_flag) {
                l_step_size = l_step;
            }
            else {
                l_step_size = -l_step;
            }
        }

        // !! Calculate derivative using RK method

        // t_now must be updated for each loop of s in order to make the diffeq method calls.
        // But we need to return to its original value later on. Store in temp variable.
        const double time_tmp = this->t_now_ptr[0];

        for (unsigned int s = 1; s < this->len_C; s++) {
            // Find the current time based on the old time and the step size.
            this->t_now_ptr[0] = l_t_old + l_C_ptr[s] * l_step;
            const unsigned int stride_A = s * this->len_Acols;

            for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
            {
                double temp_double;
                const unsigned int stride_K = y_i * this->n_stages_p1;
                // Dot Product (K, a) * step
                switch (s)
                {
                case(1):
                    // Set the first column of K
                    temp_double = l_dy_old_ptr[y_i];
                    // K[0, :] == first part of the array
                    l_K_ptr[stride_K] = temp_double;
                    temp_double *= A_at_10;
                    break;
                case(2):
                    // Loop through (j = 0; j < s; j++)
                    // j = 0
                    temp_double  = l_A_ptr[stride_A]     * l_K_ptr[stride_K];
                    // j = 1
                    temp_double += l_A_ptr[stride_A + 1] * l_K_ptr[stride_K + 1];
                    break;
                case(3):
                    // Loop through (j = 0; j < s; j++)
                    // j = 0
                    temp_double  = l_A_ptr[stride_A]     * l_K_ptr[stride_K];
                    // j = 1
                    temp_double += l_A_ptr[stride_A + 1] * l_K_ptr[stride_K + 1];
                    // j = 2
                    temp_double += l_A_ptr[stride_A + 2] * l_K_ptr[stride_K + 2];
                    break;
                case(4):
                    // Loop through (j = 0; j < s; j++)
                    // j = 0
                    temp_double  = l_A_ptr[stride_A]     * l_K_ptr[stride_K];
                    // j = 1
                    temp_double += l_A_ptr[stride_A + 1] * l_K_ptr[stride_K + 1];
                    // j = 2
                    temp_double += l_A_ptr[stride_A + 2] * l_K_ptr[stride_K + 2];
                    // j = 3
                    temp_double += l_A_ptr[stride_A + 3] * l_K_ptr[stride_K + 3];
                    break;
                case(5):
                    // Loop through (j = 0; j < s; j++)
                    // j = 0
                    temp_double  = l_A_ptr[stride_A]     * l_K_ptr[stride_K];
                    // j = 1
                    temp_double += l_A_ptr[stride_A + 1] * l_K_ptr[stride_K + 1];
                    // j = 2
                    temp_double += l_A_ptr[stride_A + 2] * l_K_ptr[stride_K + 2];
                    // j = 3
                    temp_double += l_A_ptr[stride_A + 3] * l_K_ptr[stride_K + 3];
                    // j = 4
                    temp_double += l_A_ptr[stride_A + 4] * l_K_ptr[stride_K + 4];
                    break;
                case(6):
                    // Loop through (j = 0; j < s; j++)
                    // j = 0
                    temp_double  = l_A_ptr[stride_A]     * l_K_ptr[stride_K];
                    // j = 1
                    temp_double += l_A_ptr[stride_A + 1] * l_K_ptr[stride_K + 1];
                    // j = 2
                    temp_double += l_A_ptr[stride_A + 2] * l_K_ptr[stride_K + 2];
                    // j = 3
                    temp_double += l_A_ptr[stride_A + 3] * l_K_ptr[stride_K + 3];
                    // j = 4
                    temp_double += l_A_ptr[stride_A + 4] * l_K_ptr[stride_K + 4];
                    // j = 5
                    temp_double += l_A_ptr[stride_A + 5] * l_K_ptr[stride_K + 5];
                    break;
                case(7):
                    // Loop through (j = 0; j < s; j++)
                    // j = 0
                    temp_double  = l_A_ptr[stride_A]     * l_K_ptr[stride_K];
                    // j = 1
                    temp_double += l_A_ptr[stride_A + 1] * l_K_ptr[stride_K + 1];
                    // j = 2
                    temp_double += l_A_ptr[stride_A + 2] * l_K_ptr[stride_K + 2];
                    // j = 3
                    temp_double += l_A_ptr[stride_A + 3] * l_K_ptr[stride_K + 3];
                    // j = 4
                    temp_double += l_A_ptr[stride_A + 4] * l_K_ptr[stride_K + 4];
                    // j = 5
                    temp_double += l_A_ptr[stride_A + 5] * l_K_ptr[stride_K + 5];
                    // j = 6
                    temp_double += l_A_ptr[stride_A + 6] * l_K_ptr[stride_K + 6];
                    break;
                case(8):
                    // Loop through (j = 0; j < s; j++)
                    // j = 0
                    temp_double  = l_A_ptr[stride_A]     * l_K_ptr[stride_K];
                    // j = 1
                    temp_double += l_A_ptr[stride_A + 1] * l_K_ptr[stride_K + 1];
                    // j = 2
                    temp_double += l_A_ptr[stride_A + 2] * l_K_ptr[stride_K + 2];
                    // j = 3
                    temp_double += l_A_ptr[stride_A + 3] * l_K_ptr[stride_K + 3];
                    // j = 4
                    temp_double += l_A_ptr[stride_A + 4] * l_K_ptr[stride_K + 4];
                    // j = 5
                    temp_double += l_A_ptr[stride_A + 5] * l_K_ptr[stride_K + 5];
                    // j = 6
                    temp_double += l_A_ptr[stride_A + 6] * l_K_ptr[stride_K + 6];
                    // j = 7
                    temp_double += l_A_ptr[stride_A + 7] * l_K_ptr[stride_K + 7];
                    break;
                case(9):
                    // Loop through (j = 0; j < s; j++)
                    // j = 0
                    temp_double  = l_A_ptr[stride_A]     * l_K_ptr[stride_K];
                    // j = 1
                    temp_double += l_A_ptr[stride_A + 1] * l_K_ptr[stride_K + 1];
                    // j = 2
                    temp_double += l_A_ptr[stride_A + 2] * l_K_ptr[stride_K + 2];
                    // j = 3
                    temp_double += l_A_ptr[stride_A + 3] * l_K_ptr[stride_K + 3];
                    // j = 4
                    temp_double += l_A_ptr[stride_A + 4] * l_K_ptr[stride_K + 4];
                    // j = 5
                    temp_double += l_A_ptr[stride_A + 5] * l_K_ptr[stride_K + 5];
                    // j = 6
                    temp_double += l_A_ptr[stride_A + 6] * l_K_ptr[stride_K + 6];
                    // j = 7
                    temp_double += l_A_ptr[stride_A + 7] * l_K_ptr[stride_K + 7];
                    // j = 8
                    temp_double += l_A_ptr[stride_A + 8] * l_K_ptr[stride_K + 8];
                    break;
                case(10):
                    // Loop through (j = 0; j < s; j++)
                    // j = 0
                    temp_double  = l_A_ptr[stride_A]     * l_K_ptr[stride_K];
                    // j = 1
                    temp_double += l_A_ptr[stride_A + 1] * l_K_ptr[stride_K + 1];
                    // j = 2
                    temp_double += l_A_ptr[stride_A + 2] * l_K_ptr[stride_K + 2];
                    // j = 3
                    temp_double += l_A_ptr[stride_A + 3] * l_K_ptr[stride_K + 3];
                    // j = 4
                    temp_double += l_A_ptr[stride_A + 4] * l_K_ptr[stride_K + 4];
                    // j = 5
                    temp_double += l_A_ptr[stride_A + 5] * l_K_ptr[stride_K + 5];
                    // j = 6
                    temp_double += l_A_ptr[stride_A + 6] * l_K_ptr[stride_K + 6];
                    // j = 7
                    temp_double += l_A_ptr[stride_A + 7] * l_K_ptr[stride_K + 7];
                    // j = 8
                    temp_double += l_A_ptr[stride_A + 8] * l_K_ptr[stride_K + 8];
                    // j = 9
                    temp_double += l_A_ptr[stride_A + 9] * l_K_ptr[stride_K + 9];
                    break;
                case(11):
                    // Loop through (j = 0; j < s; j++)
                    // j = 0
                    temp_double  = l_A_ptr[stride_A]      * l_K_ptr[stride_K];
                    // j = 1
                    temp_double += l_A_ptr[stride_A + 1]  * l_K_ptr[stride_K + 1];
                    // j = 2
                    temp_double += l_A_ptr[stride_A + 2]  * l_K_ptr[stride_K + 2];
                    // j = 3
                    temp_double += l_A_ptr[stride_A + 3]  * l_K_ptr[stride_K + 3];
                    // j = 4
                    temp_double += l_A_ptr[stride_A + 4]  * l_K_ptr[stride_K + 4];
                    // j = 5
                    temp_double += l_A_ptr[stride_A + 5]  * l_K_ptr[stride_K + 5];
                    // j = 6
                    temp_double += l_A_ptr[stride_A + 6]  * l_K_ptr[stride_K + 6];
                    // j = 7
                    temp_double += l_A_ptr[stride_A + 7]  * l_K_ptr[stride_K + 7];
                    // j = 8
                    temp_double += l_A_ptr[stride_A + 8]  * l_K_ptr[stride_K + 8];
                    // j = 9
                    temp_double += l_A_ptr[stride_A + 9]  * l_K_ptr[stride_K + 9];
                    // j = 10
                    temp_double += l_A_ptr[stride_A + 10] * l_K_ptr[stride_K + 10];
                    break;
                default:
                    // Resort to regular rolled loops
                    // Initialize
                    temp_double = 0.0;

                    for (unsigned int j = 0; j < s; j++)
                    {
                        temp_double += l_A_ptr[stride_A + j] * l_K_ptr[stride_K + j];
                    }
                    break;
                }
                // Update value of y_now
                l_y_now_ptr[y_i] = l_y_old_ptr[y_i] + (l_step * temp_double);
            }
            // Call diffeq method to update K with the new dydt
            // This will use the now updated dy_now_ptr based on the values of y_now_ptr and t_now_ptr.
            this->diffeq(this);

            // Update K based on the new dy values.
            for (unsigned int y_i = 0; y_i < this->num_y; y_i++) {
                const unsigned int stride_K = y_i * this->n_stages_p1;
                l_K_ptr[stride_K + s] = l_dy_now_ptr[y_i];
            }
        }

        // Restore t_now to its previous value.
        this->t_now_ptr[0] = time_tmp;

        // Dot Product (K, B) * step
        for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
        {
            double temp_double;
            const unsigned int stride_K = y_i * this->n_stages_p1;

            switch (this->n_stages)
            {
            case(3):
                // RK23
                temp_double  = l_B_ptr[0] * l_K_ptr[stride_K];
                temp_double += l_B_ptr[1] * l_K_ptr[stride_K + 1];
                temp_double += l_B_ptr[2] * l_K_ptr[stride_K + 2];
                break;
            case(6):
                //RK45
                temp_double  = l_B_ptr[0] * l_K_ptr[stride_K];
                temp_double += l_B_ptr[1] * l_K_ptr[stride_K + 1];
                temp_double += l_B_ptr[2] * l_K_ptr[stride_K + 2];
                temp_double += l_B_ptr[3] * l_K_ptr[stride_K + 3];
                temp_double += l_B_ptr[4] * l_K_ptr[stride_K + 4];
                temp_double += l_B_ptr[5] * l_K_ptr[stride_K + 5];
                break;
            case(12):
                //DOP853
                temp_double  = l_B_ptr[0]  * l_K_ptr[stride_K];
                temp_double += l_B_ptr[1]  * l_K_ptr[stride_K + 1];
                temp_double += l_B_ptr[2]  * l_K_ptr[stride_K + 2];
                temp_double += l_B_ptr[3]  * l_K_ptr[stride_K + 3];
                temp_double += l_B_ptr[4]  * l_K_ptr[stride_K + 4];
                temp_double += l_B_ptr[5]  * l_K_ptr[stride_K + 5];
                temp_double += l_B_ptr[6]  * l_K_ptr[stride_K + 6];
                temp_double += l_B_ptr[7]  * l_K_ptr[stride_K + 7];
                temp_double += l_B_ptr[8]  * l_K_ptr[stride_K + 8];
                temp_double += l_B_ptr[9]  * l_K_ptr[stride_K + 9];
                temp_double += l_B_ptr[10] * l_K_ptr[stride_K + 10];
                temp_double += l_B_ptr[11] * l_K_ptr[stride_K + 11];
                break;
            default:
                // Resort to rolled loops
                // Initialize
                temp_double = 0.0;

                for (unsigned int j = 0; j < this->n_stages; j++)
                {
                    temp_double += l_B_ptr[j] * l_K_ptr[stride_K + j];
                }
                break;
            }
            // Update y_now
            l_y_now_ptr[y_i] = l_y_old_ptr[y_i] + (l_step * temp_double);
        }

        // Find final dydt for this timestep
        // This will use the now updated dy_now_ptr based on the values of y_now_ptr and t_now_ptr.
        this->diffeq(this);

        // Set last column of K equal to dydt. K has size num_y * (n_stages + 1) so the last column is at n_stages
        for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
        {
            const unsigned int stride_K = y_i * this->n_stages_p1;
            l_K_ptr[stride_K + this->n_stages] = l_dy_now_ptr[y_i];
        }

        // Check how well this step performed by calculating its error.
        this->p_estimate_error();

        // Check the size of the error
        if (this->error_norm < 1.0) {
            // We found our step size because the error is low!
            // Update this step for the next time loop
            double step_factor = this->max_step_factor;
            if (this->error_norm != 0.0)
            {
                // Estimate a new step size based on the error.
                const double error_safe = this->error_safety / std::pow(this->error_norm, this->error_exponent);
                step_factor = std::min<double>(this->max_step_factor, error_safe);
            }

            if (step_rejected) {
                // There were problems with this step size on the previous step loop. Make sure factor does
                //   not exasperate them.
                step_factor = std::min<double>(step_factor, 1.0);
            }

            // Update step size
            l_step_size *= step_factor;
            step_accepted = true;
        }
        else {
            // Error is still large. Keep searching for a better step size.
            const double error_safe = this->error_safety / std::pow(this->error_norm, this->error_exponent);
            l_step_size *= std::max<double>(this->min_step_factor, error_safe);
            step_rejected = true;
        }
    }

    // Update status depending if there were any errors.
    if (step_error) {
        // Issue with step convergence
        this->status = -1;
    }
    else if (!step_accepted) {
        // Issue with step convergence
        this->status = -7;
    }

    // End of RK step. 
    // Update state variables
    this->step_size = l_step_size;
    this->step      = l_step;
}


// Public methods
void RKSolver::reset()
{
    // Update stride information
    this->nstages_numy = this->n_stages * this->num_y;
    this->n_stages_p1  = this->n_stages + 1;

    // It is important to initialize the K variable with zeros
    std::fill(this->K_ptr, this->K_ptr + (this->num_y * this->n_stages_p1), 0.0);

    // Call base class reset after K is established but before first step size is calculated.
    CySolverBase::reset();

    // Update initial step size
    if (this->user_provided_first_step_size == 0.0)
    {
        // User did not provide a step size. Try to find a good guess.
        this->calc_first_step_size();
    }
    else {
        this->step_size = this->user_provided_first_step_size;
    }
}

void RKSolver::calc_first_step_size()
{
    /*
        Select an initial step size based on the differential equation.
        .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
            Equations I: Nonstiff Problems", Sec. II.4.
    */

    if (this->num_y == 0)
    {
        this->step_size = INF;
    }
    else {
        double d0, d1, d2, d0_abs, d1_abs, d2_abs, scale;
        double h0, h0_direction, h1;
        double y_old_tmp;

        // Initialize tolerances to the 0 place. If `use_array_rtols` (or atols) is set then this will change in the loop.
        double rtol = this->rtols_ptr[0];
        double atol = this->atols_ptr[0];

        // Find the norm for d0 and d1
        d0 = 0.0;
        d1 = 0.0;
        for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
        {
            if (this->use_array_rtols)
            {
                rtol = this->rtols_ptr[y_i];
            }
            if (this->use_array_atols)
            {
                atol = this->atols_ptr[y_i];
            }

            y_old_tmp = this->y_old_ptr[y_i];
            scale     = atol + std::fabs(y_old_tmp) * rtol;
            d0_abs    = std::fabs(y_old_tmp / scale);
            d1_abs    = std::fabs(dy_old_ptr[y_i] / scale);
            d0       += (d0_abs * d0_abs);
            d1       += (d1_abs * d1_abs);
        }

        d0 = std::sqrt(d0) / this->num_y_sqrt;
        d1 = std::sqrt(d1) / this->num_y_sqrt;

        if ((d0 < 1.0e-5) || (d1 < 1.0e-5))
        {
            h0 = 1.0e-6;
        }
        else {
            h0 = 0.01 * d0 / d1;
        }

        h0_direction = this->direction_flag ? h0 : -h0;

        this->t_now_ptr[0] = this->t_old + h0_direction;
        for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
        {
            this->y_now_ptr[y_i] = this->y_old_ptr[y_i] + h0_direction * this->dy_old_ptr[y_i];
        }

        // Update dy
        this->diffeq(this);

        // Find the norm for d2
        d2 = 0.0;
        for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
        {
            if (this->use_array_rtols)
            {
                rtol = this->rtols_ptr[y_i];
            }
            if (this->use_array_atols)
            {
                atol = this->atols_ptr[y_i];
            }

            scale  = atol + std::fabs(this->y_old_ptr[y_i]) * rtol;
            d2_abs = std::fabs((this->dy_now_ptr[y_i] - this->dy_old_ptr[y_i]) / scale);
            d2    += (d2_abs * d2_abs);
        }

        d2 = std::sqrt(d2) / (h0 * this->num_y_sqrt);

        if ((d1 <= 1.0e-15) && (d2 <= 1.0e-15))
        {
            h1 = std::fmax(1.0e-6, h0 * 1.0e-3);
        }
        else {
            h1 = std::pow((0.01 / std::fmax(d1, d2)), this->error_exponent);
        }
        this->step_size = std::fmax(10. * std::fabs(std::nextafter(this->t_old, this->direction_inf) - this->t_old), std::fmin(100.0 * h0, h1));
    }
}

/* Dense Output Methods */

void RKSolver::p_update_Q(double* Q_ptr) const
{
    // Q's definition depends on the integrators implementation. 
    // For default RK, it is defined by Q = K.T.dot(self.P)  K has shape of (n_stages + 1, num_y) so K.T has shape of (num_y, n_stages + 1)
    // P has shape of (4, 3) for RK23; (7, 4) for RK45.. So (n_stages + 1, Q_order)
    // So Q has shape of (num_y, num_Pcols)

    for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
    {
        const unsigned int stride_Q = y_i * this->len_Pcols;
        const unsigned int stride_K = y_i * this->n_stages_p1;

        switch (this->method_int)
        {
            double temp_double;
            unsigned int stride_P;

        case(0):
            // RK23
            // len_Pcols == 3; n_stages + 1 == 4

            // P = 0
            temp_double = this->K_ptr[stride_K] * this->P_ptr[0];
            temp_double += this->K_ptr[stride_K + 1] * this->P_ptr[1];
            temp_double += this->K_ptr[stride_K + 2] * this->P_ptr[2];
            temp_double += this->K_ptr[stride_K + 3] * this->P_ptr[3];
            Q_ptr[stride_Q] = temp_double;

            // P = 1
            stride_P = this->n_stages_p1;
            temp_double = this->K_ptr[stride_K] * this->P_ptr[stride_P];
            temp_double += this->K_ptr[stride_K + 1] * this->P_ptr[stride_P + 1];
            temp_double += this->K_ptr[stride_K + 2] * this->P_ptr[stride_P + 2];
            temp_double += this->K_ptr[stride_K + 3] * this->P_ptr[stride_P + 3];
            Q_ptr[stride_Q + 1] = temp_double;

            // P = 2
            stride_P += this->n_stages_p1;
            temp_double = this->K_ptr[stride_K] * this->P_ptr[stride_P];
            temp_double += this->K_ptr[stride_K + 1] * this->P_ptr[stride_P + 1];
            temp_double += this->K_ptr[stride_K + 2] * this->P_ptr[stride_P + 2];
            temp_double += this->K_ptr[stride_K + 3] * this->P_ptr[stride_P + 3];
            Q_ptr[stride_Q + 2] = temp_double;

            break;
        case(1):
            // RK45
            // len_Pcols == 4; n_stages + 1 == 7
            // P = 0
            temp_double = this->K_ptr[stride_K] * this->P_ptr[0];
            temp_double += this->K_ptr[stride_K + 1] * this->P_ptr[1];
            temp_double += this->K_ptr[stride_K + 2] * this->P_ptr[2];
            temp_double += this->K_ptr[stride_K + 3] * this->P_ptr[3];
            temp_double += this->K_ptr[stride_K + 4] * this->P_ptr[4];
            temp_double += this->K_ptr[stride_K + 5] * this->P_ptr[5];
            temp_double += this->K_ptr[stride_K + 6] * this->P_ptr[6];
            Q_ptr[stride_Q] = temp_double;

            // P = 1
            stride_P = this->n_stages_p1;
            temp_double = this->K_ptr[stride_K] * this->P_ptr[stride_P];
            temp_double += this->K_ptr[stride_K + 1] * this->P_ptr[stride_P + 1];
            temp_double += this->K_ptr[stride_K + 2] * this->P_ptr[stride_P + 2];
            temp_double += this->K_ptr[stride_K + 3] * this->P_ptr[stride_P + 3];
            temp_double += this->K_ptr[stride_K + 4] * this->P_ptr[stride_P + 4];
            temp_double += this->K_ptr[stride_K + 5] * this->P_ptr[stride_P + 5];
            temp_double += this->K_ptr[stride_K + 6] * this->P_ptr[stride_P + 6];
            Q_ptr[stride_Q + 1] = temp_double;

            // P = 2
            stride_P += this->n_stages_p1;
            temp_double = this->K_ptr[stride_K] * this->P_ptr[stride_P];
            temp_double += this->K_ptr[stride_K + 1] * this->P_ptr[stride_P + 1];
            temp_double += this->K_ptr[stride_K + 2] * this->P_ptr[stride_P + 2];
            temp_double += this->K_ptr[stride_K + 3] * this->P_ptr[stride_P + 3];
            temp_double += this->K_ptr[stride_K + 4] * this->P_ptr[stride_P + 4];
            temp_double += this->K_ptr[stride_K + 5] * this->P_ptr[stride_P + 5];
            temp_double += this->K_ptr[stride_K + 6] * this->P_ptr[stride_P + 6];
            Q_ptr[stride_Q + 2] = temp_double;

            // P = 3
            stride_P += this->n_stages_p1;
            temp_double = this->K_ptr[stride_K] * this->P_ptr[stride_P];
            temp_double += this->K_ptr[stride_K + 1] * this->P_ptr[stride_P + 1];
            temp_double += this->K_ptr[stride_K + 2] * this->P_ptr[stride_P + 2];
            temp_double += this->K_ptr[stride_K + 3] * this->P_ptr[stride_P + 3];
            temp_double += this->K_ptr[stride_K + 4] * this->P_ptr[stride_P + 4];
            temp_double += this->K_ptr[stride_K + 5] * this->P_ptr[stride_P + 5];
            temp_double += this->K_ptr[stride_K + 6] * this->P_ptr[stride_P + 6];
            Q_ptr[stride_Q + 3] = temp_double;

            break;
        case(2):
            // DOP853
            // TODO
            break;
        default:
            for (unsigned int P_i = 0; P_i < this->len_Pcols; P_i++)
            {
                const unsigned int stride_P = P_i * this->n_stages_p1;

                // Initialize dot product
                double temp_double = 0.0;


                for (unsigned int n_i = 0; n_i < this->n_stages_p1; n_i++)
                {
                    temp_double += this->K_ptr[stride_K + n_i] * this->P_ptr[stride_P + n_i];
                }

                // Set equal to Q
                Q_ptr[stride_Q + P_i] = temp_double;
            }
            break;
        }
    }
}

CySolverDense* RKSolver::p_dense_output_heap()
{
    // Build dense output object instance
    RKDenseOutput* dense_output = new RKDenseOutput(this->t_old, this->t_now_ptr[0], this->y_old_ptr, this->num_y, this->len_Pcols);

    // Update Q
    this->p_update_Q(dense_output->Q_ptr);

    return dense_output;
}

CySolverDense RKSolver::p_dense_output_stack()
{
    // Build dense output object instance
    RKDenseOutput dense_output = RKDenseOutput(this->t_old, this->t_now_ptr[0], this->y_old_ptr, this->num_y, this->len_Pcols);

    // Update Q
    this->p_update_Q(dense_output.Q_ptr);

    return dense_output;
}


// ########################################################################################################################
// Explicit Runge - Kutta 2(3)
// ########################################################################################################################
void RK23::reset()
{
    // Setup RK constants before calling the base class reset
    this->C_ptr     = RK23_C_ptr;
    this->A_ptr     = RK23_A_ptr;
    this->B_ptr     = RK23_B_ptr;
    this->E_ptr     = RK23_E_ptr;
    this->E3_ptr    = nullptr;       // Not used for RK23
    this->E5_ptr    = nullptr;       // Not used for RK23
    this->P_ptr     = RK23_P_ptr;       
    this->D_ptr     = nullptr;       // Not used for RK23
    this->K_ptr     = &this->K[0];
    this->order     = RK23_order;
    this->n_stages  = RK23_n_stages;
    this->len_Acols = RK23_len_Acols;
    this->len_C     = RK23_len_C;
    this->len_Pcols = RK23_len_Pcols;
    this->error_estimator_order = RK23_error_estimator_order;
    this->error_exponent = RK23_error_exponent;
    this->method_int = RK23_METHOD_INT;

    RKSolver::reset();
}


// ########################################################################################################################
// Explicit Runge - Kutta 4(5)
// ########################################################################################################################
void RK45::reset()
{
    // Setup RK constants before calling the base class reset
    this->C_ptr     = RK45_C_ptr;
    this->A_ptr     = RK45_A_ptr;
    this->B_ptr     = RK45_B_ptr;
    this->E_ptr     = RK45_E_ptr;
    this->E3_ptr    = nullptr;       // Not used for RK45
    this->E5_ptr    = nullptr;       // Not used for RK45
    this->P_ptr     = RK45_P_ptr;
    this->D_ptr     = nullptr;       // Not used for RK45
    this->K_ptr     = &this->K[0];
    this->order     = RK45_order;
    this->n_stages  = RK45_n_stages;
    this->len_Acols = RK45_len_Acols;
    this->len_C     = RK45_len_C;
    this->len_Pcols = RK45_len_Pcols;
    this->error_estimator_order = RK45_error_estimator_order;
    this->error_exponent = RK45_error_exponent;
    this->method_int = RK45_METHOD_INT;

    RKSolver::reset();
}


// ########################################################################################################################
// Explicit Runge-Kutta Method of order 8(5,3) due Dormand & Prince
// ########################################################################################################################
void DOP853::reset()
{
    // Setup RK constants before calling the base class reset
    this->C_ptr     = DOP853_C_ptr;
    this->A_ptr     = DOP853_A_ptr;
    this->B_ptr     = DOP853_B_ptr;
    this->E_ptr     = nullptr;        // Not used for RK23
    this->E3_ptr    = DOP853_E3_ptr;
    this->E5_ptr    = DOP853_E5_ptr;
    this->P_ptr     = nullptr;        // TODO: Not implemented
    this->D_ptr     = nullptr;        // TODO: Not implemented
    this->K_ptr     = &this->K[0];
    this->order     = DOP853_order;
    this->n_stages  = DOP853_n_stages;
    this->len_Acols = DOP853_A_cols;
    this->len_C     = DOP853_len_C;
    this->error_estimator_order = DOP853_error_estimator_order;
    this->error_exponent = DOP853_error_exponent;
    this->method_int = DOP853_METHOD_INT;

    RKSolver::reset();
}


void DOP853::p_estimate_error()
{
    double error_norm3 = 0.0;
    double error_norm5 = 0.0;

    // Initialize rtol and atol
    double rtol = this->rtols_ptr[0];
    double atol = this->atols_ptr[0];

    for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
    {
        if (this->use_array_rtols)
        {
            rtol = this->rtols_ptr[y_i];
        }

        if (this->use_array_atols)
        {
            atol = this->atols_ptr[y_i];
        }

        const unsigned int stride_K = y_i * this->n_stages_p1;

        // Dot product between K and E3 & E5 (sum over n_stages + 1; for DOP853 n_stages = 12
        // n = 0
        double temp_double = this->K_ptr[stride_K];
        double error_dot3  = this->E3_ptr[0] * temp_double;
        double error_dot5  = this->E5_ptr[0] * temp_double;

        // n = 1
        temp_double = this->K_ptr[stride_K + 1];
        error_dot3 += this->E3_ptr[1] * temp_double;
        error_dot5 += this->E5_ptr[1] * temp_double;

        // n = 2
        temp_double = this->K_ptr[stride_K + 2];
        error_dot3 += this->E3_ptr[2] * temp_double;
        error_dot5 += this->E5_ptr[2] * temp_double;

        // n = 3
        temp_double = this->K_ptr[stride_K + 3];
        error_dot3 += this->E3_ptr[3] * temp_double;
        error_dot5 += this->E5_ptr[3] * temp_double;

        // n = 4
        temp_double = this->K_ptr[stride_K + 4];
        error_dot3 += this->E3_ptr[4] * temp_double;
        error_dot5 += this->E5_ptr[4] * temp_double;

        // n = 5
        temp_double = this->K_ptr[stride_K + 5];
        error_dot3 += this->E3_ptr[5] * temp_double;
        error_dot5 += this->E5_ptr[5] * temp_double;

        // n = 6
        temp_double = this->K_ptr[stride_K + 6];
        error_dot3 += this->E3_ptr[6] * temp_double;
        error_dot5 += this->E5_ptr[6] * temp_double;

        // n = 7
        temp_double = this->K_ptr[stride_K + 7];
        error_dot3 += this->E3_ptr[7] * temp_double;
        error_dot5 += this->E5_ptr[7] * temp_double;

        // n = 8
        temp_double = this->K_ptr[stride_K + 8];
        error_dot3 += this->E3_ptr[8] * temp_double;
        error_dot5 += this->E5_ptr[8] * temp_double;

        // n = 9
        temp_double = this->K_ptr[stride_K + 9];
        error_dot3 += this->E3_ptr[9] * temp_double;
        error_dot5 += this->E5_ptr[9] * temp_double;

        // n = 10
        temp_double = this->K_ptr[stride_K + 10];
        error_dot3 += this->E3_ptr[10] * temp_double;
        error_dot5 += this->E5_ptr[10] * temp_double;

        // n = 11
        temp_double = this->K_ptr[stride_K + 11];
        error_dot3 += this->E3_ptr[11] * temp_double;
        error_dot5 += this->E5_ptr[11] * temp_double;

        // n = 12
        temp_double = this->K_ptr[stride_K + 12];
        error_dot3 += this->E3_ptr[12] * temp_double;
        error_dot5 += this->E5_ptr[12] * temp_double;

        // Find scale of y for error calculations
        const double scale_inv = 1.0 / (atol + std::fmax(std::fabs(this->y_old_ptr[y_i]), std::fabs(this->y_now_ptr[y_i])) * rtol);

        // We need the absolute value but since we are taking the square, it is guaranteed to be positive.
        // TODO: This will need to change if CySolver ever accepts complex numbers
        // error_norm_abs = fabs(error_dot_1)
        error_dot3 *= scale_inv;
        error_dot5 *= scale_inv;

        error_norm3 += (error_dot3 * error_dot3);
        error_norm5 += (error_dot5 * error_dot5);
    }

    // Check if errors are zero
    if ((error_norm5 == 0.0) && (error_norm3) == 0.0)
    {
        this->error_norm = 0.0;
    }
    else
    {
        double error_denom = error_norm5 + 0.01 * error_norm3;
        this->error_norm = this->step_size * error_norm5 / std::sqrt(error_denom * this->num_y_dbl);
    }
}


// ########################################################################################################################
// Dense Output Implementations
// ########################################################################################################################
RKDenseOutput::RKDenseOutput(double t_old, double t_now, double* y_in_ptr, unsigned int num_y, unsigned int Q_order) :
        CySolverDense(t_old, t_now, y_in_ptr, num_y),
        Q_order(Q_order)
{
    this->step = t_now - t_old;
}

void RKDenseOutput::call(double t_interp, double* y_interped)
{
    double step_factor = (t_interp - this->t_old) / this->step;

    // SciPy Step:: p = np.tile(x, self.order + 1) (scipy order is Q_order - 1)
    double cumlative_prod = step_factor;

    // y = y_old + Q dot p.
    // Q has shape of (n_stages + 1, num_y)
    for (unsigned int y_i = 0; y_i < this->num_y; y_i++)
    {
        unsigned int Q_stride = this->Q_order * y_i;

        // Initialize dot product
        double temp_double = this->Q_ptr[Q_stride + 0] * cumlative_prod;

        for (unsigned int i = 1; i < (this->Q_order); i++)
        {
            cumlative_prod *= step_factor;
            temp_double += this->Q_ptr[Q_stride + i] * cumlative_prod;
        }
        y_interped[y_i] = this->y_stored_ptr[y_i] + this->step * temp_double;
    }
}