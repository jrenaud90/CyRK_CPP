#include "cysolver.hpp"


void find_max_num_steps(
    size_t num_y,
    size_t num_extra,
    size_t max_num_steps,
    size_t max_ram_MB,
    bool capture_extra,
    bool* user_provided_max_num_steps,
    size_t* max_num_steps_touse) {

    // Determine the maximum number of steps permitted during integration run.
    double max_num_steps_ram_dbl;
    max_num_steps_ram_dbl = max_ram_MB * (1000 * 1000);

    // Divide by number of dependnet and extra variables that will be stored. The extra "1" is for the time domain.
    if (capture_extra)
    {
        max_num_steps_ram_dbl /= (sizeof(double) * (1.0 + num_y + num_extra));
    }
    else {
        max_num_steps_ram_dbl /= (sizeof(double) * (1.0 + num_y));
    }
    size_t max_num_steps_ram = std::floor(max_num_steps_ram_dbl);

    // Parse user-provided max number of steps
    user_provided_max_num_steps[0] = false;
    if (max_num_steps == 0)
    {
        // No user input; use ram-based value
        max_num_steps_touse[0] = max_num_steps_ram;
    }
    else {
        if (max_num_steps > max_num_steps_ram)
        {
            max_num_steps_touse[0] = max_num_steps_ram;
        }
        else {
            user_provided_max_num_steps[0] = true;
            max_num_steps_touse[0] = max_num_steps;
        }
    }

    // Make sure that max number of steps does not exceed size_t limit
    if (max_num_steps_touse[0] > (MAX_SIZET_SIZE / 10))
    {
        max_num_steps_touse[0] = (MAX_SIZET_SIZE / 10);
    }
}


void CySolverResult::expand_storage()
{
    double new_storage_size_dbl = std::floor(DYNAMIC_GROWTH_RATE * this->storage_capacity);
    // Check if this new size is okay.
    if ((new_storage_size_dbl / this->num_dy_dbl) > SIZE_MAX_DBL)
    {
        this->error_code = -4;
        std::strcpy(this->message_ptr, "Attempted to expand solution storage beyond the limits of `size_t`.");
    }
    else {
        this->storage_capacity = new_storage_size_dbl;
        this->time_domain_ptr->reserve(this->storage_capacity);
        this->solution_ptr->reserve(this->storage_capacity * this->num_dy);
    }
}
CySolverResult::CySolverResult() {}
CySolverResult::CySolverResult(size_t num_y, size_t num_dy, size_t expected_size)
{
    // Set pointers to our storage
    this->message_ptr = &this->message[0];
    this->time_domain_ptr = &this->time_domain;
    this->solution_ptr = &this->solution;

    // Pull out known information
    this->num_y = num_y;
    this->num_dy = num_dy;
    this->num_dy_dbl = num_dy;

    // Initialize other parameters
    this->size = 0;
    this->success = false;
    this->error_code = 0;
    std::strcpy(this->message_ptr, "CySolverResult Initialized.");

    // Initialize storage
    this->storage_capacity = expected_size;

    // Initialize the storage arrays to the expected size
    this->time_domain_ptr->reserve(expected_size);
    this->solution_ptr->reserve(expected_size * this->num_dy);
}
void CySolverResult::save_data(double new_t, double* new_solution_y, double* new_solution_dy, bool success, int error_code, char* message)
{
    if (success)
    {
        // Update size and check if there is enough storage in the vectors.
        this->size++;

        if (this->size > this->storage_capacity)
        {
            this->expand_storage();
        }

        if (this->error_code == 0) // Checks to make sure there was no issue with the storage expansion.
        {
            this->time_domain_ptr->push_back(new_t);

            // We use num_dy instead of num_y in the following loop because the user may have requested to
            // capture additional output from the differential equation.
            for (size_t dy_i = 0; dy_i < this->num_dy; dy_i++)
            {
                if (dy_i < this->num_y)
                {
                    // Store dependent results
                    this->solution_ptr->push_back(new_solution_y[dy_i]);
                }
                else {
                    // Store extra output
                    this->solution_ptr->push_back(new_solution_dy[dy_i]);
                }
            }
        }
        else {
            // Nothing was saved so decrement counter;
            this->size--;
        }

    }
    else
    {
        // Something failed. Don't save this data set and update the code and message
        this->success = false;
        this->error_code = error_code;
        std::sprintf(this->message_ptr, "Integration failed at time step %d. Message:\n %s", this->size + 1, message);
    }
}

void CySolverBase::step_implementation() { // Overwritten by subclasses.};
}
CySolverBase::CySolverBase() {}
CySolverBase::CySolverBase(
    // Input variables
    DiffeqFuncType diffeq_ptr,
    double t_start,
    double t_end,
    double* y0_ptr,
    size_t num_y,
    bool capture_extra,
    size_t num_extra,
    double* args_ptr,
    size_t max_num_steps,
    size_t max_ram_MB)
{
    // Assume no errors for now.
    this->error_code = 0;
    this->reset_called = false;
    this->status = 0;
    this->message_ptr = &message[0];
    std::strcpy(this->message_ptr, "CySolverBase Initializing.");

    // Check for errors
    if (capture_extra)
    {
        if (num_extra == 0)
        {
            this->error_code = -1;
            std::strcpy(this->message_ptr, "CySolverBase Attribute Error: `capture_extra` set to True, but `num_extra` set to 0.");
        }
        else if (num_extra > 50)
        {
            this->error_code = -1;
            std::strcpy(this->message_ptr, "CySolverBase Attribute Error: `num_extra` exceeds the maximum supported value of 50.");
        }

    }
    else if (num_extra > 0)
    {
        this->error_code = -1;
        std::strcpy(this->message_ptr, "CySolverBase Attribute Error: `capture_extra` set to False, but `num_extra` > 0.");
    }

    if (num_y > 50)
    {
        this->error_code = -1;
        std::strcpy(this->message_ptr, "CySolverBase Attribute Error: `num_y` exceeds the maximum supported value of 50.");
    }
    else if (num_y == 0)
    {
        this->error_code = -1;
        std::strcpy(this->message_ptr, "CySolverBase Attribute Error: Integration completed. `num_y` = 0 so nothing to integrate.");
    }

    // Parse differential equation
    this->diffeq_ptr = diffeq_ptr;
    this->args_ptr = args_ptr;

    // Parse capture extra information
    this->capture_extra = capture_extra;
    this->num_extra = num_extra;

    // Parse y values
    this->num_y = num_y;
    this->num_y_dbl = this->num_y;
    this->num_y_sqrt = std::sqrt(this->num_y_dbl);
    this->num_dy = this->num_y + this->num_extra;
    // Make a copy of y0
    for (size_t y_i = 0; y_i < this->num_y; y_i++)
    {
        this->y0_ptr[y_i] = y0_ptr[y_i];
    }

    // Parse time information
    this->t_start = t_start;
    this->t_end = t_end;
    this->t_delta = t_end - t_start;
    this->t_delta_abs = std::fabs(this->t_delta);
    if (this->t_delta >= 0.0)
    {
        // Forward integration
        this->direction_flag = true;
        this->direction_inf = INF;
    }
    else {
        // Backward integration
        this->direction_flag = false;
        this->direction_inf = -INF;
    }

    // Parse maximum number of steps
    find_max_num_steps(
        this->num_y,
        num_extra,
        max_num_steps,
        max_ram_MB,
        capture_extra,
        &this->user_provided_max_num_steps,
        &this->max_num_steps);
}

void CySolverBase::diffeq()
{
    // Call differential equation
    this->diffeq_ptr(this->dy_now_ptr, this->t_now, this->y_now_ptr, this->args_ptr);
}

void CySolverBase::reset()
{
    double temp_double;

    this->reset_called = false;

    // Reset time
    this->t_now = this->t_start;
    this->t_old = this->t_start;
    this->len_t = 1;

    // Reset ys
    for (size_t y_i = 0; y_i < this->num_y; y_i++)
    {
        temp_double = this->y0_ptr[y_i];
        this->y_now_ptr[y_i] = temp_double;
        this->y_old_ptr[y_i] = temp_double;
    }

    // Call differential equation to set dy0
    this->diffeq();

    // Update dys
    for (size_t dy_i = 0; dy_i < this->num_dy; dy_i++)
    {
        this->dy_old_ptr[dy_i] = this->dy_now_ptr[dy_i];
    }

    this->reset_called = true;

}

void CySolverBase::take_step()
{

    if (!this->reset_called)
    {
        // Reset must be called first.
        this->reset();
    }

    if (this->status == 0)
    {
        if (this->t_now == this->t_end)
        {
            this->t_old = this->t_end;
            this->status = 1;
        }
        else if (this->len_t > this->max_num_steps)
        {
            if (user_provided_max_num_steps)
            {
                // Maximum number of steps reached (as set by user).
                this->status = -2;
                this->error_code = -5;
                std::strcpy(this->message_ptr, "Maximum number of steps (set by user) exceeded during integration.\n");
            }
            else {
                // Maximum number of steps reached (as set by RAM limitations).
                this->status = -3;
                this->error_code = -5;
                std::strcpy(this->message_ptr, "Maximum number of steps (set by system architecture) exceeded during integration.\n");
            }
        }
        else
        {
            // ** Make call to solver's step implementation **
            this->step_implementation();
            this->len_t++;
        }
    }
    else
    {
        this->error_code = -5;
        std::strcpy(this->message_ptr, "Warning: Step called when status != 0.");
    }
}

void RKSolver::estimate_error()
{

    size_t stride;
    double error_dot, scale;

    // Initialize rtol and atol
    double rtol = this->rtols_ptr[0];
    double atol = this->atols_ptr[0];

    this->error_norm = 0.0;

    for (size_t y_i = 0; y_i < this->num_y; y_i++)
    {
        if (this->use_array_rtols)
        {
            rtol = this->rtols_ptr[y_i];
        }

        if (this->use_array_atols)
        {
            atol = this->atols_ptr[y_i];
        }

        // Find scale of y for error calculations
        scale = atol + std::fmax(std::fabs(this->y_old_ptr[y_i]), std::fabs(this->y_now_ptr[y_i])) * rtol;

        // Set last array of K equal to dydt
         this->K_ptr[this->nstages_numy + y_i] = this->dy_now_ptr[y_i];

        // Initialize
        error_dot = 0.0;

        for (size_t j = 0; j < (this->n_stages + 1); j++) {
            error_dot += this->K_ptr[j * this->num_y + y_i] * this->E_ptr[j];
        }

        // We need the absolute value but since we are taking the square, it is guaranteed to be positive.
        // TODO: This will need to change if CySolver ever accepts complex numbers
        // error_norm_abs = fabs(error_dot_1)
        error_dot *= (this->step / scale);

        this->error_norm += (error_dot * error_dot);
    }
    this->error_norm = std::sqrt(this->error_norm) / this->num_y_sqrt;
}

void RKSolver::step_implementation()
{
    // Initialize step variables
    size_t stride_1, stride_2, stride_3;
    double step_factor, time_tmp, t_delta_check;
    double temp_double, error_pow;

    // Initialize tolerances to the 0 place. If `use_array_rtols` (or atols) is set then this will change in the loop.
    double rtol = this->rtols_ptr[0];
    double atol = this->atols_ptr[0];

    // Run RK integration step
    // Determine step size based on previous loop
    // Find minimum step size based on the value of t (less floating point numbers between numbers when t is large)
    const double min_step_size = 10. * std::fabs(std::nextafter(this->t_old, this->direction_inf) - this->t_old);
    // Look for over/undershoots in previous step size
    if (this->step_size > this->max_step_size) {
        this->step_size = this->max_step_size;
    }
    else if (this->step_size < min_step_size) {
        this->step_size = min_step_size;
    }

    // Optimization variables
    // Define a very specific A (Row 1; Col 0) now since it is called consistently and does not change.
    const double A_at_10 = this->A_ptr[1 * this->len_Acols + 0];

    // Determine new step size
    bool step_accepted = false;
    bool step_rejected = false;
    bool step_error = false;

    // !! Step Loop
    while (!step_accepted) {

        // Check if step size is too small
        // This will cause integration to fail: step size smaller than spacing between numbers
        if (this->step_size < min_step_size) {
            step_error = true;
            this->status = -1;
            break;
        }

        // Move time forward for this particular step size
        if (this->direction_flag) {
            this->step = this->step_size;
            this->t_now = this->t_old + this->step;
            t_delta_check = this->t_now - this->t_end;
        }
        else {
            this->step = -this->step_size;
            this->t_now = this->t_old + this->step;
            t_delta_check = this->t_end - this->t_now;
        }

        // Check that we are not at the end of integration with that move
        if (t_delta_check > 0.0) {
            this->t_now = this->t_end;

            // If we are, correct the step so that it just hits the end of integration.
            this->step = this->t_now - this->t_old;
            if (this->direction_flag) {
                this->step_size = this->step;
            }
            else {
                this->step_size = -this->step;
            }
        }

        // !! Calculate derivative using RK method

        // t_now must be updated for each loop of s in order to make the diffeq method calls.
        // But we need to return to its original value later on. Store in temp variable.
        time_tmp = this->t_now;

        for (size_t s = 1; s < this->len_C; s++) {
            // Find the current time based on the old time and the step size.
            this->t_now = this->t_old + this->C_ptr[s] * this->step;
            stride_1 = s * this->num_y;
            stride_3 = s * this->len_Acols;

            // Dot Product (K, a) * step
            if (s == 1) {
                for (size_t y_i = 0; y_i < this->num_y; y_i++) {
                    // Set the first column of K
                    temp_double = this->dy_old_ptr[y_i];
                    // K[0, :] == first part of the array
                    this->K_ptr[y_i] = temp_double;

                    // Calculate y_new for s==1
                    this->y_now_ptr[y_i] = this->y_old_ptr[y_i] + (temp_double * A_at_10 * this->step);
                }
            }
            else {
                for (size_t j = 0; j < s; j++) {
                    temp_double = this->A_ptr[stride_3 + j] * this->step;
                    stride_2 = j * this->num_y;
                    for (size_t y_i = 0; y_i < this->num_y; y_i++) {
                        if (j == 0) {
                            // Initialize
                            this->y_now_ptr[y_i] = this->y_old_ptr[y_i];
                        }
                        this->y_now_ptr[y_i] += this->K_ptr[stride_2 + y_i] * temp_double;
                    }
                }
            }
            // Call diffeq method to update K with the new dydt
            // This will use the now updated dy_now_ptr based on the values of y_now_ptr and t_now_ptr.
            this->diffeq();

            // Update K based on the new dy values.
            for (size_t y_i = 0; y_i < this->num_y; y_i++) {
                this->K_ptr[stride_1 + y_i] = this->dy_now_ptr[y_i];
            }
        }

        // Restore t_now to its previous value.
        this->t_now = time_tmp;

        // Dot Product (K, B) * step
        for (size_t j = 0; j < this->n_stages; j++) {
            temp_double = this->B_ptr[j] * this->step;
            stride_1 = j * this->num_y;
            // We do not use rk_n_stages_plus1 here because we are chopping off the last row of K to match
            //  the shape of B.
            for (size_t y_i = 0; y_i < this->num_y; y_i++) {
                if (j == 0) {
                    // Initialize
                    this->y_now_ptr[y_i] = this->y_old_ptr[y_i];
                }
                this->y_now_ptr[y_i] += this->K_ptr[stride_1 + y_i] * temp_double;
            }
        }

        // Find final dydt for this timestep
        // This will use the now updated dy_now_ptr based on the values of y_now_ptr and t_now_ptr.
        this->diffeq();

        // Check how well this step performed by calculating its error.
        this->estimate_error();

        // Check the size of the error
        if (this->error_norm < 1.0) {
            // We found our step size because the error is low!
            // Update this step for the next time loop
            if (this->error_norm == 0.0) {
                step_factor = this->max_step_factor;
            }
            else {
                error_pow = std::pow(this->error_norm, -this->error_exponent);
                step_factor = std::fmin(this->max_step_factor, this->error_safety * error_pow);
            }

            if (step_rejected) {
                // There were problems with this step size on the previous step loop. Make sure factor does
                //   not exasperate them.
                step_factor = std::fmin(step_factor, 1.);
            }

            // Update step size
            this->step_size *= step_factor;
            step_accepted = true;
        }
        else {
            // Error is still large. Keep searching for a better step size.
            error_pow = std::pow(this->error_norm, -this->error_exponent);
            this->step_size *= std::fmax(this->min_step_factor, this->error_safety * error_pow);
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
    // Update "old" pointers
    this->t_old = t_now;
    for (size_t y_i = 0; y_i < this->num_y; y_i++) {
        this->y_old_ptr[y_i] = this->y_now_ptr[y_i];
        this->dy_old_ptr[y_i] = this->dy_now_ptr[y_i];
    }
}

RKSolver::RKSolver() {}
RKSolver::RKSolver(
    // Input variables
    DiffeqFuncType diffeq_ptr,
    double t_start,
    double t_end,
    double* y0_ptr,
    size_t num_y,
    bool capture_extra,
    size_t num_extra,
    double* args_ptr,
    size_t max_num_steps,
    size_t max_ram_MB,
    double rtol,
    double atol,
    double* rtols_ptr,
    double* atols_ptr,
    double max_step_size,
    double first_step_size) : CySolverBase(diffeq_ptr, t_start, t_end, y0_ptr, num_y, capture_extra, num_extra, args_ptr, max_num_steps, max_ram_MB)
{
    // Check for errors
    if (first_step_size != 0.0)
    {
        if (first_step_size < 0.)
        {
            this->error_code = -1;
            std::strcpy(this->message_ptr, "User-provided initial step size must be a positive number.");
        }
        else if (first_step_size > (this->t_delta_abs * 0.5))
        {
            this->error_code = -1;
            std::strcpy(this->message_ptr, "User-provided initial step size must be smaller than 50% of the time span size.");
        }
    }

    // Setup tolerances
    // User can provide an array of relative tolerances, one for each y value.
    // The length of the pointer array must be the same as y0 (and <= 50).
    this->use_array_rtols = false;
    this->use_array_atols = false;

    double temp_double;
    double min_rtol = INF;
    if (rtols_ptr)
    {
        // rtol for each y
        use_array_rtols = true;
        for (size_t y_i = 0; y_i < this->num_y; y_i++)
        {
            temp_double = rtols_ptr[y_i];
            if (temp_double < EPS_100)
            {
                temp_double = EPS_100;
            }
            min_rtol = std::fmin(min_rtol, temp_double);
            this->rtols_ptr[y_i] = temp_double;
        }
    }
    else {
        // only one rtol
        temp_double = rtol;
        if (temp_double < EPS_100)
        {
            temp_double = EPS_100;
        }
        min_rtol = std::fmin(min_rtol, temp_double);
        this->rtols_ptr[0] = temp_double;
    }

    if (atols_ptr)
    {
        // atol for each y
        use_array_atols = true;
        for (size_t y_i = 0; y_i < this->num_y; y_i++)
        {
            this->atols_ptr[y_i] = atols_ptr[y_i];
        }
    }
    else {
        // only one atol
        this->atols_ptr[0] = atol;
    }

    // Setup step size
    this->max_step_size = max_step_size;
    this->user_provided_first_step_size = first_step_size;
}

void RKSolver::reset()
{
    // Call base class reset.
    CySolverBase::reset();
    // Update stride information
    this->nstages_numy = this->n_stages * this->num_y;

    // Update initial step size
    if (this->user_provided_first_step_size == 0)
    {
        // User did not provide a step size. Try to find a good guess.
        this->calc_first_step_size();
    }
    else {
        this->step_size = this->user_provided_first_step_size;
        this->step_size_old = this->step_size;
    }

    // It is important to initialize the K variable with zeros
    for (size_t i = 0; i < (this->n_stages + 1); i++)
    {
        for (size_t y_i = 0; y_i < this->num_y; y_i++)
        {
            this->K_ptr[i * this->num_y + y_i] = 0.0;
        }
    }
}

void RKSolver::calc_first_step_size()
{
    /*
        Select an initial step size based on the differential equation.
        .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
            Equations I: Nonstiff Problems", Sec. II.4.
    */

    double d0, d1, d2, d0_abs, d1_abs, d2_abs, scale;
    double h0, h0_direction, h1;
    double y_old_tmp;

    // Initialize tolerances to the 0 place. If `use_array_rtols` (or atols) is set then this will change in the loop.
    double rtol = this->rtols_ptr[0];
    double atol = this->atols_ptr[0];

    if (this->num_y == 0)
    {
        this->step_size = INF;
    }
    else {
        // Find the norm for d0 and d1
        d0 = 0.0;
        d1 = 0.0;
        for (size_t y_i = 0; y_i < this->num_y; y_i++)
        {
            y_old_tmp = this->y_old_ptr[y_i];

            if (this->use_array_rtols)
            {
                rtol = this->rtols_ptr[y_i];
            }
            if (this->use_array_atols)
            {
                atol = this->atols_ptr[y_i];
            }

            scale = atol + std::fabs(y_old_tmp) * rtol;
            d0_abs = std::fabs(y_old_tmp / scale);
            d1_abs = std::fabs(dy_old_ptr[y_i] / scale);
            d0 += (d0_abs * d0_abs);
            d1 += (d1_abs * d1_abs);
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

        if (this->direction_flag)
        {
            h0_direction = h0;
        }
        else {
            h0_direction = -h0;
        }

        this->t_now = this->t_old + h0_direction;
        for (size_t y_i = 0; y_i < this->num_y; y_i++)
        {
            this->y_now_ptr[y_i] = this->y_old_ptr[y_i] + h0_direction * this->dy_old_ptr[y_i];
        }

        // Update dy
        this->diffeq();

        // Find the norm for d2
        d2 = 0.0;
        for (size_t y_i = 0; y_i < this->num_y; y_i++)
        {
            if (this->use_array_rtols)
            {
                rtol = this->rtols_ptr[y_i];
            }
            if (this->use_array_atols)
            {
                atol = this->atols_ptr[y_i];
            }

            scale = atol + std::fabs(this->y_old_ptr[y_i]) * rtol;
            d2_abs = std::fabs((this->dy_now_ptr[y_i] - this->dy_old_ptr[y_i]) / scale);
            d2 += (d2_abs * d2_abs);
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
        this->step_size_old = this->step_size;
    }
}


void RK23::reset()
{
    // Setup RK constants before calling the base class reset
    this->C_ptr = &RK23_C[0];
    this->A_ptr = &RK23_A[0];
    this->B_ptr = &RK23_B[0];
    this->E_ptr = &RK23_E[0];
    this->E3_ptr = nullptr;  // Not used for RK23
    this->E5_ptr = nullptr;  // Not used for RK23
    this->P_ptr = nullptr;  // TODO: Not implemented
    this->D_ptr = nullptr;  // TODO: Not implemented
    this->K_ptr = &this->K[0];
    this->order = RK23_order;
    this->error_estimator_order = RK23_error_estimator_order;
    this->error_exponent = RK23_error_exponent;
    this->n_stages = RK23_n_stages;
    this->len_Acols = RK23_len_Acols;
    this->len_C = RK23_len_C;

    RKSolver::reset();
}


void RK45::reset()
{
    // Setup RK constants before calling the base class reset
    this->C_ptr = &RK45_C[0];
    this->A_ptr = &RK45_A[0];
    this->B_ptr = &RK45_B[0];
    this->E_ptr = &RK45_E[0];
    this->E3_ptr = nullptr;  // Not used for RK23
    this->E5_ptr = nullptr;  // Not used for RK23
    this->P_ptr = nullptr;  // TODO: Not implemented
    this->D_ptr = nullptr;  // TODO: Not implemented
    this->K_ptr = &this->K[0];
    this->order = RK45_order;
    this->error_estimator_order = RK45_error_estimator_order;
    this->error_exponent = RK45_error_exponent;
    this->n_stages = RK45_n_stages;
    this->len_Acols = RK45_len_Acols;
    this->len_C = RK45_len_C;

    RKSolver::reset();
}
