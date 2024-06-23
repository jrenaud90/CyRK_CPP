#include <vector>
#include <limits>
#include <cmath>
#include <cstring>
#include <cstdio>

// Integration Constants
// Multiply steps computed from asymptotic behaviour of errors by this.
const double SAFETY = 0.9;  // Error coefficient factor (1 == use calculated error; < 1 means be conservative).
const double MIN_FACTOR = 0.2;  // Minimum allowed decrease in a step size.
const double MAX_FACTOR = 10.;  // Maximum allowed increase in a step size.
const double INF = std::numeric_limits<double>::infinity();
const double MAX_STEP = INF;
const double EPS = std::numeric_limits<double>::epsilon();
const double EPS_10 = EPS * 10.0;
const double EPS_100 = EPS * 100.0;
const size_t MAX_SIZET_SIZE = std::numeric_limits<size_t>::max();
const size_t MAX_INT_SIZE = std::numeric_limits<int>::max();

// Solution parameters
const double DYNAMIC_GROWTH_RATE = 1.618;
// To learn why the golden ratio is used, read this:
// https://stackoverflow.com/questions/1100311/what-is-the-ideal-growth-rate-for-a-dynamically-allocated-array
const double SIZE_MAX_DBL = 0.99 * SIZE_MAX;


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


typedef void (*DiffeqFuncType)(double*, double, double*, double*);

class CySolverResult {

protected:
    // ** Attributes **
    // Message storage
    char message[512];

    // Metadata
    size_t num_y;
    size_t num_dy;
    double num_dy_dbl;

    // Current storage information
    size_t storage_capacity;

    // ** Methods **
    void expand_storage() {

        double new_storage_size_dbl = std::floor(DYNAMIC_GROWTH_RATE * this->storage_capacity);
        // Check if this new size is okay.
        if ((new_storage_size_dbl / this->num_dy_dbl) > SIZE_MAX_DBL)
        {
            this->error_code = -4;
            strcpy(this->message_ptr, "Attempted to expand solution storage beyond the limits of `size_t`.");
        }
        else {
            this->storage_capacity = new_storage_size_dbl;
            this->time_domain_ptr->reserve(this->storage_capacity);
            this->solution_ptr->reserve(this->storage_capacity * this->num_dy);
        }
    }

public:
    // ** Attributes **
    // Status information
    char* message_ptr;
    bool success;
    int error_code;
    size_t size;

    // Storage for arrays
    std::vector<double> time_domain;
    std::vector<double> solution;
    std::vector<double>* time_domain_ptr;
    std::vector<double>* solution_ptr;

    // ** Methods **
    CySolverResult() {}
    CySolverResult(size_t num_y, size_t num_dy, size_t expected_size) {
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

    void save_data(double new_t, double* new_solution_y, double* new_solution_dy, bool success, int error_code, char* message) {
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
        else {
            // Something failed. Don't save this data set and update the code and message
            this->success = false;
            this->error_code = error_code;
            std::sprintf(this->message_ptr, "Integration failed at time step %d. Message:\n %s", this->size + 1, message);
        }

    }
};

class CySolverBase {

protected:
    // ** Attributes **
    // Status variables
    bool reset_called;
    char message[512];

    // Time variables
    double t_start;
    double t_end;
    double t_old;
    double t_delta;
    double t_delta_abs;
    double direction_inf;
    bool direction_flag;

    // Dependent variables
    size_t num_y;
    double num_y_dbl;
    double num_y_sqrt;
    size_t num_dy;
    // The size of the stack allocated tracking arrays is equal to the maximum allowed `num_y` (50).
    double y0[50];
    double y_old[50];
    double y_now[50];
    // For dy, both the dy/dt and any extra outputs are stored. So the maximum size is `num_y` (50) + `num_extra` (50)
    double dy_old[100];
    double dy_now[100];

    // dy_now_ptr and y_now_ptr are declared in public.
    double* y0_ptr = &y0[0];
    double* y_old_ptr = &y_old[0];
    double* dy_old_ptr = &dy_old[0];

    // Integration step information
    size_t len_t;
    size_t max_num_steps;
    bool user_provided_max_num_steps;

    // Information on capturing extra information during integration.
    bool capture_extra;
    size_t num_extra;

    // Differential equation information
    double* args_ptr;
    DiffeqFuncType diffeq_ptr;

    // ** Methods **
    virtual void step_implementation() {
        // Overwritten by subclasses.
    }

public:

    // ** Attributes **
    int status;
    // Error Codes:
    // -1 : Attribute Error
    // -2 : Value Error
    // -3 : Not Implemented Error
    // -4 : <Reserved>
    // -5 : Runtime error
    int error_code;
    char* message_ptr;

    // State properties
    double t_now;
    double* y_now_ptr = &y_now[0];
    double* dy_now_ptr = &dy_now[0];

    // ** Methods **
    // Constructors
    CySolverBase() {}
    CySolverBase(
        // Input variables
        DiffeqFuncType diffeq_ptr,
        double t_start,
        double t_end,
        double* y0_ptr,
        size_t num_y,
        bool capture_extra = false,
        size_t num_extra = 0,
        double* args_ptr = nullptr,
        size_t max_num_steps = 0,
        size_t max_ram_MB = 2000
    ) {
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

    void diffeq() {
        // Call differential equation
        this->diffeq_ptr(this->dy_now_ptr, this->t_now, this->y_now_ptr, this->args_ptr);
    }

    virtual void reset() {

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

    void take_step() {

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
                    strcpy(this->message_ptr, "Maximum number of steps (set by user) exceeded during integration.\n");
                }
                else {
                    // Maximum number of steps reached (as set by RAM limitations).
                    this->status = -3;
                    this->error_code = -5;
                    strcpy(this->message_ptr, "Maximum number of steps (set by system architecture) exceeded during integration.\n");
                }
            }
            else {
                // ** Make call to solver's step implementation **
                this->step_implementation();
                this->len_t++;
            }
        }
        else {
            this->error_code = -5;
            strcpy(this->message_ptr, "Warning: Step called when status != 0.");
        }
    }
};


class RKSolver : public CySolverBase {

protected:
    // ** Attributes **

    // Step globals
    double error_safety = MIN_FACTOR;
    double min_step_factor = MAX_FACTOR;
    double max_step_factor = SAFETY;

    // RK constants
    int order = 0;
    int error_estimator_order = 0;
    double error_exponent = std::nan("");
    size_t n_stages = 0;
    size_t len_Acols = 0;
    size_t len_C = 0;

    // The maximum size for each container is chosen.
    // When a container cares about the number of y values, as is the case for K, then the max supported num_y of 50
    // is used.
    double C[1] = { std::nan("") };
    double A[1] = { std::nan("") };
    double B[1] = { std::nan("") };
    double E3[1] = { std::nan("") };
    double E5[1] = { std::nan("") };
    double P[1] = { std::nan("") };
    double D[1] = { std::nan("") };

    // K is not const
    double K[1] = { 0.0 };

    double* C_ptr;
    double* A_ptr;
    double* B_ptr;
    double* E_ptr;
    double* E3_ptr;
    double* E5_ptr;
    double* P_ptr;
    double* D_ptr;

    // K is not const
    double* K_ptr;

    // Tolerances
    // For the same reason num_y is limited, the total number of tolerances are limited.
    double rtols[50];
    double atols[50];
    double* rtols_ptr = &rtols[0];
    double* atols_ptr = &atols[0];
    bool use_array_rtols;
    bool use_array_atols;

    // Step size parameters
    double step;
    double step_size_old;
    double step_size;
    double max_step_size;
    double user_provided_first_step_size;

    // Error estimate
    double error_norm;

    // ** Methods **
    void estimate_error() {

        double error_dot, scale;

        // Initialize rtol and atol
        double rtol = this->rtols_ptr[0];
        double atol = this->atols_ptr[0];

        this->error_norm = 0.0;
        for (size_t y_i = 0; y_i < this->num_y; y_i++) {
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
            this->K_ptr[this->n_stages * this->num_y + y_i] = this->dy_now_ptr[y_i];

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


    void step_implementation() override {
        // Initialize step variables
        double step_factor, time_tmp, t_delta_check;
        double temp_double, error_pow;

        // Initialize tolerances to the 0 place. If `use_array_rtols` (or atols) is set then this will change in the loop.
        double rtol = this->rtols_ptr[0];
        double atol = this->atols_ptr[0];

        // Run RK integration step
        // Determine step size based on previous loop
        // Find minimum step size based on the value of t (less floating point numbers between numbers when t is large)
        double min_step_size = 10. * std::fabs(std::nextafter(this->t_old, this->direction_inf) - this->t_old);
        // Look for over/undershoots in previous step size
        if (this->step_size > this->max_step_size) {
            this->step_size = this->max_step_size;
        }
        else if (this->step_size < min_step_size) {
            this->step_size = min_step_size;
        }

        // Optimization variables
        // Define a very specific A (Row 1; Col 0) now since it is called consistently and does not change.
        double A_at_10 = this->A_ptr[1 * this->len_Acols + 0];

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
                        temp_double = this->A_ptr[s * this->len_Acols + j] * this->step;
                        for (size_t y_i = 0; y_i < this->num_y; y_i++) {
                            if (j == 0) {
                                // Initialize
                                this->y_now_ptr[y_i] = this->y_old_ptr[y_i];
                            }
                            this->y_now_ptr[y_i] += this->K_ptr[j * this->num_y + y_i] * temp_double;
                        }
                    }
                }
                // Call diffeq method to update K with the new dydt
                // This will use the now updated dy_now_ptr based on the values of y_now_ptr and t_now_ptr.
                this->diffeq();

                // Update K based on the new dy values.
                for (size_t y_i = 0; y_i < this->num_y; y_i++) {
                    this->K_ptr[s * this->num_y + y_i] = this->dy_now_ptr[y_i];
                }
            }

            // Restore t_now to its previous value.
            this->t_now = time_tmp;

            // Dot Product (K, B) * step
            for (size_t j = 0; j < this->n_stages; j++) {
                temp_double = this->B_ptr[j] * this->step;
                // We do not use rk_n_stages_plus1 here because we are chopping off the last row of K to match
                //  the shape of B.
                for (size_t y_i = 0; y_i < this->num_y; y_i++) {
                    if (j == 0) {
                        // Initialize
                        this->y_now_ptr[y_i] = this->y_old_ptr[y_i];
                    }
                    this->y_now_ptr[y_i] += this->K_ptr[j * this->num_y + y_i] * temp_double;
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

public:
    // ** Attributes **

    // ** Methods **
    // Constructors
    RKSolver() {}
    RKSolver(
        // Input variables
        DiffeqFuncType diffeq_ptr,
        double t_start,
        double t_end,
        double* y0_ptr,
        size_t num_y,
        bool capture_extra = false,
        size_t num_extra = 0,
        double* args_ptr = nullptr,
        size_t max_num_steps = 0,
        size_t max_ram_MB = 2000,
        double rtol = 1.0e-3,
        double atol = 1.0e-6,
        double* rtols_ptr = nullptr,
        double* atols_ptr = nullptr,
        double max_step_size = MAX_STEP,
        double first_step_size = 0.0
    ) : CySolverBase(diffeq_ptr, t_start, t_end, y0_ptr, num_y, capture_extra, num_extra, args_ptr, max_num_steps, max_ram_MB) {

        // Check for errors
        if (first_step_size != 0.0)
        {
            if (first_step_size < 0.)
            {
                this->error_code = -1;
                strcpy(this->message_ptr, "User-provided initial step size must be a positive number.");
            }
            else if (first_step_size > (this->t_delta_abs * 0.5))
            {
                this->error_code = -1;
                strcpy(this->message_ptr, "User-provided initial step size must be smaller than 50% of the time span size.");
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

    void reset() override {
        // Call base class reset.
        CySolverBase::reset();

        // Setup RK pointers
        this->C_ptr = &this->C[0];
        this->A_ptr = &this->A[0];
        this->B_ptr = &this->B[0];
        // E and E3 point to the same memory. E3 is only used for DOP853 method other RK methods use E.
        this->E3_ptr = &this->E3[0];
        this->E_ptr = &this->E3[0];
        this->E5_ptr = &this->E5[0];
        this->P_ptr = &this->P[0];
        this->D_ptr = &this->D[0];
        this->K_ptr = &this->K[0];

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

    void calc_first_step_size() {
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
};


class RK23 : public RKSolver {

protected:
    int order = 3;
    size_t n_stages = 3;
    size_t len_Acols = 3;
    size_t len_C = 3;
    int error_estimator_order = 2;
    double error_exponent = 1.0 / (2.0 + 1.0);  // Defined as 1 / (error_order + 1)

    static double A[9];
    double B[3] = {
        2.0 / 9.0,
        1.0 / 3.0,
        4.0 / 9.0
    };
    double C[3] = {
        0.0,
        1.0 / 2.0,
        3.0 / 4.0
    };
    double E3[4] = {
        5.0 / 72.0,
        -1.0 / 12.0,
        -1.0 / 9.0,
        1.0 / 8.0
    };
    double K[4 * 50];

public:
    // Copy over base class constructors
    using RKSolver::RKSolver;

};

double RK23::A = {
    // A - Row 0
    0.0,
    0.0,
    0.0,
    // A - Row 1
    1.0 / 2.0,
    0.0,
    0.0,
    // A - Row 2
    0.0,
    3.0 / 4.0,
    0.0
};

class RK45 : public RKSolver {

protected:
    int order = 5;
    size_t n_stages = 6;
    size_t len_Acols = 5;
    size_t len_C = 6;
    int error_estimator_order = 4;
    double error_exponent = 1.0 / (4.0 + 1.0);  // Defined as 1 / (error_order + 1)

    double A[30] = {
        // A - Row 0
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        // A - Row 1
        1.0 / 5.0,
        0.0,
        0.0,
        0.0,
        0.0,
        // A - Row 2
        3.0 / 40.0,
        9.0 / 40.0,
        0.0,
        0.0,
        0.0,
        // A - Row 3
        44.0 / 45.0,
        -56.0 / 15.0,
        32.0 / 9.0,
        0.0,
        0.0,
        // A - Row 4
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        // A - Row 5
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0
    };
    double B[6] = {
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0
    };
    double C[6] = {
        0.0,
        1.0 / 5.0,
        3.0 / 10.0,
        4.0 / 5.0,
        8.0 / 9.0,
        1.0
    };
    double E3[7] = {
        -71.0 / 57600.0,
        0.0,
        71.0 / 16695.0,
        -71.0 / 1920.0,
        17253.0 / 339200.0,
        -22.0 / 525.0,
        1.0 / 40.0
    };
    double K[7 * 50];

public:
    // Copy over base class constructors
    using RKSolver::RKSolver;
};