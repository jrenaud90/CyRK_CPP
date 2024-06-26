#include "cysolution.hpp"


// Constructors
CySolverResult::CySolverResult() {}
CySolverResult::CySolverResult(size_t num_y, size_t num_dy, size_t expected_size)
{
    // Pull out known information
    this->num_y = num_y;
    this->num_dy = num_dy;
    this->num_dy_dbl = num_dy;

    // Initialize other parameters
    this->size = 0;
    this->success = false;
    this->reset_called = false;
    this->error_code = 0;
    this->update_message("CySolverResult Initialized.");

    // Initialize storage
    this->original_expected_size = expected_size;
    this->storage_capacity = expected_size;
}


// Deconstructors
CySolverResult::~CySolverResult()
{
    // Vector header inforamtion is stack allocated, no need to delete them.
}


// Protected methods
void CySolverResult::p_expand_storage()
{
    double new_storage_size_dbl = std::floor(DYNAMIC_GROWTH_RATE * this->storage_capacity);
    // Check if this new size is okay.
    if ((new_storage_size_dbl / this->num_dy_dbl) > SIZE_MAX_DBL)
    {
        this->error_code = -11;
        this->update_message("Value Error: Requested new vector size is larger than the limits set by the system (specifically the max of size_t).");
    }
    else {
        this->storage_capacity = new_storage_size_dbl;
        try
        {
            this->time_domain_ref.reserve(this->storage_capacity);
            this->solution_ref.reserve(this->storage_capacity * this->num_dy);
        }
        catch (std::bad_alloc const&)
        {
            this->error_code = -12;
            this->update_message("Memory Error: Malloc failed when reserving additional memory for storage vectors.");
        }
    }
}



// Public methods
void CySolverResult::reset()
{
    // Inititalize the storage array
    if (this->reset_called)
    {
        // The storage array may have already been set. Delete any data and reset it.
        this->time_domain_ref.clear();
        this->solution_ref.clear();
    }

    // Set the storage size to the original expected size.
    this->storage_capacity = this->original_expected_size;

    // Reserve the memory for the vectors
    try
    {
        this->time_domain_ref.reserve(this->storage_capacity);
        this->solution_ref.reserve(this->storage_capacity * this->num_dy);
    }
    catch (std::bad_alloc const&)
    {
        this->error_code = -12;
        this->update_message("Memory Error: Malloc failed when reserving initial memory for storage vectors.\n");
    }

    this->reset_called = true;
}

void CySolverResult::save_data(double new_t, double* new_solution_y, double* new_solution_dy)
{
    this->size++;

    if (this->size > this->storage_capacity)
    {
        // There is not enough room in the storage vectors. Expand them.
        this->p_expand_storage();
    }

    this->time_domain_ref.push_back(new_t);

    // We use num_dy instead of num_y in the following loop because the user may have requested to
    // capture additional output from the differential equation.
    for (size_t dy_i = 0; dy_i < this->num_dy; dy_i++)
    {
        if (dy_i < this->num_y)
        {
            // Store dependent results
            this->solution_ref.push_back(new_solution_y[dy_i]);
        }
        else {
            // Store extra output
            this->solution_ref.push_back(new_solution_dy[dy_i]);
        }
    }
}

void CySolverResult::update_message(const char* new_message)
{
    std::strcpy(this->message_ptr, new_message);
}
