/* Methods to store and retrieve data saved by CySolver */

#include "cysolution.hpp"


inline void round_to_2(size_t &initital_value)
{
    /* Rounds the initial value to the nearest power of 2 */
    // Method is the fastest for 64-bit numbers
    initital_value--;
    initital_value |= initital_value >> 1;
    initital_value |= initital_value >> 2;
    initital_value |= initital_value >> 4;
    initital_value |= initital_value >> 8;
    initital_value |= initital_value >> 16;
    initital_value |= initital_value >> 32;
    initital_value++;
}

// Constructors
CySolverResult::CySolverResult()
{

}

CySolverResult::CySolverResult(const int num_y, const int num_extra, const size_t expected_size) :
        num_y(num_y),
        num_extra(num_extra),
        error_code(0)        
{
    // Round expected size and store it.
    this->original_expected_size = expected_size;
    round_to_2(this->original_expected_size);

    // num_dy will be larger than num_y if the user wishes to capture extra output during integration.
    this->capture_extra = this->num_extra > 0;
    this->num_dy        = this->num_y + this->num_extra;
    this->num_dy_dbl    = (double)this->num_dy;

    // Initialize other parameters
    this->update_message("CySolverResult Initialized.");
}


// Deconstructors
CySolverResult::~CySolverResult()
{
    // Vector header inforamtion is stack allocated, no need to delete them.
    // The data itself is heap allocated but the vector class will handle that.
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
    else
    {
        this->storage_capacity = (size_t)new_storage_size_dbl;
        round_to_2(this->storage_capacity);
        try
        {
            this->time_domain.reserve(this->storage_capacity);
            this->solution.reserve(this->storage_capacity * this->num_dy);
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
        this->time_domain.clear();
        this->solution.clear();
    }

    // Set the storage size to the original expected size.
    this->storage_capacity = this->original_expected_size;

    // Reserve the memory for the vectors
    try
    {
        this->time_domain.reserve(this->storage_capacity);
        this->solution.reserve(this->storage_capacity * this->num_dy);
    }
    catch (std::bad_alloc const&)
    {
        this->error_code = -12;
        this->update_message("Memory Error: Malloc failed when reserving initial memory for storage vectors.\n");
    }

    this->reset_called = true;
}

void CySolverResult::save_data(const double new_t, double* const new_solution_y_ptr, double* const new_solution_dy_ptr)
{
    this->size++;

    if (this->size > this->storage_capacity)
    {
        // There is not enough room in the storage vectors. Expand them.
        this->p_expand_storage();
    }

    this->time_domain.push_back(new_t);

    // Store y values (use push back if num_y less than 5 otherwise use insert)
    switch (this->num_y)
    {
        case 0:
            break;
        case 1:
            this->solution.push_back(new_solution_y_ptr[0]);
            break;
        case 2:
            this->solution.push_back(new_solution_y_ptr[0]);
            this->solution.push_back(new_solution_y_ptr[1]);
            break;
        case 3:
            this->solution.push_back(new_solution_y_ptr[0]);
            this->solution.push_back(new_solution_y_ptr[1]);
            this->solution.push_back(new_solution_y_ptr[2]);
            break;
        case 4:
            this->solution.push_back(new_solution_y_ptr[0]);
            this->solution.push_back(new_solution_y_ptr[1]);
            this->solution.push_back(new_solution_y_ptr[2]);
            this->solution.push_back(new_solution_y_ptr[3]);
            this->solution.push_back(new_solution_y_ptr[4]);
            break;
        case 5:
            this->solution.push_back(new_solution_y_ptr[0]);
            this->solution.push_back(new_solution_y_ptr[1]);
            this->solution.push_back(new_solution_y_ptr[2]);
            this->solution.push_back(new_solution_y_ptr[3]);
            this->solution.push_back(new_solution_y_ptr[4]);
            break;
        case 6:
            this->solution.push_back(new_solution_y_ptr[0]);
            this->solution.push_back(new_solution_y_ptr[1]);
            this->solution.push_back(new_solution_y_ptr[2]);
            this->solution.push_back(new_solution_y_ptr[3]);
            this->solution.push_back(new_solution_y_ptr[4]);
            this->solution.push_back(new_solution_y_ptr[5]);
            break;
        default:
            this->solution.insert(this->solution.end(), new_solution_y_ptr, new_solution_y_ptr + this->num_y);
            break;
    }

    // Repeak for any extra information that is captured.
    switch (this->num_extra)
    {
        case 0:
            // Not capturing extra. do nothing.
            break;
        case 1:
            this->solution.push_back(new_solution_dy_ptr[this->num_y]);
            break;
        case 2:
            this->solution.push_back(new_solution_dy_ptr[this->num_y    ]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 1]);
            break;
        case 3:
            this->solution.push_back(new_solution_dy_ptr[this->num_y    ]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 1]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 2]);
            break;
        case 4:
            this->solution.push_back(new_solution_dy_ptr[this->num_y    ]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 1]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 2]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 3]);
            break;
        case 5:
            this->solution.push_back(new_solution_dy_ptr[this->num_y    ]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 1]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 2]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 3]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 4]);
            break;
        case 6:
            this->solution.push_back(new_solution_dy_ptr[this->num_y    ]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 1]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 2]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 3]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 4]);
            this->solution.push_back(new_solution_dy_ptr[this->num_y + 5]);
            break;
        default:
            this->solution.insert(this->solution.end(), new_solution_dy_ptr[this->num_y], new_solution_dy_ptr[this->num_y] + this->num_extra);
            break;
    }
}

void CySolverResult::update_message(const char* const new_message_ptr)
{
    std::strcpy(this->message_ptr, new_message_ptr);
}
