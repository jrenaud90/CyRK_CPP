#pragma once


#include "common.hpp"
#include "cysolver.hpp"
#include "dense.hpp"

// ########################################################################################################################
// Dense Output Implementations
// ########################################################################################################################
class RKDenseOutput : public CySolverDense
{
protected:
    double step = 0.0;

    // Q is defined by Q = K.T.dot(self.P)  K has shape of (n_stages + 1, num_y) so K.T has shape of (num_y, n_stages + 1)
    // P has shape of (4, 3) for RK23; (7, 4) for RK45.. So (n_stages + 1, Q_order)
    // So Q has shape of (num_y, n_stages + 1)
    // Let's change it to (n_stages + 1, num_y)
    // The max size of Q is (7-3) * 16 = 64. Lets assume this is the max size and stack allocate Q.
    double Q[64] = { };
    unsigned int Q_order = 0;

public:
    double* Q_ptr = &Q[0];

protected:

public:
    virtual ~RKDenseOutput() {};
    RKDenseOutput() : CySolverDense() {};
    RKDenseOutput(double t_old, double t_now, double* y_in_ptr, unsigned int num_y, unsigned int Q_order);
    virtual void call(double t_interp, double* y_interped) override;
};


// ########################################################################################################################
// RK Integrators
// ########################################################################################################################

// #####################################################################################################################
// Runge - Kutta 2(3)
// #####################################################################################################################
const unsigned int RK23_METHOD_INT            = 0;
const unsigned int RK23_order                 = 3;
const unsigned int RK23_n_stages              = 3;
const unsigned int RK23_len_Arows             = 3;
const unsigned int RK23_len_Acols             = 3;
const unsigned int RK23_len_C                 = 3;
const unsigned int RK23_len_Pcols             = 3;
const unsigned int RK23_error_estimator_order = 2;
const double RK23_error_exponent = 1.0 / (2.0 + 1.0);  // Defined as 1 / (error_order + 1)

const double RK23_A[9] = {
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
const double* const RK23_A_ptr = &RK23_A[0];

const double RK23_B[3] = {
    2.0 / 9.0,
    1.0 / 3.0,
    4.0 / 9.0
};
const double* const RK23_B_ptr = &RK23_B[0];

const double RK23_C[3] = {
    0.0,
    1.0 / 2.0,
    3.0 / 4.0
};
const double* const RK23_C_ptr = &RK23_C[0];

const double RK23_E[4] = {
    5.0 / 72.0,
    -1.0 / 12.0,
    -1.0 / 9.0,
    1.0 / 8.0
};
const double* const RK23_E_ptr = &RK23_E[0];

// P is the transpose of the one scipy uses.
const double RK23_P[12] = {
    // Column 1
    1.0,
    0.0,
    0.0,
    0.0,

    // Column 2
    -4.0 / 3.0,
    1.0,
    4.0 / 3.0,
    -1.0,

    // Column 3
    5.0 / 9.0,
    -2.0 / 3.0,
    -8.0 / 9.0,
    1.0
};
const double* const RK23_P_ptr = &RK23_P[0];

// #####################################################################################################################
// Runge - Kutta 4(5)
// #####################################################################################################################
const unsigned int RK45_METHOD_INT            = 1;
const unsigned int RK45_order                 = 5;
const unsigned int RK45_n_stages              = 6;
const unsigned int RK45_len_Arows             = 6;
const unsigned int RK45_len_Acols             = 5;
const unsigned int RK45_len_C                 = 6;
const unsigned int RK45_len_Pcols             = 4;
const unsigned int RK45_error_estimator_order = 4;
const double RK45_error_exponent = 1.0 / (4.0 + 1.0);  // Defined as 1 / (error_order + 1)

const double RK45_A[30] = {
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
const double* const RK45_A_ptr = &RK45_A[0];

const double RK45_B[6] = {
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0
};
const double* const RK45_B_ptr = &RK45_B[0];

const double RK45_C[6] = {
    0.0,
    1.0 / 5.0,
    3.0 / 10.0,
    4.0 / 5.0,
    8.0 / 9.0,
    1.0
};
const double* const RK45_C_ptr = &RK45_C[0];

const double RK45_E[7] = {
    -71.0 / 57600.0,
    0.0,
    71.0 / 16695.0,
    -71.0 / 1920.0,
    17253.0 / 339200.0,
    -22.0 / 525.0,
    1.0 / 40.0
};
const double* const RK45_E_ptr = &RK45_E[0];


// P is the transpose of the one scipy uses.
const double RK45_P[28] = {
    // Column 1
    1.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,

    // Column 2
    -8048581381.0 / 2820520608.0,
    0.0,
    131558114200.0 / 32700410799.0,
    -1754552775.0 / 470086768.0,
    127303824393.0 / 49829197408.0,
    -282668133.0 / 205662961.0,
    40617522.0 / 29380423.0,

    // Column 3
    8663915743.0 / 2820520608.0,
    0.0,
    -68118460800.0 / 10900136933.0,
    14199869525.0 / 1410260304.0,
    -318862633887.0 / 49829197408.0,
    2019193451.0 / 616988883.0,
    -110615467.0 / 29380423.0,

    // Column 4
    -12715105075.0 / 11282082432.0,
    0.0,
    87487479700.0 / 32700410799.0,
    -10690763975.0 / 1880347072.0,
    701980252875.0 / 199316789632.0,
    -1453857185.0 / 822651844.0,
    69997945.0 / 29380423.0
};
const double* const RK45_P_ptr = &RK45_P[0];


// #####################################################################################################################
// Runge - Kutta DOP 8(5; 3)
// #####################################################################################################################
const unsigned int DOP853_METHOD_INT            = 2;
const unsigned int DOP853_order                 = 8;
const unsigned int DOP853_n_stages              = 12;
const unsigned int DOP853_A_rows                = 12;
const unsigned int DOP853_A_cols                = 12;
const unsigned int DOP853_len_C                 = 12;
const unsigned int DOP853_error_estimator_order = 7;
const double DOP853_error_exponent              = 1.0 / (7.0 + 1.0);  // Defined as 1 / (error_order + 1)

// Note both A and C are the _reduced_ versions.The full A and C are not shown.
const double DOP853_A[144] = {
    // A - Row 0
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    // A - Row 1
    5.26001519587677318785587544488e-2,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    // A - Row 2
    1.97250569845378994544595329183e-2,
    5.91751709536136983633785987549e-2,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    // A - Row 3
    2.95875854768068491816892993775e-2,
    0.0,
    8.87627564304205475450678981324e-2,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    // A - Row 4
    2.41365134159266685502369798665e-1,
    0.0,
    -8.84549479328286085344864962717e-1,
    9.24834003261792003115737966543e-1,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    // A - Row 5
    3.7037037037037037037037037037e-2,
    0.0,
    0.0,
    1.70828608729473871279604482173e-1,
    1.25467687566822425016691814123e-1,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0, // Nice
    0.0,
    0.0,
    // A - Row 6
    3.7109375e-2,
    0.0,
    0.0,
    1.70252211019544039314978060272e-1,
    6.02165389804559606850219397283e-2,
    -1.7578125e-2,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    // A - Row 7
    3.70920001185047927108779319836e-2,
    0.0,
    0.0,
    1.70383925712239993810214054705e-1,
    1.07262030446373284651809199168e-1,
    -1.53194377486244017527936158236e-2,
    8.27378916381402288758473766002e-3,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    // A - Row 8
    6.24110958716075717114429577812e-1,
    0.0,
    0.0,
    -3.36089262944694129406857109825,
    -8.68219346841726006818189891453e-1,
    2.75920996994467083049415600797e1,
    2.01540675504778934086186788979e1,
    -4.34898841810699588477366255144e1,
    0.0,
    0.0,
    0.0,
    0.0,
    // A - Row 9
    4.77662536438264365890433908527e-1,
    0.0,
    0.0,
    -2.48811461997166764192642586468,
    -5.90290826836842996371446475743e-1,
    2.12300514481811942347288949897e1,
    1.52792336328824235832596922938e1,
    -3.32882109689848629194453265587e1,
    -2.03312017085086261358222928593e-2,
    0.0,
    0.0,
    0.0,
    // A - Row 10
    -9.3714243008598732571704021658e-1,
    0.0,
    0.0,
    5.18637242884406370830023853209,
    1.09143734899672957818500254654,
    -8.14978701074692612513997267357,
    -1.85200656599969598641566180701e1,
    2.27394870993505042818970056734e1,
    2.49360555267965238987089396762,
    -3.0467644718982195003823669022,
    0.0,
    0.0,
    // A - Row 11
    2.27331014751653820792359768449,
    0.0,
    0.0,
    -1.05344954667372501984066689879e1,
    -2.00087205822486249909675718444,
    -1.79589318631187989172765950534e1,
    2.79488845294199600508499808837e1,
    -2.85899827713502369474065508674,
    -8.87285693353062954433549289258,
    1.23605671757943030647266201528e1,
    6.43392746015763530355970484046e-1,
    0.0
};
const double* const DOP853_A_ptr = &DOP853_A[0];

// Note: B is equal to the 13th row of the expanded version of A(which we do not define above)
const double DOP853_B[12] = {
    5.42937341165687622380535766363e-2,
    0.0,
    0.0,
    0.0,
    0.0,
    4.45031289275240888144113950566,
    1.89151789931450038304281599044,
    -5.8012039600105847814672114227,
    3.1116436695781989440891606237e-1,
    -1.52160949662516078556178806805e-1,
    2.01365400804030348374776537501e-1,
    4.47106157277725905176885569043e-2
};
const double* const DOP853_B_ptr = &DOP853_B[0];


// Note this is the reduced C array.The expanded version is not shown.
const double DOP853_C[12] = {
    0.0,
    0.526001519587677318785587544488e-01,
    0.789002279381515978178381316732e-01,
    0.118350341907227396726757197510,
    0.281649658092772603273242802490,
    0.333333333333333333333333333333,
    0.25,
    0.307692307692307692307692307692,
    0.651282051282051282051282051282,
    0.6,
    0.857142857142857142857142857142,
    1.0
};
const double* const DOP853_C_ptr = &DOP853_C[0];

// All except last value equals B(B length is one less than E3).
const double DOP853_E3[13] = {
    5.42937341165687622380535766363e-2 - 0.244094488188976377952755905512,
    0.0,
    0.0,
    0.0,
    0.0,
    4.45031289275240888144113950566,
    1.89151789931450038304281599044,
    -5.8012039600105847814672114227,
    3.1116436695781989440891606237e-1 - 0.733846688281611857341361741547,
    -1.52160949662516078556178806805e-1,
    2.01365400804030348374776537501e-1,
    4.47106157277725905176885569043e-2 - 0.220588235294117647058823529412e-1,
    0.0
};
const double* const DOP853_E3_ptr = &DOP853_E3[0];

const double DOP853_E5[13] = {
    0.1312004499419488073250102996e-1,
    0.0,
    0.0,
    0.0,
    0.0,
    -0.1225156446376204440720569753e+1,
    -0.4957589496572501915214079952,
    0.1664377182454986536961530415e+1,
    -0.3503288487499736816886487290,
    0.3341791187130174790297318841,
    0.8192320648511571246570742613e-1,
    -0.2235530786388629525884427845e-1,
    0.
};
const double* const DOP853_E5_ptr = &DOP853_E5[0];

// ########################################################################################################################
// Classes
// ########################################################################################################################
class RKSolver : public CySolverBase {

// Attributes
protected:
    // Step globals
    const double error_safety    = SAFETY;
    const double min_step_factor = MIN_FACTOR;
    const double max_step_factor = MAX_FACTOR;

    // RK constants
    unsigned int order = 0;
    unsigned int error_estimator_order = 0;
    unsigned int n_stages     = 0;
    unsigned int n_stages_p1  = 0;
    unsigned int len_Acols    = 0;
    unsigned int len_C        = 0;
    unsigned int len_Pcols    = 0;
    unsigned int nstages_numy = 0;
    double error_exponent     = 0.0;

    // Pointers to RK constant arrays
    const double* C_ptr  = nullptr;
    const double* A_ptr  = nullptr;
    const double* B_ptr  = nullptr;
    const double* E_ptr  = nullptr;
    const double* E3_ptr = nullptr;
    const double* E5_ptr = nullptr;
    const double* P_ptr  = nullptr;
    const double* D_ptr  = nullptr;

    // K is not const. Its values are stored in an array that is held by this class.
    double K[1]   = { std::nan("") };
    double* K_ptr = &this->K[0];

    // Tolerances
    // For the same reason num_y is limited, the total number of tolerances are limited.
    double rtols[Y_LIMIT] = { std::nan("") };
    double atols[Y_LIMIT] = { std::nan("") };
    double* rtols_ptr     = &rtols[0];
    double* atols_ptr     = &atols[0];
    bool use_array_rtols  = false;
    bool use_array_atols  = false;

    // Step size parameters
    double user_provided_first_step_size = 0.0;
    double step          = 0.0;
    double step_size     = 0.0;
    double max_step_size = 0.0;

    // Error estimate
    double error_norm = 0.0;


// Methods
protected:
    virtual void p_estimate_error() override;
    virtual void p_step_implementation() override;
    virtual CySolverDense* p_dense_output() override;

public:
    RKSolver();
    virtual ~RKSolver() override;
    RKSolver(
        // Base Class input arguments
        DiffeqFuncType diffeq_ptr,
        std::shared_ptr<CySolverResult> const storage_ptr,
        const double t_start,
        const double t_end,
        const double* y0_ptr,
        const unsigned int num_y,
        const unsigned int num_extra = 0,
        const double* args_ptr = nullptr,
        const size_t max_num_steps = 0,
        const size_t max_ram_MB = 2000,
        const bool dense_output = false,
        const double* t_eval = nullptr,
        const size_t len_t_eval = 0,
        // RKSolver input arguments
        const double rtol = 1.0e-3,
        const double atol = 1.0e-6,
        const double* rtols_ptr = nullptr,
        const double* atols_ptr = nullptr,
        const double max_step_size = MAX_STEP,
        const double first_step_size = 0.0
    );
    virtual void reset() override;
    virtual void calc_first_step_size() override;
};



class RK23 : public RKSolver {

protected:
    double K[4 * Y_LIMIT] = { 0.0 };

public:
    // Copy over base class constructors
    using RKSolver::RKSolver;
    virtual void reset() override;
};

class RK45 : public RKSolver {

protected:
    double K[7 * Y_LIMIT] = { 0.0 };

public:
    // Copy over base class constructors
    using RKSolver::RKSolver;
    virtual void reset() override;
};

class DOP853 : public RKSolver {

protected:
    double K[13 * Y_LIMIT] = { 0.0 };

public:
    // Copy over base class constructors
    using RKSolver::RKSolver;
    virtual void reset() override;
    virtual void p_estimate_error() override;
};
