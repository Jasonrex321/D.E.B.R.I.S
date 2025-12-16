/**
 * @file state_vector.c
 * @brief Implementation of core state vector operations
 * 
 * @see state_vector.h for documentation
 */

#include "core/state_vector.h"
#include "core/error_handling.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* ============================================================================
 * PRIVATE CONSTANTS & HELPER FUNCTIONS
 * ========================================================================== */

/** Earth radius at poles (WGS84) in km */
static const double EARTH_RADIUS_POLAR = 6356.752314245;

/** Small value to avoid division by zero */
static const double EPS = 1e-15;

/** Maximum iterations for singularity handling */
static const int MAX_SINGULARITY_ITERATIONS = 10;

/**
 * @brief Kahan-Babushka-Neumaier summation for 3 elements
 * 
 * Compensated summation to minimize round-off error.
 * 
 * @param x Array of 3 doubles
 * @return double Compensated sum
 */
static double kbn_sum_3(const double x[3]) {
    double sum = 0.0;
    double c = 0.0;  // Compensation term
    
    for (int i = 0; i < 3; i++) {
        double y = x[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    
    return sum;
}

/**
 * @brief Safe inverse sine with domain checking
 * 
 * @param x Input in [-1, 1]
 * @return double arcsin(x) with domain protection
 */
static double safe_asin(double x) {
    if (x >= 1.0 - FLOAT_EPSILON) {
        return M_PI_2;
    } else if (x <= -1.0 + FLOAT_EPSILON) {
        return -M_PI_2;
    } else {
        return asin(x);
    }
}

/**
 * @brief Safe inverse cosine with domain checking
 * 
 * @param x Input in [-1, 1]
 * @return double arccos(x) with domain protection
 */
static double safe_acos(double x) {
    if (x >= 1.0 - FLOAT_EPSILON) {
        return 0.0;
    } else if (x <= -1.0 + FLOAT_EPSILON) {
        return M_PI;
    } else {
        return acos(x);
    }
}

/**
 * @brief Wrap angle to [0, 2π)
 * 
 * @param angle Angle in radians
 * @return double Wrapped angle
 */
static double wrap_to_2pi(double angle) {
    angle = fmod(angle, 2.0 * M_PI);
    if (angle < 0.0) {
        angle += 2.0 * M_PI;
    }
    return angle;
}

/**
 * @brief Wrap angle to [-π, π)
 * 
 * @param angle Angle in radians
 * @return double Wrapped angle
 */
static double wrap_to_pi(double angle) {
    angle = fmod(angle + M_PI, 2.0 * M_PI);
    if (angle < 0.0) {
        angle += 2.0 * M_PI;
    }
    return angle - M_PI;
}

/* ============================================================================
 * LIFE CYCLE IMPLEMENTATIONS
 * ========================================================================== */

StateVector_t* create_state_vector(double epoch_jd, uint64_t object_id) {
    // Allocate with calloc for zero-initialization
    StateVector_t* state = (StateVector_t*)calloc(1, sizeof(StateVector_t));
    if (state == NULL) {
        return NULL;
    }
    
    state->epoch_jd = epoch_jd;
    state->object_id = object_id;
    state->requires_cleanup = true;
    state->reflectance = 1.0;  // Default reflectance
    
    // Position and velocity are already zeroed by calloc
    // mass_kg and cross_section_m2 are 0.0
    
    return state;
}

StateVector_t* create_state_vector_from_arrays(double epoch_jd, 
                                              const double pos_km[3],
                                              const double vel_kms[3],
                                              uint64_t object_id) {
    StateVector_t* state = create_state_vector(epoch_jd, object_id);
    if (state == NULL) {
        return NULL;
    }
    
    memcpy(state->position_km, pos_km, 3 * sizeof(double));
    memcpy(state->velocity_kms, vel_kms, 3 * sizeof(double));
    
    return state;
}

void destroy_state_vector(StateVector_t** state_pointer) {
    if (state_pointer == NULL || *state_pointer == NULL) {
        return;
    }
    
    StateVector_t* state = *state_pointer;
    
    // Log for debugging (compile-time option)
#ifdef DEBUG_MEMORY
    printf("[DEBUG] Freeing StateVector at %p (ID: %lu)\n", 
           (void*)state, state->object_id);
#endif
    
    // Free the state
    free(state);
    
    // NULLify to prevent dangling pointer
    *state_pointer = NULL;
}

StateVector_t* copy_state_vector(const StateVector_t* source) {
    if (source == NULL) {
        return NULL;
    }
    
    StateVector_t* copy = create_state_vector(source->epoch_jd, source->object_id);
    if (copy == NULL) {
        return NULL;
    }
    
    // Copy all fields
    memcpy(copy->position_km, source->position_km, 3 * sizeof(double));
    memcpy(copy->velocity_kms, source->velocity_kms, 3 * sizeof(double));
    copy->mass_kg = source->mass_kg;
    copy->cross_section_m2 = source->cross_section_m2;
    copy->reflectance = source->reflectance;
    copy->is_debris = source->is_debris;
    
    return copy;
}

Result_t initialize_state_vector(StateVector_t* state, 
                                double epoch_jd, 
                                uint64_t object_id) {
    if (state == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    state->epoch_jd = epoch_jd;
    state->object_id = object_id;
    
    // Zero all arrays and fields
    for (int i = 0; i < 3; i++) {
        state->position_km[i] = 0.0;
        state->velocity_kms[i] = 0.0;
    }
    
    state->mass_kg = 0.0;
    state->cross_section_m2 = 0.0;
    state->reflectance = 1.0;
    state->requires_cleanup = false;  // Stack-allocated
    state->is_debris = false;
    
    return RESULT_OK;
}

/* ============================================================================
 * VECTOR OPERATIONS IMPLEMENTATIONS (WITH KBN)
 * ========================================================================== */

double vector_norm_kbn(const double vec[3]) {
    if (vec == NULL) {
        return NAN;
    }
    
    // Compute squared components
    double squares[3];
    for (int i = 0; i < 3; i++) {
        squares[i] = vec[i] * vec[i];
    }
    
    // Sum with KBN
    double sum = kbn_sum_3(squares);
    
    // Return sqrt (safe for zero)
    return sqrt(fmax(sum, 0.0));
}

double vector_dot_product_kbn(const double a[3], const double b[3]) {
    if (a == NULL || b == NULL) {
        return NAN;
    }
    
    // Compute component-wise products
    double products[3];
    for (int i = 0; i < 3; i++) {
        products[i] = a[i] * b[i];
    }
    
    // Sum with KBN
    return kbn_sum_3(products);
}

void vector_cross_product(const double a[3], const double b[3], double result[3]) {
    if (a == NULL || b == NULL || result == NULL) {
        return;
    }
    
    // Check for aliasing
    if (result == a || result == b) {
        // Need temporary storage
        double temp[3];
        temp[0] = a[1] * b[2] - a[2] * b[1];
        temp[1] = a[2] * b[0] - a[0] * b[2];
        temp[2] = a[0] * b[1] - a[1] * b[0];
        
        memcpy(result, temp, 3 * sizeof(double));
    } else {
        result[0] = a[1] * b[2] - a[2] * b[1];
        result[1] = a[2] * b[0] - a[0] * b[2];
        result[2] = a[0] * b[1] - a[1] * b[0];
    }
}

void vector_add_scaled(double result[3], 
                      const double a[3], 
                      double scale, 
                      const double b[3]) {
    if (result == NULL || a == NULL || b == NULL) {
        return;
    }
    
    // Handle aliasing cases
    if (result == a) {
        // result is a, so update in-place
        result[0] += scale * b[0];
        result[1] += scale * b[1];
        result[2] += scale * b[2];
    } else if (result == b) {
        // result is b, need temporary
        double temp[3];
        temp[0] = a[0] + scale * b[0];
        temp[1] = a[1] + scale * b[1];
        temp[2] = a[2] + scale * b[2];
        
        memcpy(result, temp, 3 * sizeof(double));
    } else {
        // No aliasing, simple case
        result[0] = a[0] + scale * b[0];
        result[1] = a[1] + scale * b[1];
        result[2] = a[2] + scale * b[2];
    }
}

double vector_normalize(double vec[3]) {
    if (vec == NULL) {
        return 0.0;
    }
    
    double norm = vector_norm_kbn(vec);
    
    if (norm > EPS) {
        double inv_norm = 1.0 / norm;
        vec[0] *= inv_norm;
        vec[1] *= inv_norm;
        vec[2] *= inv_norm;
    } else {
        // Zero vector, set to zero (don't crash)
        vec[0] = vec[1] = vec[2] = 0.0;
    }
    
    return norm;
}

double vector_angle_between(const double a[3], const double b[3]) {
    if (a == NULL || b == NULL) {
        return NAN;
    }
    
    double norm_a = vector_norm_kbn(a);
    double norm_b = vector_norm_kbn(b);
    
    if (norm_a < EPS || norm_b < EPS) {
        return 0.0;  // Zero vector has undefined angle, return 0
    }
    
    double dot = vector_dot_product_kbn(a, b);
    double cos_angle = dot / (norm_a * norm_b);
    
    // Clamp to valid range for acos
    cos_angle = fmax(fmin(cos_angle, 1.0), -1.0);
    
    return safe_acos(cos_angle);
}

/* ============================================================================
 * PHYSICAL PROPERTIES IMPLEMENTATIONS
 * ========================================================================== */

double get_altitude_km(const StateVector_t* state) {
    if (state == NULL) {
        return NAN;
    }
    
    // Compute geocentric distance
    double r = vector_norm_kbn(state->position_km);
    
    // Simple spherical Earth approximation for now
    // TODO: Implement WGS84 ellipsoid altitude
    return r - EARTH_RADIUS_EQUATORIAL;
}

double compute_specific_energy(const StateVector_t* state) {
    if (state == NULL) {
        return NAN;
    }
    
    double v_sq = vector_dot_product_kbn(state->velocity_kms, state->velocity_kms);
    double r = vector_norm_kbn(state->position_km);
    
    if (r < EPS) {
        return INFINITY;  // At Earth center - singular
    }
    
    return 0.5 * v_sq - EARTH_MU / r;
}

Result_t compute_angular_momentum_vector(const StateVector_t* state, double h_vec[3]) {
    if (state == NULL || h_vec == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    vector_cross_product(state->position_km, state->velocity_kms, h_vec);
    return RESULT_OK;
}

double compute_angular_momentum_magnitude(const StateVector_t* state) {
    if (state == NULL) {
        return NAN;
    }
    
    double h_vec[3];
    if (compute_angular_momentum_vector(state, h_vec) != RESULT_OK) {
        return NAN;
    }
    
    return vector_norm_kbn(h_vec);
}

double compute_orbital_period(const StateVector_t* state) {
    if (state == NULL) {
        return NAN;
    }
    
    double energy = compute_specific_energy(state);
    
    // Only elliptical orbits have finite period
    if (energy >= 0.0) {
        return INFINITY;  // Parabolic or hyperbolic
    }
    
    double a = -EARTH_MU / (2.0 * energy);
    if (a <= 0.0) {
        return INFINITY;  // Invalid
    }
    
    return 2.0 * M_PI * sqrt(a * a * a / EARTH_MU);
}

double compute_semi_major_axis(const StateVector_t* state) {
    if (state == NULL) {
        return NAN;
    }
    
    double energy = compute_specific_energy(state);
    
    if (fabs(energy) < EPS) {
        return INFINITY;  // Parabolic
    }
    
    return -EARTH_MU / (2.0 * energy);
}

Result_t compute_eccentricity_vector(const StateVector_t* state, double e_vec[3]) {
    if (state == NULL || e_vec == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    // Compute angular momentum
    double h_vec[3];
    compute_angular_momentum_vector(state, h_vec);
    
    // Compute v × h
    double v_cross_h[3];
    vector_cross_product(state->velocity_kms, h_vec, v_cross_h);
    
    // First term: (v × h)/μ
    for (int i = 0; i < 3; i++) {
        e_vec[i] = v_cross_h[i] / EARTH_MU;
    }
    
    // Second term: r/|r|
    double r = vector_norm_kbn(state->position_km);
    
    if (r < EPS) {
        return RESULT_NUMERICAL_ERROR;  // Division by zero
    }
    
    double inv_r = 1.0 / r;
    for (int i = 0; i < 3; i++) {
        e_vec[i] -= state->position_km[i] * inv_r;
    }
    
    return RESULT_OK;
}

double compute_eccentricity(const StateVector_t* state) {
    if (state == NULL) {
        return NAN;
    }
    
    double e_vec[3];
    if (compute_eccentricity_vector(state, e_vec) != RESULT_OK) {
        return NAN;
    }
    
    return vector_norm_kbn(e_vec);
}

/* ============================================================================
 * COORDINATE TRANSFORMATIONS IMPLEMENTATIONS
 * ========================================================================== */

Result_t state_vector_to_elements(const StateVector_t* state, 
                                 KeplerianElements_t* elements) {
    if (state == NULL || elements == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    // Validate input state
    Result_t valid = validate_state_vector(state);
    if (valid != RESULT_OK) {
        return valid;
    }
    
    // Compute fundamental vectors
    double r_vec[3];
    memcpy(r_vec, state->position_km, 3 * sizeof(double));
    
    double v_vec[3];
    memcpy(v_vec, state->velocity_kms, 3 * sizeof(double));
    
    double r = vector_norm_kbn(r_vec);
    double v = vector_norm_kbn(v_vec);
    
    // Compute angular momentum vector
    double h_vec[3];
    vector_cross_product(r_vec, v_vec, h_vec);
    double h = vector_norm_kbn(h_vec);
    
    // Check for degenerate orbit (collinear position and velocity)
    if (h < EPS) {
        return RESULT_SINGULAR_MATRIX;  // Radial orbit
    }
    
    // Compute eccentricity vector
    double e_vec[3];
    Result_t e_result = compute_eccentricity_vector(state, e_vec);
    if (e_result != RESULT_OK) {
        return e_result;
    }
    
    double e = vector_norm_kbn(e_vec);
    
    // Compute semi-major axis
    double energy = 0.5 * v * v - EARTH_MU / r;
    double a = -EARTH_MU / (2.0 * energy);
    
    // Compute inclination
    double hz = h_vec[2];
    double cos_i = hz / h;
    
    // Clamp to valid range for acos
    cos_i = fmax(fmin(cos_i, 1.0), -1.0);
    double i = safe_acos(cos_i);
    
    // Compute RAAN (Ω)
    double nx = -h_vec[1];
    double ny = h_vec[0];
    double n = sqrt(nx * nx + ny * ny);
    
    double raan;
    if (n < EPS) {
        // Equatorial orbit, RAAN undefined, set to 0
        raan = 0.0;
    } else {
        raan = atan2(ny, nx);
        raan = wrap_to_2pi(raan);
    }
    
    // Compute argument of perigee (ω)
    double arg_perigee;
    if (e < EPS) {
        // Circular orbit, argument of perigee undefined, set to 0
        arg_perigee = 0.0;
    } else {
        double ne_dot = nx * e_vec[0] + ny * e_vec[1];
        double e_cos_omega = ne_dot / (n * e);
        
        // Clamp to valid range
        e_cos_omega = fmax(fmin(e_cos_omega, 1.0), -1.0);
        
        // Check if orbit is prograde or retrograde
        double ez = e_vec[2];
        double e_sin_omega = ez / (e * sin(i + EPS));
        
        arg_perigee = atan2(e_sin_omega, e_cos_omega);
        arg_perigee = wrap_to_2pi(arg_perigee);
    }
    
    // Compute true anomaly (ν)
    double true_anomaly;
    if (e < EPS) {
        // Circular orbit, true anomaly = argument of latitude
        double sin_u = r_vec[2] / (r * sin(i + EPS));
        double cos_u = (r_vec[0] * cos(raan) + r_vec[1] * sin(raan)) / r;
        true_anomaly = atan2(sin_u, cos_u);
        true_anomaly = wrap_to_2pi(true_anomaly);
    } else {
        double r_dot_v = vector_dot_product_kbn(r_vec, v_vec);
        double sin_nu = r_dot_v / (e * sqrt(EARTH_MU * a + EPS));
        
        double cos_nu = (a * (1.0 - e * e) - r) / (e * r + EPS);
        
        true_anomaly = atan2(sin_nu, cos_nu);
        true_anomaly = wrap_to_2pi(true_anomaly);
    }
    
    // Populate output structure
    elements->semi_major_axis = a;
    elements->eccentricity = e;
    elements->inclination = i;
    elements->raan = raan;
    elements->arg_perigee = arg_perigee;
    elements->true_anomaly = true_anomaly;
    elements->epoch_jd = state->epoch_jd;
    
    return RESULT_OK;
}

Result_t elements_to_state_vector(const KeplerianElements_t* elements,
                                 StateVector_t* state) {
    if (elements == NULL || state == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    // Validate elements
    if (elements->semi_major_axis <= 0.0) {
        return RESULT_INVALID_ORBIT;
    }
    if (elements->eccentricity < 0.0) {
        return RESULT_INVALID_ARGUMENT;
    }
    
    double a = elements->semi_major_axis;
    double e = elements->eccentricity;
    double i = elements->inclination;
    double raan = elements->raan;
    double arg_perigee = elements->arg_perigee;
    double nu = elements->true_anomaly;
    
    // Compute parameter p (semi-latus rectum)
    double p = a * (1.0 - e * e);
    
    if (p <= 0.0) {
        return RESULT_INVALID_ORBIT;
    }
    
    // Compute radius at true anomaly
    double r = p / (1.0 + e * cos(nu));
    
    // Position in perifocal frame (PQW)
    double pos_pqw[3];
    pos_pqw[0] = r * cos(nu);
    pos_pqw[1] = r * sin(nu);
    pos_pqw[2] = 0.0;
    
    // Velocity in perifocal frame
    double sqrt_mu_p = sqrt(EARTH_MU / p);
    double vel_pqw[3];
    vel_pqw[0] = -sqrt_mu_p * sin(nu);
    vel_pqw[1] = sqrt_mu_p * (e + cos(nu));
    vel_pqw[2] = 0.0;
    
    // Rotation matrices
    double cos_raan = cos(raan);
    double sin_raan = sin(raan);
    double cos_i = cos(i);
    double sin_i = sin(i);
    double cos_argp = cos(arg_perigee);
    double sin_argp = sin(arg_perigee);
    
    // Combined rotation matrix R = R_z(Ω) * R_x(i) * R_z(ω)
    double R[3][3];
    
    // Column 1
    R[0][0] = cos_raan * cos_argp - sin_raan * sin_argp * cos_i;
    R[1][0] = sin_raan * cos_argp + cos_raan * sin_argp * cos_i;
    R[2][0] = sin_argp * sin_i;
    
    // Column 2
    R[0][1] = -cos_raan * sin_argp - sin_raan * cos_argp * cos_i;
    R[1][1] = -sin_raan * sin_argp + cos_raan * cos_argp * cos_i;
    R[2][1] = cos_argp * sin_i;
    
    // Column 3
    R[0][2] = sin_raan * sin_i;
    R[1][2] = -cos_raan * sin_i;
    R[2][2] = cos_i;
    
    // Transform from perifocal to ECI
    for (int row = 0; row < 3; row++) {
        state->position_km[row] = 0.0;
        state->velocity_kms[row] = 0.0;
        
        for (int col = 0; col < 3; col++) {
            state->position_km[row] += R[row][col] * pos_pqw[col];
            state->velocity_kms[row] += R[row][col] * vel_pqw[col];
        }
    }
    
    // Set other fields
    state->epoch_jd = elements->epoch_jd;
    state->object_id = 0;  // Unknown object ID from elements alone
    
    return RESULT_OK;
}

Result_t state_vector_to_equinoctial(const StateVector_t* state,
                                    EquinoctialElements_t* equinoctial) {
    if (state == NULL || equinoctial == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    // First get Keplerian elements
    KeplerianElements_t keplerian;
    Result_t result = state_vector_to_elements(state, &keplerian);
    if (result != RESULT_OK) {
        return result;
    }
    
    double a = keplerian.semi_major_axis;
    double e = keplerian.eccentricity;
    double i = keplerian.inclination;
    double omega = keplerian.arg_perigee;
    double raan = keplerian.raan;
    double nu = keplerian.true_anomaly;
    
    // Compute equinoctial elements
    equinoctial->p = a * (1.0 - e * e);
    equinoctial->f = e * cos(omega + raan);
    equinoctial->g = e * sin(omega + raan);
    
    double tan_i_over_2 = tan(i / 2.0);
    equinoctial->h = tan_i_over_2 * cos(raan);
    equinoctial->k = tan_i_over_2 * sin(raan);
    
    // True longitude: L = Ω + ω + ν
    equinoctial->L = wrap_to_2pi(raan + omega + nu);
    equinoctial->epoch_jd = state->epoch_jd;
    
    return RESULT_OK;
}

Result_t compute_rtn_frame(const StateVector_t* state, RTNFrame_t* rtn) {
    if (state == NULL || rtn == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    // Radial direction: R = r/|r|
    double r = vector_norm_kbn(state->position_km);
    if (r < EPS) {
        return RESULT_NUMERICAL_ERROR;
    }
    
    for (int i = 0; i < 3; i++) {
        rtn->radial[i] = state->position_km[i] / r;
    }
    
    // Normal direction: N = (r × v)/|r × v|
    double h_vec[3];
    vector_cross_product(state->position_km, state->velocity_kms, h_vec);
    double h = vector_norm_kbn(h_vec);
    
    if (h < EPS) {
        // Radial orbit, use arbitrary normal (z-axis)
        rtn->normal[0] = rtn->normal[1] = 0.0;
        rtn->normal[2] = 1.0;
    } else {
        for (int i = 0; i < 3; i++) {
            rtn->normal[i] = h_vec[i] / h;
        }
    }
    
    // Tangential direction: T = N × R
    vector_cross_product(rtn->normal, rtn->radial, rtn->tangential);
    
    return RESULT_OK;
}

Result_t transform_eci_to_rtn(const StateVector_t* state,
                             const double vec_eci[3],
                             double vec_rtn[3]) {
    if (state == NULL || vec_eci == NULL || vec_rtn == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    RTNFrame_t rtn;
    Result_t result = compute_rtn_frame(state, &rtn);
    if (result != RESULT_OK) {
        return result;
    }
    
    // Transform by dotting with each axis
    vec_rtn[0] = vector_dot_product_kbn(vec_eci, rtn.radial);
    vec_rtn[1] = vector_dot_product_kbn(vec_eci, rtn.tangential);
    vec_rtn[2] = vector_dot_product_kbn(vec_eci, rtn.normal);
    
    return RESULT_OK;
}

Result_t transform_rtn_to_eci(const StateVector_t* state,
                             const double vec_rtn[3],
                             double vec_eci[3]) {
    if (state == NULL || vec_rtn == NULL || vec_eci == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    RTNFrame_t rtn;
    Result_t result = compute_rtn_frame(state, &rtn);
    if (result != RESULT_OK) {
        return result;
    }
    
    // Reconstruct: vec_eci = R*vec_rtn[0] + T*vec_rtn[1] + N*vec_rtn[2]
    for (int i = 0; i < 3; i++) {
        vec_eci[i] = rtn.radial[i] * vec_rtn[0] +
                    rtn.tangential[i] * vec_rtn[1] +
                    rtn.normal[i] * vec_rtn[2];
    }
    
    return RESULT_OK;
}

/* ============================================================================
 * VALIDATION & SAFETY IMPLEMENTATIONS
 * ========================================================================== */

Result_t validate_state_vector(const StateVector_t* state) {
    if (state == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    // Check for NaN/Inf
    if (check_for_nan_inf(state)) {
        return RESULT_NUMERICAL_ERROR;
    }
    
    // Check epoch
    if (!isfinite(state->epoch_jd)) {
        return RESULT_INVALID_ARGUMENT;
    }
    
    // Check position magnitude (not too large)
    double r = vector_norm_kbn(state->position_km);
    if (!isfinite(r) || r > 1e6) {  // More than 1 million km from Earth
        return RESULT_OUT_OF_RANGE;
    }
    
    // Check velocity magnitude (not too large)
    double v = vector_norm_kbn(state->velocity_kms);
    if (!isfinite(v) || v > 100.0) {  // More than 100 km/s
        return RESULT_OUT_OF_RANGE;
    }
    
    // Check if orbit is physically possible
    double energy = compute_specific_energy(state);
    if (!isfinite(energy)) {
        return RESULT_NUMERICAL_ERROR;
    }
    
    // Check altitude
    double altitude = get_altitude_km(state);
    if (altitude < -100.0) {  // 100 km below surface
        return RESULT_INVALID_ORBIT;
    }
    
    return RESULT_OK;
}

bool check_for_nan_inf(const StateVector_t* state) {
    if (state == NULL) {
        return true;
    }
    
    // Check all double fields
    if (!isfinite(state->epoch_jd)) return true;
    if (!isfinite(state->mass_kg)) return true;
    if (!isfinite(state->cross_section_m2)) return true;
    if (!isfinite(state->reflectance)) return true;
    
    // Check arrays
    for (int i = 0; i < 3; i++) {
        if (!isfinite(state->position_km[i])) return true;
        if (!isfinite(state->velocity_kms[i])) return true;
    }
    
    return false;
}

bool is_orbit_stable(const StateVector_t* state) {
    if (state == NULL) {
        return false;
    }
    
    // Check eccentricity < 1 (elliptical)
    double e = compute_eccentricity(state);
    if (e >= 1.0 - FLOAT_EPSILON) {
        return false;
    }
    
    // Check perigee altitude > 0
    double altitude = get_altitude_km(state);
    if (altitude < 0.0) {
        return false;
    }
    
    // Check not escaping
    double energy = compute_specific_energy(state);
    if (energy >= 0.0) {
        return false;
    }
    
    return true;
}

bool is_decaying_orbit(const StateVector_t* state, double threshold_km) {
    if (state == NULL) {
        return false;
    }
    
    // Compute perigee distance
    double a = compute_semi_major_axis(state);
    double e = compute_eccentricity(state);
    
    if (e >= 1.0) {
        return false;  // Hyperbolic, not decaying
    }
    
    double perigee = a * (1.0 - e);
    double perigee_altitude = perigee - EARTH_RADIUS_EQUATORIAL;
    
    return perigee_altitude < threshold_km;
}

bool is_escape_trajectory(const StateVector_t* state) {
    if (state == NULL) {
        return false;
    }
    
    double energy = compute_specific_energy(state);
    return energy >= 0.0;
}

bool is_circular_orbit(const StateVector_t* state, double tolerance) {
    if (state == NULL) {
        return false;
    }
    
    double e = compute_eccentricity(state);
    return e < tolerance;
}

/* ============================================================================
 * MEMORY ARENA IMPLEMENTATIONS
 * ========================================================================== */

Result_t initialize_particle_arena(ParticleArena_t* arena, size_t size_mb) {
    if (arena == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    if (size_mb == 0) {
        return RESULT_INVALID_ARGUMENT;
    }
    
    size_t size_bytes = size_mb * 1024 * 1024;
    
    // Allocate aligned memory for better performance
#ifdef _WIN32
    arena->memory_block = (uint8_t*)_aligned_malloc(size_bytes, 64);
#else
    if (posix_memalign((void**)&arena->memory_block, 64, size_bytes) != 0) {
        arena->memory_block = NULL;
    }
#endif
    
    if (arena->memory_block == NULL) {
        return RESULT_OUT_OF_MEMORY;
    }
    
    arena->block_size = size_bytes;
    arena->current_offset = 0;
    arena->allocation_count = 0;
    arena->high_water_mark = 0;
    
    // Zero-initialize the memory
    memset(arena->memory_block, 0, size_bytes);
    
    return RESULT_OK;
}

void* arena_allocate_particles(ParticleArena_t* arena, 
                              size_t particle_count, 
                              size_t particle_size) {
    if (arena == NULL || arena->memory_block == NULL) {
        return NULL;
    }
    
    if (particle_count == 0 || particle_size == 0) {
        return NULL;
    }
    
    size_t required_bytes = particle_count * particle_size;
    
    // Check bounds
    if (arena->current_offset + required_bytes > arena->block_size) {
        return NULL;
    }
    
    void* allocation = arena->memory_block + arena->current_offset;
    arena->current_offset += required_bytes;
    arena->allocation_count++;
    
    // Update high water mark
    if (arena->current_offset > arena->high_water_mark) {
        arena->high_water_mark = arena->current_offset;
    }
    
    // Memory is already zeroed from initialization
    return allocation;
}

void arena_reset(ParticleArena_t* arena) {
    if (arena == NULL || arena->memory_block == NULL) {
        return;
    }
    
    // Just reset offset, don't zero memory (will be overwritten)
    arena->current_offset = 0;
    arena->allocation_count = 0;
}

void destroy_particle_arena(ParticleArena_t* arena) {
    if (arena == NULL) {
        return;
    }
    
    if (arena->memory_block != NULL) {
#ifdef _WIN32
        _aligned_free(arena->memory_block);
#else
        free(arena->memory_block);
#endif
        arena->memory_block = NULL;
    }
    
    arena->block_size = 0;
    arena->current_offset = 0;
    arena->allocation_count = 0;
    arena->high_water_mark = 0;
}

/* ============================================================================
 * UTILITY & DEBUG IMPLEMENTATIONS
 * ========================================================================== */

void print_state_vector(const StateVector_t* state, const char* label) {
    if (state == NULL) {
        printf("%s: NULL state vector\n", label ? label : "State");
        return;
    }
    
    if (label) {
        printf("=== %s ===\n", label);
    }
    
    printf("Epoch JD: %.6f\n", state->epoch_jd);
    printf("Object ID: %lu\n", state->object_id);
    printf("Position (km): [%12.6f, %12.6f, %12.6f]\n",
           state->position_km[0], state->position_km[1], state->position_km[2]);
    printf("Velocity (km/s): [%12.6f, %12.6f, %12.6f]\n",
           state->velocity_kms[0], state->velocity_kms[1], state->velocity_kms[2]);
    
    if (state->mass_kg > 0.0) {
        printf("Mass: %.3f kg\n", state->mass_kg);
    }
    
    if (state->cross_section_m2 > 0.0) {
        printf("Cross-section: %.6f m²\n", state->cross_section_m2);
    }
    
    // Compute derived quantities
    double altitude = get_altitude_km(state);
    double energy = compute_specific_energy(state);
    double period = compute_orbital_period(state);
    
    printf("Altitude: %.3f km\n", altitude);
    printf("Specific energy: %.6f km²/s²\n", energy);
    
    if (isfinite(period)) {
        printf("Orbital period: %.1f minutes\n", period / 60.0);
    } else {
        printf("Orbital period: Infinite (escape trajectory)\n");
    }
    
    printf("Is debris: %s\n", state->is_debris ? "Yes" : "No");
    printf("\n");
}

void print_keplerian_elements(const KeplerianElements_t* elements, const char* label) {
    if (elements == NULL) {
        printf("%s: NULL elements\n", label ? label : "Elements");
        return;
    }
    
    if (label) {
        printf("=== %s ===\n", label);
    }
    
    printf("Epoch JD: %.6f\n", elements->epoch_jd);
    printf("Semi-major axis: %.6f km\n", elements->semi_major_axis);
    printf("Eccentricity: %.6f\n", elements->eccentricity);
    printf("Inclination: %.3f°\n", elements->inclination * RAD_TO_DEG);
    printf("RAAN: %.3f°\n", elements->raan * RAD_TO_DEG);
    printf("Argument of perigee: %.3f°\n", elements->arg_perigee * RAD_TO_DEG);
    printf("True anomaly: %.3f°\n", elements->true_anomaly * RAD_TO_DEG);
    
    // Compute derived parameters
    double perigee = elements->semi_major_axis * (1.0 - elements->eccentricity);
    double apogee = elements->semi_major_axis * (1.0 + elements->eccentricity);
    
    printf("Perigee altitude: %.3f km\n", perigee - EARTH_RADIUS_EQUATORIAL);
    printf("Apogee altitude: %.3f km\n", apogee - EARTH_RADIUS_EQUATORIAL);
    
    if (elements->eccentricity < 1.0) {
        double period = 2.0 * M_PI * sqrt(pow(elements->semi_major_axis, 3) / EARTH_MU);
        printf("Orbital period: %.1f minutes\n", period / 60.0);
    }
    
    printf("\n");
}

uint32_t compute_state_checksum(const StateVector_t* state) {
    if (state == NULL) {
        return 0;
    }
    
    // Simple Fletcher-32 checksum
    uint32_t sum1 = 0;
    uint32_t sum2 = 0;
    
    // Hash epoch
    uint8_t* epoch_bytes = (uint8_t*)&state->epoch_jd;
    for (size_t i = 0; i < sizeof(double); i++) {
        sum1 = (sum1 + epoch_bytes[i]) % 65535;
        sum2 = (sum2 + sum1) % 65535;
    }
    
    // Hash position and velocity
    for (int i = 0; i < 3; i++) {
        uint8_t* pos_bytes = (uint8_t*)&state->position_km[i];
        uint8_t* vel_bytes = (uint8_t*)&state->velocity_kms[i];
        
        for (size_t j = 0; j < sizeof(double); j++) {
            sum1 = (sum1 + pos_bytes[j]) % 65535;
            sum2 = (sum2 + sum1) % 65535;
            
            sum1 = (sum1 + vel_bytes[j]) % 65535;
            sum2 = (sum2 + sum1) % 65535;
        }
    }
    
    return (sum2 << 16) | sum1;
}

bool state_vectors_equal(const StateVector_t* a, 
                        const StateVector_t* b,
                        double position_tol_km,
                        double velocity_tol_kms) {
    if (a == NULL || b == NULL) {
        return a == b;  // Both NULL are equal
    }
    
    // Check epoch
    if (fabs(a->epoch_jd - b->epoch_jd) > FLOAT_EPSILON) {
        return false;
    }
    
    // Check object ID
    if (a->object_id != b->object_id) {
        return false;
    }
    
    // Check position
    for (int i = 0; i < 3; i++) {
        if (fabs(a->position_km[i] - b->position_km[i]) > position_tol_km) {
            return false;
        }
    }
    
    // Check velocity
    for (int i = 0; i < 3; i++) {
        if (fabs(a->velocity_kms[i] - b->velocity_kms[i]) > velocity_tol_kms) {
            return false;
        }
    }
    
    // Check other fields (with looser tolerance)
    if (fabs(a->mass_kg - b->mass_kg) > 1e-6) {
        return false;
    }
    
    if (fabs(a->cross_section_m2 - b->cross_section_m2) > 1e-9) {
        return false;
    }
    
    if (a->is_debris != b->is_debris) {
        return false;
    }
    
    return true;
}