/**
 * @file state_vector.h
 * @brief Core state vector operations for orbital propagation
 * 
 * This module implements the fundamental StateVector_t structure and associated
 * operations that form the backbone of the DEBRIS propagation system. All orbital
 * states in the system are represented using this consistent structure.
 * 
 * @section conventions Conventions
 * - Positions in kilometers (km) in ECI (J2000) frame
 * - Velocities in kilometers per second (km/s)
 * - Time in Julian days (J2000 epoch)
 * - Right-handed coordinate system
 * 
 * @section precision Numerical Precision
 * - Double precision (64-bit) throughout
 * - Kahan-Babushka-Neumaier summation for all accumulations
 * - Relative tolerance: 1e-12 for unit tests
 * 
 * @author DEBRIS Engineering Team
 * @date 2025/12/16
 * @version 1.0
 * @copyright MIT License
 */

#ifndef STATE_VECTOR_H
#define STATE_VECTOR_H

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include "core/error_handling.h"

/* ============================================================================
 * CONSTANTS & MACROS
 * ========================================================================== */

#ifndef M_PI
#define M_PI 3.14159265358979323846 
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/** Earth gravitational parameter (μ) in km³/s² */
#define EARTH_MU 398600.4418

/** WGS84 Earth equatorial radius in kilometers */
#define EARTH_RADIUS_EQUATORIAL 6378.137

/** WGS84 Earth flattening factor (1/f) */
#define EARTH_FLATTENING 298.257223563

/** Maximum allowable eccentricity for stable orbit */
#define MAX_STABLE_ECCENTRICITY 0.999

/** Minimum perigee altitude for non-decaying orbit (km) */
#define MIN_NON_DECAYING_ALTITUDE 200.0

/** Conversion: degrees to radians */
#define DEG_TO_RAD 0.017453292519943295

/** Conversion: radians to degrees */
#define RAD_TO_DEG 57.29577951308232

/** Safe floating-point comparison tolerance */
#define FLOAT_EPSILON 1e-12

/** Check if value is approximately zero */
#define IS_ZERO(x) (fabs(x) < FLOAT_EPSILON)

/** Check if two floats are approximately equal */
#define FLOAT_EQ(a, b) (fabs((a) - (b)) < FLOAT_EPSILON)


/* ============================================================================
 * CORE DATA STRUCTURES
 * ========================================================================== */

/**
 * @brief Orbital state vector in Cartesian coordinates
 * 
 * The fundamental data structure representing an object's orbital state
 * at a specific epoch. Used throughout the DEBRIS system for propagation,
 * conjunction analysis, and visualization.
 */
typedef struct {
    double epoch_jd;           ///< Julian date (J2000) of this state
    double position_km[3];     ///< ECI J2000 position (km)
    double velocity_kms[3];    ///< ECI J2000 velocity (km/s)
    double mass_kg;            ///< Mass in kilograms (0.0 if unknown)
    double cross_section_m2;   ///< Cross-sectional area (m²) for drag/SRP
    double reflectance;        ///< Reflectivity coefficient (0.0 to 2.0)
    uint64_t object_id;        ///< Unique identifier for the object
    bool requires_cleanup;     ///< Internal: tracks allocation state
    bool is_debris;            ///< True if object is debris (different force models)
} StateVector_t;

/**
 * @brief Classical Keplerian orbital elements
 * 
 * Alternative representation of orbital state using Keplerian elements.
 * Used for analytical solutions and regime detection.
 */
typedef struct {
    double semi_major_axis;    ///< a: Semi-major axis (km)
    double eccentricity;       ///< e: Eccentricity (0 ≤ e < ∞)
    double inclination;        ///< i: Inclination (radians, 0 to π)
    double raan;               ///< Ω: Right Ascension of Ascending Node (radians)
    double arg_perigee;        ///< ω: Argument of perigee (radians)
    double true_anomaly;       ///< ν: True anomaly (radians)
    double epoch_jd;           ///< Epoch of these elements
} KeplerianElements_t;

/**
 * @brief Modified equinoctial elements for near-circular/polar orbits
 * 
 * Alternative element set that avoids singularities at zero eccentricity
 * and zero inclination. Used for special-case transformations.
 */
typedef struct {
    double p;                  ///< Semi-latus rectum (km)
    double f;                  ///< e*cos(ω+Ω)
    double g;                  ///< e*sin(ω+Ω)
    double h;                  ///< tan(i/2)*cos(Ω)
    double k;                  ///< tan(i/2)*sin(Ω)
    double L;                  ///< True longitude (radians)
    double epoch_jd;           ///< Epoch of these elements
} EquinoctialElements_t;

/**
 * @brief Orbital plane coordinate system (Radial, Tangential, Normal)
 * 
 * Local orbital frame used for relative motion analysis.
 */
typedef struct {
    double radial[3];          ///< R: Unit vector from Earth center to satellite
    double tangential[3];      ///< T: Unit vector in velocity direction (in-plane)
    double normal[3];          ///< N: Unit vector normal to orbital plane
} RTNFrame_t;

/* ============================================================================
 * MEMORY ARENA FOR PARTICLE MANAGEMENT
 * ========================================================================== */

/**
 * @brief Arena allocator for efficient particle memory management
 * 
 * Pre-allocates a large block of memory for debris particles to avoid
 * fragmentation and improve cache locality during cloud propagation.
 */
typedef struct {
    uint8_t* memory_block;     ///< Base pointer to allocated memory
    size_t block_size;         ///< Total size of memory block (bytes)
    size_t current_offset;     ///< Current allocation offset
    size_t allocation_count;   ///< Number of allocations made
    size_t high_water_mark;    ///< Maximum memory used
} ParticleArena_t;

/* ============================================================================
 * LIFE CYCLE FUNCTIONS
 * ========================================================================== */

/**
 * @brief Create a new StateVector on the heap
 * 
 * Allocates and initializes a StateVector_t. The returned pointer
 * must be freed with destroy_state_vector().
 * 
 * @param epoch_jd Julian date for the state vector
 * @param object_id Unique identifier for the object
 * @return StateVector_t* New state vector, or NULL on allocation failure
 * 
 * @warning Always check return value for NULL
 * @post State vector is zero-initialized except for epoch and object_id
 * 
 * @see destroy_state_vector()
 */
StateVector_t* create_state_vector(double epoch_jd, uint64_t object_id);

/**
 * @brief Create a StateVector from raw position/velocity
 * 
 * Convenience function to create a state vector with specified
 * position and velocity. Mass and cross-section are set to defaults.
 * 
 * @param epoch_jd Julian date
 * @param pos_km Position vector (km) in ECI J2000
 * @param vel_kms Velocity vector (km/s) in ECI J2000
 * @param object_id Unique object identifier
 * @return StateVector_t* New state vector, or NULL on failure
 */
StateVector_t* create_state_vector_from_arrays(double epoch_jd, 
                                               const double pos_km[3],
                                               const double vel_kms[3],
                                               uint64_t object_id);

/**
 * @brief Safely destroy a StateVector
 * 
 * Frees all resources associated with the state vector and sets
 * the pointer to NULL to prevent dangling references.
 * 
 * @param state_pointer Pointer to the state vector pointer
 * 
 * @note Safe to call with NULL pointer or pointer to NULL
 * @post *state_pointer == NULL
 * 
 * @see create_state_vector()
 */
void destroy_state_vector(StateVector_t** state_pointer);

/**
 * @brief Create a deep copy of a StateVector
 * 
 * Allocates a new StateVector_t and copies all fields from source.
 * 
 * @param source State vector to copy
 * @return StateVector_t* New copy, or NULL on allocation failure
 */
StateVector_t* copy_state_vector(const StateVector_t* source);

/**
 * @brief Initialize an existing StateVector with default values
 * 
 * Sets all fields to safe defaults. Useful for stack-allocated vectors.
 * 
 * @param state State vector to initialize (must not be NULL)
 * @param epoch_jd Julian date
 * @param object_id Unique object identifier
 * @return Result_t RESULT_OK on success, error code on failure
 */
Result_t initialize_state_vector(StateVector_t* state, 
                                double epoch_jd, 
                                uint64_t object_id);

/* ============================================================================
 * VECTOR OPERATIONS (WITH KBN SUMMATION)
 * ========================================================================== */

/**
 * @brief Compute vector norm with Kahan-Babushka-Neumaier summation
 * 
 * Computes the Euclidean norm of a 3D vector with compensated summation
 * to minimize numerical error.
 * 
 * @param vec Input vector (3 elements)
 * @return double Norm of the vector
 * 
 * @note Uses KBN summation for numerical stability
 * @warning Assumes vec is not NULL and has at least 3 elements
 */
double vector_norm_kbn(const double vec[3]);

/**
 * @brief Compute dot product with KBN summation
 * 
 * Computes the dot product of two vectors with compensated summation.
 * 
 * @param a First vector (3 elements)
 * @param b Second vector (3 elements)
 * @return double Dot product a·b
 * 
 * @note Uses KBN summation for numerical stability
 */
double vector_dot_product_kbn(const double a[3], const double b[3]);

/**
 * @brief Compute cross product
 * 
 * Computes the cross product a × b.
 * 
 * @param a First vector (input)
 * @param b Second vector (input)
 * @param result Cross product (output, must have space for 3 elements)
 * 
 * @warning result must not alias a or b
 */
void vector_cross_product(const double a[3], const double b[3], double result[3]);

/**
 * @brief Add scaled vector: result = a + scale * b
 * 
 * Performs linear combination with scaling. Useful for integrator steps.
 * 
 * @param result Output vector (can alias a or b)
 * @param a First vector
 * @param scale Scaling factor
 * @param b Second vector
 */
void vector_add_scaled(double result[3], 
                      const double a[3], 
                      double scale, 
                      const double b[3]);

/**
 * @brief Normalize a vector in-place
 * 
 * Scales the vector to unit length. If the vector has zero norm,
 * it is set to zero (no crash).
 * 
 * @param vec Vector to normalize (modified in-place)
 * @return double Original norm before normalization
 */
double vector_normalize(double vec[3]);

/**
 * @brief Compute angle between two vectors
 * 
 * Computes the angle between vectors a and b using dot product.
 * Returns angle in radians [0, π].
 * 
 * @param a First vector
 * @param b Second vector
 * @return double Angle in radians
 */
double vector_angle_between(const double a[3], const double b[3]);

/* ============================================================================
 * PHYSICAL PROPERTIES & ORBITAL PARAMETERS
 * ========================================================================== */

/**
 * @brief Compute altitude above Earth's surface (WGS84)
 * 
 * Calculates the distance from the Earth's surface (WGS84 ellipsoid)
 * to the satellite position.
 * 
 * @param state State vector to compute altitude for
 * @return double Altitude in kilometers, or NAN if state is NULL
 * 
 * @note Uses WGS84 Earth radius (6378.137 km)
 * @warning Returns NAN for positions below Earth's surface
 */
double get_altitude_km(const StateVector_t* state);

/**
 * @brief Compute specific orbital energy (energy per unit mass)
 * 
 * ε = v²/2 - μ/r
 * Negative for elliptical orbits, zero for parabolic, positive for hyperbolic.
 * 
 * @param state State vector
 * @return double Specific energy (km²/s²)
 */
double compute_specific_energy(const StateVector_t* state);

/**
 * @brief Compute specific angular momentum vector
 * 
 * h = r × v
 * 
 * @param state State vector
 * @param h_vec Output: angular momentum vector (km²/s)
 * @return Result_t RESULT_OK on success
 */
Result_t compute_angular_momentum_vector(const StateVector_t* state, double h_vec[3]);

/**
 * @brief Compute magnitude of specific angular momentum
 * 
 * h = |r × v|
 * 
 * @param state State vector
 * @return double Angular momentum magnitude (km²/s)
 */
double compute_angular_momentum_magnitude(const StateVector_t* state);

/**
 * @brief Compute orbital period (for elliptical orbits)
 * 
 * T = 2π√(a³/μ)
 * 
 * @param state State vector
 * @return double Period in seconds, or INFINITY for hyperbolic orbits
 */
double compute_orbital_period(const StateVector_t* state);

/**
 * @brief Compute semi-major axis from state vector
 * 
 * a = -μ/(2ε)
 * 
 * @param state State vector
 * @return double Semi-major axis (km), positive for ellipses, negative for hyperbolas
 */
double compute_semi_major_axis(const StateVector_t* state);

/**
 * @brief Compute eccentricity vector
 * 
 * e = (v × h)/μ - r/|r|
 * 
 * @param state State vector
 * @param e_vec Output: eccentricity vector (unitless)
 * @return Result_t RESULT_OK on success
 */
Result_t compute_eccentricity_vector(const StateVector_t* state, double e_vec[3]);

/**
 * @brief Compute eccentricity magnitude
 * 
 * @param state State vector
 * @return double Eccentricity (0 ≤ e < ∞)
 */
double compute_eccentricity(const StateVector_t* state);

/* ============================================================================
 * COORDINATE TRANSFORMATIONS
 * ========================================================================== */

/**
 * @brief Convert Cartesian state to Keplerian elements
 * 
 * Transforms from Cartesian (position, velocity) to Keplerian elements.
 * Handles special cases: circular, equatorial, and parabolic orbits.
 * 
 * @param state Input state vector (must be valid)
 * @param elements Output Keplerian elements
 * @return Result_t RESULT_OK on success, error code on failure
 * 
 * @note Uses modified algorithm to handle singularities
 * @warning May return RESULT_SINGULAR_MATRIX for degenerate cases
 */
Result_t state_vector_to_elements(const StateVector_t* state, 
                                 KeplerianElements_t* elements);

/**
 * @brief Convert Keplerian elements to Cartesian state
 * 
 * Transforms from Keplerian elements to Cartesian (position, velocity).
 * 
 * @param elements Input Keplerian elements
 * @param state Output state vector (allocated by caller)
 * @return Result_t RESULT_OK on success, error code on failure
 */
Result_t elements_to_state_vector(const KeplerianElements_t* elements,
                                 StateVector_t* state);

/**
 * @brief Convert to equinoctial elements (non-singular)
 * 
 * Transforms to equinoctial elements, which avoid singularities
 * at zero eccentricity and zero inclination.
 * 
 * @param state Input state vector
 * @param equinoctial Output equinoctial elements
 * @return Result_t RESULT_OK on success
 */
Result_t state_vector_to_equinoctial(const StateVector_t* state,
                                    EquinoctialElements_t* equinoctial);

/**
 * @brief Compute RTN (Radial, Tangential, Normal) frame
 * 
 * Computes the local orbital coordinate system.
 * 
 * @param state State vector
 * @param rtn Output RTN frame
 * @return Result_t RESULT_OK on success
 */
Result_t compute_rtn_frame(const StateVector_t* state, RTNFrame_t* rtn);

/**
 * @brief Transform vector from ECI to RTN frame
 * 
 * @param state Reference state vector defining RTN frame
 * @param vec_eci Vector in ECI frame
 * @param vec_rtn Output vector in RTN frame
 * @return Result_t RESULT_OK on success
 */
Result_t transform_eci_to_rtn(const StateVector_t* state,
                             const double vec_eci[3],
                             double vec_rtn[3]);

/**
 * @brief Transform vector from RTN to ECI frame
 * 
 * @param state Reference state vector defining RTN frame
 * @param vec_rtn Vector in RTN frame
 * @param vec_eci Output vector in ECI frame
 * @return Result_t RESULT_OK on success
 */
Result_t transform_rtn_to_eci(const StateVector_t* state,
                             const double vec_rtn[3],
                             double vec_eci[3]);

/* ============================================================================
 * VALIDATION & SAFETY CHECKS
 * ========================================================================== */

/**
 * @brief Validate state vector for physical plausibility
 * 
 * Checks for NaN/Inf, valid position/velocity magnitudes,
 * and basic orbital sanity.
 * 
 * @param state State vector to validate
 * @return Result_t RESULT_OK if valid, error code otherwise
 */
Result_t validate_state_vector(const StateVector_t* state);

/**
 * @brief Check for NaN or Inf values in state vector
 * 
 * @param state State vector to check
 * @return true If any field is NaN or Inf
 * @return false If all fields are finite
 */
bool check_for_nan_inf(const StateVector_t* state);

/**
 * @brief Check if orbit is stable (elliptical and above surface)
 * 
 * An orbit is considered stable if:
 * - Eccentricity < 1.0 (elliptical)
 * - Perigee altitude > 0.0 km
 * - Not on escape trajectory
 * 
 * @param state State vector
 * @return true Orbit is stable
 * @return false Orbit is unstable or decaying
 */
bool is_orbit_stable(const StateVector_t* state);

/**
 * @brief Check if orbit is decaying (perigee < threshold)
 * 
 * @param state State vector
 * @param threshold_km Altitude threshold in km (default: 200.0)
 * @return true Orbit is decaying (will re-enter)
 * @return false Orbit is not decaying
 */
bool is_decaying_orbit(const StateVector_t* state, double threshold_km);

/**
 * @brief Check if object has reached escape velocity
 * 
 * @param state State vector
 * @return true Object is on escape trajectory
 * @return false Object is bound to Earth
 */
bool is_escape_trajectory(const StateVector_t* state);

/**
 * @brief Check if orbit is circular (within tolerance)
 * 
 * @param state State vector
 * @param tolerance Maximum eccentricity for "circular"
 * @return true Orbit is circular (e < tolerance)
 * @return false Orbit is not circular
 */
bool is_circular_orbit(const StateVector_t* state, double tolerance);

/* ============================================================================
 * MEMORY ARENA FUNCTIONS
 * ========================================================================== */

/**
 * @brief Initialize a particle memory arena
 * 
 * Allocates a large contiguous block of memory for particle storage.
 * 
 * @param arena Arena to initialize
 * @param size_mb Size of arena in megabytes
 * @return Result_t RESULT_OK on success
 */
Result_t initialize_particle_arena(ParticleArena_t* arena, size_t size_mb);

/**
 * @brief Allocate particles from arena
 * 
 * @param arena Memory arena
 * @param particle_count Number of particles to allocate
 * @param particle_size Size of each particle in bytes
 * @return void* Pointer to allocated particles, or NULL if arena full
 */
void* arena_allocate_particles(ParticleArena_t* arena, 
                              size_t particle_count, 
                              size_t particle_size);

/**
 * @brief Reset arena (free all allocations without freeing memory block)
 * 
 * @param arena Arena to reset
 */
void arena_reset(ParticleArena_t* arena);

/**
 * @brief Destroy arena and free all memory
 * 
 * @param arena Arena to destroy
 */
void destroy_particle_arena(ParticleArena_t* arena);

/* ============================================================================
 * UTILITY & DEBUG FUNCTIONS
 * ========================================================================== */

/**
 * @brief Print state vector to stdout for debugging
 * 
 * @param state State vector to print
 * @param label Optional label for output
 */
void print_state_vector(const StateVector_t* state, const char* label);

/**
 * @brief Print Keplerian elements to stdout
 * 
 * @param elements Elements to print
 * @param label Optional label for output
 */
void print_keplerian_elements(const KeplerianElements_t* elements, const char* label);

/**
 * @brief Compute checksum of state vector for validation
 * 
 * @param state State vector
 * @return uint32_t 32-bit checksum
 */
uint32_t compute_state_checksum(const StateVector_t* state);

/**
 * @brief Compare two state vectors for equality within tolerance
 * 
 * @param a First state vector
 * @param b Second state vector
 * @param position_tol_km Position tolerance (km)
 * @param velocity_tol_kms Velocity tolerance (km/s)
 * @return true Vectors are equal within tolerance
 * @return false Vectors differ
 */
bool state_vectors_equal(const StateVector_t* a, 
                        const StateVector_t* b,
                        double position_tol_km,
                        double velocity_tol_kms);

#endif /* STATE_VECTOR_H */