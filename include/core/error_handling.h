/**
 * @file error_handling.h
 * @brief Centralized error handling and logging for DEBRIS system
 * 
 * This module provides:
 * 1. Standardized error codes (Result_t) used by all modules
 * 2. Structured logging with context (file, line, function)
 * 3. Assertion and contract checking macros
 * 4. Error propagation utilities
 * 5. Graceful degradation protocols
 * 
 * @section design Design Philosophy
 * - Never crash silently: All errors are logged with context
 * - Graceful degradation: Continue with reduced accuracy when possible
 * - Deterministic behavior: Errors are predictable and reproducible
 * - Thread-safe: Each thread has independent error context
 * 
 * @author DEBRIS Engineering Team
 * @date 2024
 * @version 1.0
 * @copyright MIT License
 */

#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

/* ============================================================================
 * COMPILE-TIME CONFIGURATION
 * ========================================================================== */

/** Enable/disable assertions in debug builds */
#ifndef DEBRIS_DEBUG
#define DEBRIS_DEBUG 1
#endif

/** Log level: 0=OFF, 1=ERROR, 2=WARNING, 3=INFO, 4=DEBUG, 5=TRACE */
#ifndef DEBRIS_LOG_LEVEL
#define DEBRIS_LOG_LEVEL 3  // Default: INFO
#endif

/** Maximum error message length */
#define MAX_ERROR_MESSAGE 1024

/** Maximum nested error context depth */
#define MAX_ERROR_DEPTH 16

/* ============================================================================
 * ERROR CODES - EXHAUSTIVE LIST
 * ========================================================================== */

/**
 * @brief Comprehensive error code enumeration for the DEBRIS system
 * 
 * Every function in DEBRIS returns a Result_t indicating success or
 * specific failure mode. Error codes are organized by category.
 */
typedef enum {
    /* ========== SUCCESS (0) ========== */
    RESULT_OK = 0,
    RESULT_OK_WITH_WARNING = 1,  ///< Operation succeeded but with warnings
    
    /* ========== INPUT ERRORS (1000-1099) ========== */
    RESULT_NULL_POINTER = 1000,       ///< Unexpected NULL pointer
    RESULT_INVALID_ARGUMENT = 1001,   ///< Argument violates precondition
    RESULT_OUT_OF_RANGE = 1002,       ///< Value outside valid range
    RESULT_INVALID_FORMAT = 1003,     ///< Input format invalid (TLE, etc.)
    RESULT_FILE_NOT_FOUND = 1004,     ///< Requested file doesn't exist
    RESULT_PERMISSION_DENIED = 1005,  ///< Insufficient permissions
    RESULT_END_OF_FILE = 1006,        ///< Unexpected end of file
    RESULT_INVALID_CHECKSUM = 1007,   ///< Checksum validation failed
    
    /* ========== NUMERICAL ERRORS (1100-1199) ========== */
    RESULT_NUMERICAL_ERROR = 1100,    ///< General numerical error
    RESULT_DIVERGENCE = 1101,         ///< Numerical method diverged
    RESULT_SINGULAR_MATRIX = 1102,    ///< Matrix is singular/non-invertible
    RESULT_MAX_ITERATIONS = 1103,     ///< Max iterations exceeded
    RESULT_TOLERANCE_NOT_MET = 1104,  ///< Convergence tolerance not met
    RESULT_OVERFLOW = 1105,           ///< Numerical overflow
    RESULT_UNDERFLOW = 1106,          ///< Numerical underflow
    RESULT_NAN_DETECTED = 1107,       ///< NaN value encountered
    RESULT_INF_DETECTED = 1108,       ///< Infinite value encountered
    
    /* ========== RESOURCE ERRORS (1200-1299) ========== */
    RESULT_OUT_OF_MEMORY = 1200,      ///< Memory allocation failed
    RESULT_OUT_OF_GPU_MEMORY = 1201,  ///< GPU memory exhausted
    RESULT_DISK_FULL = 1202,          ///< No space left on device
    RESULT_TOO_MANY_FILES = 1203,     ///< Too many open files
    RESULT_RESOURCE_BUSY = 1204,      ///< Resource temporarily unavailable
    RESULT_NETWORK_ERROR = 1205,      ///< Network operation failed
    RESULT_TIMEOUT = 1206,            ///< Operation timed out
    
    /* ========== PHYSICS/DOMAIN ERRORS (1300-1399) ========== */
    RESULT_INVALID_ORBIT = 1300,      ///< Orbit parameters invalid
    RESULT_ORBIT_DECAYED = 1301,      ///< Orbit has decayed/re-entered
    RESULT_ESCAPE_TRAJECTORY = 1302,  ///< Object on escape trajectory
    RESULT_SUBORBITAL = 1303,         ///< Object is suborbital
    RESULT_COLLISION_DETECTED = 1304, ///< Collision occurred in propagation
    RESULT_INVALID_MANEUVER = 1305,   ///< Maneuver parameters invalid
    RESULT_OUT_OF_FOV = 1306,         ///< Object outside sensor field of view
    RESULT_BELOW_HORIZON = 1307,      ///< Object below local horizon
    
    /* ========== SYSTEM/STATE ERRORS (1400-1499) ========== */
    RESULT_NOT_INITIALIZED = 1400,    ///< System/component not initialized
    RESULT_ALREADY_INITIALIZED = 1401,///< Attempt to re-initialize
    RESULT_INVALID_STATE = 1402,      ///< System in invalid state for operation
    RESULT_SHUTTING_DOWN = 1403,      ///< System is shutting down
    RESULT_CONFIG_ERROR = 1404,       ///< Configuration error
    
    /* ========== FEATURE ERRORS (1500-1599) ========== */
    RESULT_NOT_IMPLEMENTED = 1500,    ///< Feature not implemented
    RESULT_DEPRECATED = 1501,         ///< Feature deprecated
    RESULT_UNSUPPORTED = 1502,        ///< Feature not supported on this platform
    
    /* ========== INTEGRATOR ERRORS (1600-1699) ========== */
    RESULT_INTEGRATOR_STUCK = 1600,   ///< Integrator cannot make progress
    RESULT_STIFF_SYSTEM = 1601,       ///< System too stiff for current method
    RESULT_SMALL_STEPSIZE = 1602,     ///< Step size reduced below minimum
    RESULT_SYMPLECTIC_BREAKDOWN = 1603, ///< Symplectic integrator invalid
    
    /* ========== FORCE MODEL ERRORS (1700-1799) ========== */
    RESULT_MODEL_OUT_OF_RANGE = 1700, ///< Model used outside valid domain
    RESULT_ATMOSPHERE_UNAVAILABLE = 1701, ///< Atmospheric data unavailable
    RESULT_EPHEMERIS_UNAVAILABLE = 1702, ///< Ephemeris data unavailable
    
    /* ========== CLOUD/GMM ERRORS (1800-1899) ========== */
    RESULT_GMM_NOT_CONVERGED = 1800,  ///< GMM fitting didn't converge
    RESULT_INVALID_CLOUD = 1801,      ///< Cloud representation invalid
    RESULT_TOO_MANY_COMPONENTS = 1802, ///< Too many GMM components
    
    /* ========== CONJUNCTION ERRORS (1900-1999) ========== */
    RESULT_NO_CONJUNCTION = 1900,     ///< No conjunction found (not an error)
    RESULT_AMBIGUOUS_TCA = 1901,      ///< Multiple TCA candidates
    RESULT_COVARIANCE_INVALID = 1902, ///< Covariance matrix invalid
    
    /* ========== INTERNAL ERRORS (2000-2099) ========== */
    RESULT_INTERNAL_ERROR = 2000,     ///< Unexpected internal error
    RESULT_ASSERTION_FAILED = 2001,   ///< Assertion failed (debug builds)
    RESULT_INVARIANT_VIOLATED = 2002, ///< Invariant violated
    RESULT_UNREACHABLE_CODE = 2003,   ///< Code marked unreachable was reached
    
    /* ========== RESILIENCE/GRACEFUL DEGRADATION (2100-2199) ========== */
    RESULT_DEGRADED_ACCURACY = 2100,  ///< Operation completed with reduced accuracy
    RESULT_USING_FALLBACK = 2101,     ///< Using fallback method
    RESULT_CACHE_MISS = 2102,         ///< Cache miss (normal operation)
    RESULT_RECOVERED_ERROR = 2103,    ///< Error was recovered from
    
    /* ========== VALIDATION ERRORS (2200-2299) ========== */
    RESULT_VALIDATION_FAILED = 2200,  ///< Validation against reference failed
    RESULT_TEST_FAILED = 2201,        ///< Unit/integration test failed
    RESULT_BENCHMARK_FAILED = 2202,   ///< Performance benchmark failed
    
    /* ========== THREADING/CONCURRENCY (2300-2399) ========== */
    RESULT_THREAD_ERROR = 2300,       ///< Thread operation failed
    RESULT_DEADLOCK_DETECTED = 2301,  ///< Potential deadlock detected
    RESULT_RACE_CONDITION = 2302,     ///< Race condition detected
    
    /* ========== I/O AND SERIALIZATION (2400-2499) ========== */
    RESULT_SERIALIZATION_ERROR = 2400, ///< Serialization failed
    RESULT_DESERIALIZATION_ERROR = 2401, ///< Deserialization failed
    RESULT_VERSION_MISMATCH = 2402,   ///< Version mismatch in serialized data
    
    /* ========== UNKNOWN ERROR (MAX) ========== */
    RESULT_UNKNOWN_ERROR = 9999       ///< Unknown/unclassified error
} Result_t;

/**
 * @brief Error severity levels for logging and handling
 */
typedef enum {
    SEVERITY_DEBUG = 0,    ///< Debug information (verbose)
    SEVERITY_INFO = 1,     ///< Informational messages
    SEVERITY_WARNING = 2,  ///< Warnings (non-fatal issues)
    SEVERITY_ERROR = 3,    ///< Errors (operation failed)
    SEVERITY_FATAL = 4     ///< Fatal errors (cannot continue)
} ErrorSeverity_t;

/**
 * @brief Error context structure for nested error tracking
 * 
 * Maintains a stack of error contexts to preserve call chain information.
 */
typedef struct {
    Result_t code;                  ///< Error code
    const char* file;              ///< Source file where error occurred
    int line;                      ///< Line number
    const char* function;          ///< Function name
    const char* message;           ///< Human-readable message
    uint64_t timestamp_ns;         ///< Timestamp in nanoseconds
    uint32_t thread_id;            ///< Thread ID for thread-safe logging
    struct ErrorContext* previous; ///< Previous error in chain (or NULL)
} ErrorContext;

/* ============================================================================
 * LOGGING MACROS & FUNCTIONS
 * ========================================================================== */

/**
 * @brief Log a message with severity and context
 * 
 * @param severity Error severity level
 * @param file Source file name (use __FILE__)
 * @param line Line number (use __LINE__)
 * @param function Function name (use __func__)
 * @param format printf-style format string
 * @param ... Format arguments
 */
void debris_log(ErrorSeverity_t severity,
                const char* file,
                int line,
                const char* function,
                const char* format,
                ...);

/**
 * @brief Log an error with full context and error code
 * 
 * This is the primary logging function for errors. It captures:
 * - Error code
 * - Source location
 * - Thread ID
 * - Timestamp
 * - User message
 * 
 * @param code Error code
 * @param file Source file name
 * @param line Line number
 * @param function Function name
 * @param format printf-style format string
 * @param ... Format arguments
 */
void debris_log_error_context(Result_t code,
                              const char* file,
                              int line,
                              const char* function,
                              const char* format,
                              ...);

/**
 * @brief Set the log output file
 * 
 * @param file FILE* to write logs to (NULL for stderr)
 * @return Result_t RESULT_OK on success
 */
Result_t debris_log_set_output(FILE* file);

/**
 * @brief Enable or disable structured logging (JSON format)
 * 
 * @param enabled True for JSON, false for human-readable
 */
void debris_log_set_structured(bool enabled);

/**
 * @brief Set minimum log level
 * 
 * @param level Minimum severity to log
 */
void debris_log_set_level(ErrorSeverity_t level);

/* ============================================================================
 * LOGGING MACROS (Compile-time level filtering)
 * ========================================================================== */

#if DEBRIS_LOG_LEVEL >= 4  // DEBUG
#define LOG_DEBUG(format, ...) \
    debris_log(SEVERITY_DEBUG, __FILE__, __LINE__, __func__, format, ##__VA_ARGS__)
#else
#define LOG_DEBUG(format, ...) ((void)0)
#endif

#if DEBRIS_LOG_LEVEL >= 3  // INFO
#define LOG_INFO(format, ...) \
    debris_log(SEVERITY_INFO, __FILE__, __LINE__, __func__, format, ##__VA_ARGS__)
#else
#define LOG_INFO(format, ...) ((void)0)
#endif

#if DEBRIS_LOG_LEVEL >= 2  // WARNING
#define LOG_WARNING(format, ...) \
    debris_log(SEVERITY_WARNING, __FILE__, __LINE__, __func__, format, ##__VA_ARGS__)
#else
#define LOG_WARNING(format, ...) ((void)0)
#endif

#if DEBRIS_LOG_LEVEL >= 1  // ERROR
#define LOG_ERROR(format, ...) \
    debris_log(SEVERITY_ERROR, __FILE__, __LINE__, __func__, format, ##__VA_ARGS__)

/** @brief Log error with context and return error code */
#define LOG_ERROR_CONTEXT(code, format, ...) \
    do { \
        debris_log_error_context(code, __FILE__, __LINE__, __func__, format, ##__VA_ARGS__); \
        return code; \
    } while(0)

/** @brief Log error with context but don't return (for void functions) */
#define LOG_ERROR_CONTEXT_VOID(code, format, ...) \
    debris_log_error_context(code, __FILE__, __LINE__, __func__, format, ##__VA_ARGS__)

#else
#define LOG_ERROR(format, ...) ((void)0)
#define LOG_ERROR_CONTEXT(code, format, ...) return code
#define LOG_ERROR_CONTEXT_VOID(code, format, ...) ((void)0)
#endif

#if DEBRIS_LOG_LEVEL >= 0  // Always log FATAL
#define LOG_FATAL(format, ...) \
    do { \
        debris_log(SEVERITY_FATAL, __FILE__, __LINE__, __func__, format, ##__VA_ARGS__); \
        debris_handle_fatal_error(); \
    } while(0)
#else
#define LOG_FATAL(format, ...) debris_handle_fatal_error()
#endif

/* ============================================================================
 * ASSERTION & CONTRACT MACROS
 * ========================================================================== */

#if DEBRIS_DEBUG

/** @brief Assert condition is true, abort with message if false */
#define DEBRIS_ASSERT(condition) \
    do { \
        if (!(condition)) { \
            LOG_FATAL("Assertion failed: %s", #condition); \
        } \
    } while(0)

/** @brief Assert condition is true, return error code if false */
#define DEBRIS_REQUIRE(condition, error_code) \
    do { \
        if (!(condition)) { \
            LOG_ERROR_CONTEXT(error_code, "Requirement failed: %s", #condition); \
        } \
    } while(0)

/** @brief Assert condition is true for void functions */
#define DEBRIS_REQUIRE_VOID(condition, error_code) \
    do { \
        if (!(condition)) { \
            LOG_ERROR_CONTEXT_VOID(error_code, "Requirement failed: %s", #condition); \
            return; \
        } \
    } while(0)

/** @brief Check post-condition, log error if false */
#define DEBRIS_ENSURE(condition) \
    do { \
        if (!(condition)) { \
            LOG_ERROR("Post-condition failed: %s", #condition); \
        } \
    } while(0)

/** @brief Check loop invariant */
#define DEBRIS_INVARIANT(condition) \
    do { \
        if (!(condition)) { \
            LOG_ERROR("Invariant violated: %s", #condition); \
        } \
    } while(0)

/** @brief Mark code as unreachable (aborts if reached) */
#define DEBRIS_UNREACHABLE() \
    LOG_FATAL("Reached unreachable code at %s:%d", __FILE__, __LINE__)

#else
// Release builds: assertions become no-ops or simple returns
#define DEBRIS_ASSERT(condition) ((void)0)
#define DEBRIS_REQUIRE(condition, error_code) \
    do { if (!(condition)) return error_code; } while(0)
#define DEBRIS_REQUIRE_VOID(condition, error_code) \
    do { if (!(condition)) return; } while(0)
#define DEBRIS_ENSURE(condition) ((void)0)
#define DEBRIS_INVARIANT(condition) ((void)0)
#define DEBRIS_UNREACHABLE() return RESULT_UNREACHABLE_CODE
#endif

/* ============================================================================
 * ERROR PROPAGATION MACROS
 * ========================================================================== */

/**
 * @brief Propagate error: if err != RESULT_OK, return it immediately
 * 
 * Use for early return on error while preserving the error code.
 */
#define PROPAGATE_ERROR(err) \
    do { \
        Result_t _err = (err); \
        if (_err != RESULT_OK && _err != RESULT_OK_WITH_WARNING) { \
            return _err; \
        } \
    } while(0)

/**
 * @brief Try operation: execute expr, propagate error if not RESULT_OK
 * 
 * Wraps a function call with automatic error checking.
 */
#define TRY(expr) \
    do { \
        Result_t _result = (expr); \
        if (_result != RESULT_OK && _result != RESULT_OK_WITH_WARNING) { \
            return _result; \
        } \
    } while(0)

/**
 * @brief Try operation with context: log error context before propagating
 */
#define TRY_CONTEXT(expr, context_msg) \
    do { \
        Result_t _result = (expr); \
        if (_result != RESULT_OK && _result != RESULT_OK_WITH_WARNING) { \
            LOG_ERROR_CONTEXT(_result, "Failed: %s", context_msg); \
        } \
    } while(0)

/* ============================================================================
 * ERROR HANDLING FUNCTIONS
 * ========================================================================== */

/**
 * @brief Convert error code to human-readable string
 * 
 * @param code Error code
 * @return const char* String representation
 */
const char* result_to_string(Result_t code);

/**
 * @brief Get category of error code
 * 
 * @param code Error code
 * @return const char* Category name
 */
const char* result_get_category(Result_t code);

/**
 * @brief Check if error code indicates success (with or without warnings)
 * 
 * @param code Error code to check
 * @return true Code indicates success
 * @return false Code indicates failure
 */
bool result_is_success(Result_t code);

/**
 * @brief Check if error code indicates a fatal error
 * 
 * @param code Error code to check
 * @return true Code indicates fatal error
 * @return false Code is non-fatal
 */
bool result_is_fatal(Result_t code);

/**
 * @brief Check if error is recoverable (can continue with degradation)
 * 
 * @param code Error code to check
 * @return true Error is recoverable
 * @return false Error is not recoverable
 */
bool result_is_recoverable(Result_t code);

/**
 * @brief Get the last error context for current thread
 * 
 * @return const ErrorContext* Last error (or NULL if no error)
 */
const ErrorContext* debris_get_last_error(void);

/**
 * @brief Clear error context for current thread
 */
void debris_clear_error_context(void);

/**
 * @brief Set custom error handler callback
 * 
 * @param handler Function to call when errors occur (can be NULL)
 */
void debris_set_error_handler(void (*handler)(Result_t, const ErrorContext*));

/**
 * @brief Handle fatal error (called by LOG_FATAL)
 * 
 * Attempts graceful shutdown, then aborts.
 */
void debris_handle_fatal_error(void) __attribute__((noreturn));

/**
 * @brief Initialize error handling subsystem
 * 
 * Must be called once at program startup.
 * 
 * @return Result_t RESULT_OK on success
 */
Result_t debris_error_init(void);

/**
 * @brief Shutdown error handling subsystem
 * 
 * Flushes logs and cleans up resources.
 */
void debris_error_shutdown(void);

/* ============================================================================
 * GRACEFUL DEGRADATION UTILITIES
 * ========================================================================== */

/**
 * @brief Degradation level for operations
 */
typedef enum {
    DEGRADATION_NONE = 0,      ///< Full accuracy
    DEGRADATION_MINIMAL = 1,   ///< Minor accuracy reduction
    DEGRADATION_MODERATE = 2,  ///< Moderate accuracy reduction
    DEGRADATION_SEVERE = 3,    ///< Severe accuracy reduction
    DEGRADATION_CRITICAL = 4   ///< Critical - only basic functionality
} DegradationLevel_t;

/**
 * @brief Get current degradation level for a subsystem
 * 
 * @param subsystem Subsystem identifier
 * @return DegradationLevel_t Current degradation level
 */
DegradationLevel_t debris_get_degradation_level(const char* subsystem);

/**
 * @brief Set degradation level for a subsystem
 * 
 * @param subsystem Subsystem identifier
 * @param level New degradation level
 */
void debris_set_degradation_level(const char* subsystem, DegradationLevel_t level);

/**
 * @brief Check if operation should proceed given current degradation
 * 
 * @param subsystem Subsystem identifier
 * @param required_level Minimum required degradation level
 * @return true Operation can proceed
 * @return false Operation should be skipped
 */
bool debris_can_proceed(const char* subsystem, DegradationLevel_t required_level);

/* ============================================================================
 * MEMORY ERROR CHECKING
 * ========================================================================== */

/**
 * @brief Check for memory errors in allocated block
 * 
 * @param ptr Pointer to check
 * @param size Size of allocation
 * @return Result_t RESULT_OK if valid, error code otherwise
 */
Result_t debris_check_memory(const void* ptr, size_t size);

/**
 * @brief Enable/disable memory error checking
 * 
 * @param enabled True to enable checking
 */
void debris_memory_checking_enable(bool enabled);

#endif /* ERROR_HANDLING_H */