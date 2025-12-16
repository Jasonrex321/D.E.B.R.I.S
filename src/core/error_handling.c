/**
 * @file error_handling.c
 * @brief Implementation of error handling and logging system
 * 
 * @see error_handling.h for documentation
 */

#include "core/error_handling.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <errno.h>

#ifdef _WIN32
#include <windows.h>
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#include <pthread.h>
#include <sys/time.h>  // Added for gettimeofday()
#endif

/* ============================================================================
 * PRIVATE DATA STRUCTURES
 * ========================================================================== */

/** Thread-local error context */
#ifdef _WIN32
static __declspec(thread) ErrorContext* thread_error_context = NULL;
#else
static __thread ErrorContext* thread_error_context = NULL;
#endif

/** Logging configuration */
static struct {
    FILE* output_file;            ///< Log output stream (NULL = stderr)
    bool structured;              ///< True for JSON output
    ErrorSeverity_t min_level;    ///< Minimum severity to log
    bool initialized;             ///< True if subsystem initialized
    void (*error_handler)(Result_t, const ErrorContext*); ///< Custom handler
    bool memory_checking_enabled; ///< Enable memory error checking
} log_config = {
    .output_file = NULL,
    .structured = false,
    .min_level = SEVERITY_INFO,
    .initialized = false,
    .error_handler = NULL,
    .memory_checking_enabled = false
};

/** Degradation levels for subsystems */
static struct {
    const char* subsystem;
    DegradationLevel_t level;
} degradation_levels[] = {
    {"integrator", DEGRADATION_NONE},
    {"force_models", DEGRADATION_NONE},
    {"cloud_propagation", DEGRADATION_NONE},
    {"conjunction", DEGRADATION_NONE},
    {"spectral", DEGRADATION_NONE},
    {NULL, DEGRADATION_NONE}
};

/* ============================================================================
 * PRIVATE HELPER FUNCTIONS
 * ========================================================================== */

/**
 * @brief Get current timestamp in nanoseconds
 */
static uint64_t get_timestamp_ns(void) {
#ifdef _WIN32
    FILETIME ft;
    GetSystemTimeAsFileTime(&ft);
    ULARGE_INTEGER uli;
    uli.LowPart = ft.dwLowDateTime;
    uli.HighPart = ft.dwHighDateTime;
    // Convert from 100-nanosecond intervals since 1601 to nanoseconds
    return (uli.QuadPart - 116444736000000000ULL) * 100;
#else
    // Use clock_gettime if available (POSIX)
    #ifdef _POSIX_TIMERS
        #ifdef CLOCK_REALTIME
            struct timespec ts;
            if (clock_gettime(CLOCK_REALTIME, &ts) == 0) {
                return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
            }
        #endif
    #endif
    
    // Fallback to gettimeofday for older systems
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (uint64_t)tv.tv_sec * 1000000000ULL + (uint64_t)tv.tv_usec * 1000ULL;
#endif
}
/**
 * @brief Get current thread ID
 */
static uint32_t get_thread_id(void) {
#ifdef _WIN32
    return (uint32_t)GetCurrentThreadId();
#else
    return (uint32_t)pthread_self();
#endif
}

/**
 * @brief Get process ID
 */
static uint32_t get_process_id(void) {
#ifdef _WIN32
    return (uint32_t)_getpid();
#else
    return (uint32_t)getpid();
#endif
}

/**
 * @brief Format timestamp as ISO 8601 string
 */
static void format_timestamp(char* buffer, size_t size, uint64_t timestamp_ns) {
    time_t seconds = timestamp_ns / 1000000000ULL;
    int nanoseconds = timestamp_ns % 1000000000ULL;
    
    struct tm tm_time;
#ifdef _WIN32
    gmtime_s(&tm_time, &seconds);
#else
    gmtime_r(&seconds, &tm_time);
#endif
    
    snprintf(buffer, size, "%04d-%02d-%02dT%02d:%02d:%02d.%09dZ",
             tm_time.tm_year + 1900, tm_time.tm_mon + 1, tm_time.tm_mday,
             tm_time.tm_hour, tm_time.tm_min, tm_time.tm_sec, nanoseconds);
}

/**
 * @brief Create new error context
 */
static ErrorContext* create_error_context(Result_t code,
                                         const char* file,
                                         int line,
                                         const char* function,
                                         const char* message) {
    ErrorContext* context = (ErrorContext*)malloc(sizeof(ErrorContext));
    if (!context) {
        return NULL;
    }
    
    context->code = code;
    context->file = file;
    context->line = line;
    context->function = function;
    context->timestamp_ns = get_timestamp_ns();
    context->thread_id = get_thread_id();
    context->previous = thread_error_context;
    
    // Copy message (truncate if necessary)
    if (message) {
        size_t len = strlen(message);
        if (len >= MAX_ERROR_MESSAGE) {
            len = MAX_ERROR_MESSAGE - 1;
        }
        char* msg_copy = (char*)malloc(len + 1);
        if (msg_copy) {
            strncpy(msg_copy, message, len);
            msg_copy[len] = '\0';
            context->message = msg_copy;
        } else {
            context->message = "Memory allocation failed for error message";
        }
    } else {
        context->message = result_to_string(code);
    }
    
    return context;
}

/**
 * @brief Free error context chain
 */
static void free_error_context(ErrorContext* context) {
    while (context) {
        ErrorContext* next = context->previous;
        if (context->message && context->message != result_to_string(context->code)) {
            free((void*)context->message);
        }
        free(context);
        context = next;
    }
}

/**
 * @brief Update thread error context
 */
static void update_error_context(ErrorContext* context) {
    ErrorContext* old = thread_error_context;
    thread_error_context = context;
    if (old && old != context) {
        // Don't free old if it's part of the new chain
        ErrorContext* check = context;
        while (check) {
            if (check == old) {
                return;
            }
            check = check->previous;
        }
        free_error_context(old);
    }
}

/* ============================================================================
 * PUBLIC FUNCTION IMPLEMENTATIONS
 * ========================================================================== */

void debris_log(ErrorSeverity_t severity,
                const char* file,
                int line,
                const char* function,
                const char* format,
                ...) {
    // Check log level
    if (severity < log_config.min_level || !log_config.initialized) {
        return;
    }
    
    // Get output stream
    FILE* output = log_config.output_file ? log_config.output_file : stderr;
    
    // Format message
    char message[MAX_ERROR_MESSAGE];
    va_list args;
    va_start(args, format);
    vsnprintf(message, sizeof(message), format, args);
    va_end(args);
    
    if (log_config.structured) {
        // JSON format
        char timestamp[64];
        format_timestamp(timestamp, sizeof(timestamp), get_timestamp_ns());
        
        const char* severity_str;
        switch (severity) {
            case SEVERITY_DEBUG: severity_str = "DEBUG"; break;
            case SEVERITY_INFO: severity_str = "INFO"; break;
            case SEVERITY_WARNING: severity_str = "WARNING"; break;
            case SEVERITY_ERROR: severity_str = "ERROR"; break;
            case SEVERITY_FATAL: severity_str = "FATAL"; break;
            default: severity_str = "UNKNOWN";
        }
        
        fprintf(output,
                "{\"timestamp\":\"%s\",\"pid\":%u,\"tid\":%u,\"severity\":\"%s\","
                "\"file\":\"%s\",\"line\":%d,\"function\":\"%s\",\"message\":\"%s\"}\n",
                timestamp, get_process_id(), get_thread_id(), severity_str,
                file, line, function, message);
    } else {
        // Human-readable format
        const char* severity_prefix;
        switch (severity) {
            case SEVERITY_DEBUG: severity_prefix = "[DEBUG]"; break;
            case SEVERITY_INFO: severity_prefix = "[INFO]"; break;
            case SEVERITY_WARNING: severity_prefix = "[WARNING]"; break;
            case SEVERITY_ERROR: severity_prefix = "[ERROR]"; break;
            case SEVERITY_FATAL: severity_prefix = "[FATAL]"; break;
            default: severity_prefix = "[UNKNOWN]";
        }
        
        char timestamp[32];
        time_t now = time(NULL);
        struct tm* tm_info = localtime(&now);
        strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", tm_info);
        
        fprintf(output, "%s %s %s:%d (%s): %s\n",
                timestamp, severity_prefix, file, line, function, message);
    }
    
    fflush(output);
}

void debris_log_error_context(Result_t code,
                              const char* file,
                              int line,
                              const char* function,
                              const char* format,
                              ...) {
    // Format message
    char message[MAX_ERROR_MESSAGE];
    va_list args;
    va_start(args, format);
    vsnprintf(message, sizeof(message), format, args);
    va_end(args);
    
    // Create error context
    ErrorContext* context = create_error_context(code, file, line, function, message);
    if (!context) {
        // Fallback to simple logging if context creation fails
        LOG_ERROR("Failed to create error context for code %d: %s", code, message);
        return;
    }
    
    // Update thread context
    update_error_context(context);
    
    // Log the error
    const char* code_str = result_to_string(code);
    LOG_ERROR("Error %d (%s): %s", code, code_str, message);
    
    // Call custom error handler if set
    if (log_config.error_handler) {
        log_config.error_handler(code, context);
    }
}

Result_t debris_log_set_output(FILE* file) {
    if (file && file != stdout && file != stderr && file != stdin) {
        // Test if file is writable
        if (fprintf(file, "") < 0) {
            return RESULT_PERMISSION_DENIED;
        }
    }
    
    log_config.output_file = file;
    return RESULT_OK;
}

void debris_log_set_structured(bool enabled) {
    log_config.structured = enabled;
}

void debris_log_set_level(ErrorSeverity_t level) {
    log_config.min_level = level;
}

const char* result_to_string(Result_t code) {
    switch (code) {
        case RESULT_OK: return "OK";
        case RESULT_OK_WITH_WARNING: return "OK with warning";
        
        // Input errors
        case RESULT_NULL_POINTER: return "Null pointer";
        case RESULT_INVALID_ARGUMENT: return "Invalid argument";
        case RESULT_OUT_OF_RANGE: return "Out of range";
        case RESULT_INVALID_FORMAT: return "Invalid format";
        case RESULT_FILE_NOT_FOUND: return "File not found";
        case RESULT_PERMISSION_DENIED: return "Permission denied";
        case RESULT_END_OF_FILE: return "End of file";
        case RESULT_INVALID_CHECKSUM: return "Invalid checksum";
        
        // Numerical errors
        case RESULT_NUMERICAL_ERROR: return "Numerical error";
        case RESULT_DIVERGENCE: return "Divergence";
        case RESULT_SINGULAR_MATRIX: return "Singular matrix";
        case RESULT_MAX_ITERATIONS: return "Maximum iterations exceeded";
        case RESULT_TOLERANCE_NOT_MET: return "Tolerance not met";
        case RESULT_OVERFLOW: return "Overflow";
        case RESULT_UNDERFLOW: return "Underflow";
        case RESULT_NAN_DETECTED: return "NaN detected";
        case RESULT_INF_DETECTED: return "Infinity detected";
        
        // Resource errors
        case RESULT_OUT_OF_MEMORY: return "Out of memory";
        case RESULT_OUT_OF_GPU_MEMORY: return "Out of GPU memory";
        case RESULT_DISK_FULL: return "Disk full";
        case RESULT_TOO_MANY_FILES: return "Too many open files";
        case RESULT_RESOURCE_BUSY: return "Resource busy";
        case RESULT_NETWORK_ERROR: return "Network error";
        case RESULT_TIMEOUT: return "Timeout";
        
        // Physics errors
        case RESULT_INVALID_ORBIT: return "Invalid orbit";
        case RESULT_ORBIT_DECAYED: return "Orbit decayed";
        case RESULT_ESCAPE_TRAJECTORY: return "Escape trajectory";
        case RESULT_SUBORBITAL: return "Suborbital";
        case RESULT_COLLISION_DETECTED: return "Collision detected";
        case RESULT_INVALID_MANEUVER: return "Invalid maneuver";
        case RESULT_OUT_OF_FOV: return "Out of field of view";
        case RESULT_BELOW_HORIZON: return "Below horizon";
        
        // System errors
        case RESULT_NOT_INITIALIZED: return "Not initialized";
        case RESULT_ALREADY_INITIALIZED: return "Already initialized";
        case RESULT_INVALID_STATE: return "Invalid state";
        case RESULT_SHUTTING_DOWN: return "Shutting down";
        case RESULT_CONFIG_ERROR: return "Configuration error";
        
        // Feature errors
        case RESULT_NOT_IMPLEMENTED: return "Not implemented";
        case RESULT_DEPRECATED: return "Deprecated";
        case RESULT_UNSUPPORTED: return "Unsupported";
        
        // Integrator errors
        case RESULT_INTEGRATOR_STUCK: return "Integrator stuck";
        case RESULT_STIFF_SYSTEM: return "Stiff system";
        case RESULT_SMALL_STEPSIZE: return "Step size too small";
        case RESULT_SYMPLECTIC_BREAKDOWN: return "Symplectic breakdown";
        
        // Force model errors
        case RESULT_MODEL_OUT_OF_RANGE: return "Model out of range";
        case RESULT_ATMOSPHERE_UNAVAILABLE: return "Atmosphere data unavailable";
        case RESULT_EPHEMERIS_UNAVAILABLE: return "Ephemeris data unavailable";
        
        // Cloud/GMM errors
        case RESULT_GMM_NOT_CONVERGED: return "GMM not converged";
        case RESULT_INVALID_CLOUD: return "Invalid cloud";
        case RESULT_TOO_MANY_COMPONENTS: return "Too many GMM components";
        
        // Conjunction errors
        case RESULT_NO_CONJUNCTION: return "No conjunction";
        case RESULT_AMBIGUOUS_TCA: return "Ambiguous time of closest approach";
        case RESULT_COVARIANCE_INVALID: return "Covariance invalid";
        
        // Internal errors
        case RESULT_INTERNAL_ERROR: return "Internal error";
        case RESULT_ASSERTION_FAILED: return "Assertion failed";
        case RESULT_INVARIANT_VIOLATED: return "Invariant violated";
        case RESULT_UNREACHABLE_CODE: return "Unreachable code";
        
        // Resilience errors
        case RESULT_DEGRADED_ACCURACY: return "Degraded accuracy";
        case RESULT_USING_FALLBACK: return "Using fallback";
        case RESULT_CACHE_MISS: return "Cache miss";
        case RESULT_RECOVERED_ERROR: return "Recovered error";
        
        // Validation errors
        case RESULT_VALIDATION_FAILED: return "Validation failed";
        case RESULT_TEST_FAILED: return "Test failed";
        case RESULT_BENCHMARK_FAILED: return "Benchmark failed";
        
        // Threading errors
        case RESULT_THREAD_ERROR: return "Thread error";
        case RESULT_DEADLOCK_DETECTED: return "Deadlock detected";
        case RESULT_RACE_CONDITION: return "Race condition";
        
        // I/O errors
        case RESULT_SERIALIZATION_ERROR: return "Serialization error";
        case RESULT_DESERIALIZATION_ERROR: return "Deserialization error";
        case RESULT_VERSION_MISMATCH: return "Version mismatch";
        
        // Unknown
        case RESULT_UNKNOWN_ERROR: return "Unknown error";
        
        default: return "Unrecognized error code";
    }
}

const char* result_get_category(Result_t code) {
    if (code == RESULT_OK || code == RESULT_OK_WITH_WARNING) {
        return "SUCCESS";
    } else if (code >= 1000 && code <= 1099) {
        return "INPUT";
    } else if (code >= 1100 && code <= 1199) {
        return "NUMERICAL";
    } else if (code >= 1200 && code <= 1299) {
        return "RESOURCE";
    } else if (code >= 1300 && code <= 1399) {
        return "PHYSICS";
    } else if (code >= 1400 && code <= 1499) {
        return "SYSTEM";
    } else if (code >= 1500 && code <= 1599) {
        return "FEATURE";
    } else if (code >= 1600 && code <= 1699) {
        return "INTEGRATOR";
    } else if (code >= 1700 && code <= 1799) {
        return "FORCE_MODEL";
    } else if (code >= 1800 && code <= 1899) {
        return "CLOUD";
    } else if (code >= 1900 && code <= 1999) {
        return "CONJUNCTION";
    } else if (code >= 2000 && code <= 2099) {
        return "INTERNAL";
    } else if (code >= 2100 && code <= 2199) {
        return "RESILIENCE";
    } else if (code >= 2200 && code <= 2299) {
        return "VALIDATION";
    } else if (code >= 2300 && code <= 2399) {
        return "THREADING";
    } else if (code >= 2400 && code <= 2499) {
        return "IO";
    } else {
        return "UNKNOWN";
    }
}

bool result_is_success(Result_t code) {
    return code == RESULT_OK || code == RESULT_OK_WITH_WARNING;
}

bool result_is_fatal(Result_t code) {
    switch (code) {
        case RESULT_OUT_OF_MEMORY:
        case RESULT_DISK_FULL:
        case RESULT_INTERNAL_ERROR:
        case RESULT_ASSERTION_FAILED:
        case RESULT_UNREACHABLE_CODE:
            return true;
        default:
            return false;
    }
}

bool result_is_recoverable(Result_t code) {
    // Most errors are recoverable except fatal ones
    return !result_is_fatal(code);
}

const ErrorContext* debris_get_last_error(void) {
    return thread_error_context;
}

void debris_clear_error_context(void) {
    free_error_context(thread_error_context);
    thread_error_context = NULL;
}

void debris_set_error_handler(void (*handler)(Result_t, const ErrorContext*)) {
    log_config.error_handler = handler;
}

void debris_handle_fatal_error(void) {
    // Log emergency message
    fprintf(stderr, "\n=== FATAL ERROR - DEBRIS SYSTEM TERMINATING ===\n");
    
    // Print error chain if available
    if (thread_error_context) {
        const ErrorContext* ctx = thread_error_context;
        int depth = 0;
        
        fprintf(stderr, "Error chain (most recent first):\n");
        while (ctx && depth < MAX_ERROR_DEPTH) {
            char timestamp[64];
            format_timestamp(timestamp, sizeof(timestamp), ctx->timestamp_ns);
            
            fprintf(stderr, "  [%d] %s: %s\n", depth, result_to_string(ctx->code), ctx->message);
            fprintf(stderr, "      at %s:%d in %s\n", ctx->file, ctx->line, ctx->function);
            fprintf(stderr, "      thread: %u, time: %s\n", ctx->thread_id, timestamp);
            
            ctx = ctx->previous;
            depth++;
        }
    }
    
    // Attempt to flush all output
    fflush(stdout);
    fflush(stderr);
    
    // Terminate
    abort();
}

Result_t debris_error_init(void) {
    if (log_config.initialized) {
        LOG_WARNING("Error handling already initialized");
        return RESULT_ALREADY_INITIALIZED;
    }
    
    // Initialize configuration
    log_config.output_file = NULL;
    log_config.structured = false;
    log_config.min_level = SEVERITY_INFO;
    log_config.error_handler = NULL;
    log_config.memory_checking_enabled = false;
    
    // Initialize degradation levels
    for (int i = 0; degradation_levels[i].subsystem; i++) {
        degradation_levels[i].level = DEGRADATION_NONE;
    }
    
    log_config.initialized = true;
    LOG_INFO("Error handling subsystem initialized");
    
    return RESULT_OK;
}

void debris_error_shutdown(void) {
    if (!log_config.initialized) {
        return;
    }
    
    LOG_INFO("Error handling subsystem shutting down");
    
    // Clear all thread error contexts
    debris_clear_error_context();
    
    // Flush output
    if (log_config.output_file && log_config.output_file != stderr) {
        fflush(log_config.output_file);
    }
    
    log_config.initialized = false;
}

DegradationLevel_t debris_get_degradation_level(const char* subsystem) {
    for (int i = 0; degradation_levels[i].subsystem; i++) {
        if (strcmp(degradation_levels[i].subsystem, subsystem) == 0) {
            return degradation_levels[i].level;
        }
    }
    
    // Subsystem not found, create new entry
    for (int i = 0; i < sizeof(degradation_levels)/sizeof(degradation_levels[0]); i++) {
        if (degradation_levels[i].subsystem == NULL) {
            degradation_levels[i].subsystem = subsystem;
            degradation_levels[i].level = DEGRADATION_NONE;
            return DEGRADATION_NONE;
        }
    }
    
    // Table full
    return DEGRADATION_NONE;
}

void debris_set_degradation_level(const char* subsystem, DegradationLevel_t level) {
    for (int i = 0; i < sizeof(degradation_levels)/sizeof(degradation_levels[0]); i++) {
        if (degradation_levels[i].subsystem == NULL) {
            // End of table
            break;
        }
        
        if (strcmp(degradation_levels[i].subsystem, subsystem) == 0) {
            DegradationLevel_t old = degradation_levels[i].level;
            degradation_levels[i].level = level;
            
            if (old != level) {
                LOG_WARNING("Degradation level for %s changed: %d -> %d", 
                           subsystem, old, level);
            }
            return;
        }
    }
    
    // Subsystem not found, add to first empty slot
    for (int i = 0; i < sizeof(degradation_levels)/sizeof(degradation_levels[0]); i++) {
        if (degradation_levels[i].subsystem == NULL) {
            degradation_levels[i].subsystem = subsystem;
            degradation_levels[i].level = level;
            LOG_INFO("Added degradation tracking for %s at level %d", subsystem, level);
            return;
        }
    }
    
    LOG_ERROR("Degradation level table full, cannot add %s", subsystem);
}

bool debris_can_proceed(const char* subsystem, DegradationLevel_t required_level) {
    DegradationLevel_t current = debris_get_degradation_level(subsystem);
    return current <= required_level;
}

Result_t debris_check_memory(const void* ptr, size_t size) {
    if (!log_config.memory_checking_enabled) {
        return RESULT_OK;
    }
    
    if (ptr == NULL) {
        return RESULT_NULL_POINTER;
    }
    
    // Simple checks - could be expanded with guard pages, etc.
    // Check for obviously bad pointers (misaligned)
    if (((uintptr_t)ptr & 0x7) != 0) {  // Check 8-byte alignment
        LOG_WARNING("Pointer %p not 8-byte aligned", ptr);
        return RESULT_INVALID_ARGUMENT;
    }
    
    // Note: More sophisticated memory checking would require
    // platform-specific code or libraries like AddressSanitizer
    
    return RESULT_OK;
}

void debris_memory_checking_enable(bool enabled) {
    log_config.memory_checking_enabled = enabled;
    LOG_INFO("Memory checking %s", enabled ? "enabled" : "disabled");
}