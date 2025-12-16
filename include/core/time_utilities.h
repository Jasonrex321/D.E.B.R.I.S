/**
 * @file time_utilities.h
 * @brief Time and date utilities for DEBRIS system
 * 
 * This module provides high-precision time handling for astrodynamics:
 * - Julian Date (JD) and Modified Julian Date (MJD) operations
 * - Time scale conversions (UTC, TAI, TT, TDB, GPS)
 * - Leap second handling with IERS data
 * - Calendar date conversions (Gregorian, Julian, Year/Day-of-year)
 * - Time arithmetic and comparison
 * - Earth rotation and sidereal time calculations
 * 
 * @section precision Time Precision
 * - All times stored as double precision floating point
 * - Julian Dates: days since 4713 BC, fractional part represents time of day
 * - Internal precision: ~2.3e-10 days = ~20 microseconds
 * - For sub-microsecond precision, use two doubles (high/low parts)
 * 
 * @section timescales Time Scales
 * - UTC: Coordinated Universal Time (with leap seconds)
 * - TAI: International Atomic Time (no leap seconds)
 * - TT: Terrestrial Time (TDT) = TAI + 32.184s
 * - TDB: Barycentric Dynamical Time (relativistic)
 * - GPS: GPS Time = TAI - 19s (constant offset)
 * 
 * @author DEBRIS Engineering Team
 * @date 2024
 * @version 1.0
 * @copyright MIT License
 */

#ifndef TIME_UTILITIES_H
#define TIME_UTILITIES_H

#include "core/error_handling.h"
#include <stdbool.h>
#include <stdint.h>

/* ============================================================================
 * CONSTANTS & MACROS
 * ========================================================================== */

/** Modified Julian Date of J2000 epoch: 2000-01-01 12:00:00 TT */
#define MJD_J2000 51544.5

/** Julian Date of J2000 epoch: 2451545.0 */
#define JD_J2000 2451545.0

/** Julian Date of 1950 epoch: 2433282.5 */
#define JD_1950 2433282.5

/** Days per Julian century: 36525.0 */
#define DAYS_PER_JULIAN_CENTURY 36525.0

/** Days per Julian year: 365.25 */
#define DAYS_PER_JULIAN_YEAR 365.25

/** Seconds per day: 86400.0 */
#define SECONDS_PER_DAY 86400.0

/** Seconds per hour: 3600.0 */
#define SECONDS_PER_HOUR 3600.0

/** Minutes per day: 1440.0 */
#define MINUTES_PER_DAY 1440.0

/** Hours per day: 24.0 */
#define HOURS_PER_DAY 24.0

/** TT - TAI offset in seconds: 32.184 */
#define TT_MINUS_TAI 32.184

/** TAI - UTC offset at J2000 epoch: 32.0 seconds */
#define TAI_MINUS_UTC_J2000 32.0

/** GPS - TAI offset: -19.0 seconds */
#define GPS_MINUS_TAI -19.0

/** Maximum year handled by routines: 2200 */
#define MAX_YEAR 2200

/** Minimum year handled by routines: 1900 */
#define MIN_YEAR 1900

/** Microseconds per day: 86400000000.0 */
#define MICROSECONDS_PER_DAY 86400000000.0

/** Nanoseconds per day: 8.64e13 */
#define NANOSECONDS_PER_DAY 8.64e13

/** Earth rotation rate in radians per second: 7.2921151467e-5 */
#define EARTH_ROTATION_RATE 7.2921151467e-5

/** Earth rotation rate in degrees per second: 4.178074e-3 */
#define EARTH_ROTATION_RATE_DEG 4.178074e-3

/* ============================================================================
 * TIME SCALE ENUMERATIONS
 * ========================================================================== */

/**
 * @brief Time scales supported by DEBRIS
 * 
 * Important: All internal calculations use TT (Terrestrial Time)
 * unless otherwise specified.
 */
typedef enum {
    TIME_SCALE_UTC,   ///< Coordinated Universal Time (with leap seconds)
    TIME_SCALE_TAI,   ///< International Atomic Time
    TIME_SCALE_TT,    ///< Terrestrial Time (TDT) = TAI + 32.184s
    TIME_SCALE_TDB,   ///< Barycentric Dynamical Time
    TIME_SCALE_GPS,   ///< GPS Time = TAI - 19s
    TIME_SCALE_UT1,   ///< Universal Time (Earth rotation)
    TIME_SCALE_TCG,   ///< Geocentric Coordinate Time
    TIME_SCALE_TCB    ///< Barycentric Coordinate Time
} TimeScale_t;

/**
 * @brief Date representation systems
 */
typedef enum {
    DATE_SYSTEM_GREGORIAN,  ///< Gregorian calendar (after 1582-10-15)
    DATE_SYSTEM_JULIAN,     ///< Julian calendar (before 1582-10-04)
    DATE_SYSTEM_MIXED       ///< Mixed Julian/Gregorian (NASA convention)
} DateSystem_t;

/* ============================================================================
 * DATA STRUCTURES
 * ========================================================================== */

/**
 * @brief Calendar date and time structure
 * 
 * Stores date and time with microsecond precision.
 * All fields are integers for exact representation.
 */
typedef struct {
    int16_t year;        ///< Year (e.g., 2024)
    int8_t month;        ///< Month 1-12
    int8_t day;          ///< Day 1-31
    int8_t hour;         ///< Hour 0-23
    int8_t minute;       ///< Minute 0-59
    int8_t second;       ///< Second 0-60 (60 for leap second)
    int32_t microsecond; ///< Microsecond 0-999999
    TimeScale_t scale;   ///< Time scale of this date
} CalendarDateTime_t;

/**
 * @brief Year and day-of-year representation
 */
typedef struct {
    int16_t year;        ///< Year
    int16_t day_of_year; ///< Day of year 1-366
    int32_t second_of_day; ///< Seconds since midnight 0-86399
    int32_t microsecond; ///< Microsecond 0-999999
    TimeScale_t scale;   ///< Time scale
} YearDayTime_t;

/**
 * @brief Two-part Julian Date for high precision
 * 
 * For very high precision applications, split Julian Date
 * into two doubles to minimize rounding errors.
 */
typedef struct {
    double days;    ///< Integer days (or high part)
    double fraction;///< Fractional day (0.0 <= fraction < 1.0)
} TwoPartJD_t;

/**
 * @brief Leap second record
 * 
 * Stores a single leap second event from IERS.
 */
typedef struct {
    double mjd;     ///< Modified Julian Date when leap second occurs
    int8_t delta;   ///< Leap second value (+1 or -1)
    bool is_utc;    ///< True if this is a UTC leap second
} LeapSecond_t;

/**
 * @brief Earth orientation parameters (simplified)
 * 
 * For precise UT1-UTC and polar motion.
 */
typedef struct {
    double mjd;         ///< Modified Julian Date
    double ut1_minus_utc; ///< UT1 - UTC in seconds
    double x_pole;      ///< Polar motion x (arcseconds)
    double y_pole;      ///< Polar motion y (arcseconds)
    double lod;         ///< Length of day (seconds)
} EarthOrientation_t;

/**
 * @brief Time conversion context
 * 
 * Maintains state for time conversions including
 * leap second table and Earth orientation parameters.
 */
typedef struct {
    LeapSecond_t* leap_seconds;   ///< Array of leap second records
    size_t leap_count;            ///< Number of leap seconds
    EarthOrientation_t* eop_data; ///< Earth orientation parameters
    size_t eop_count;             ///< Number of EOP records
    bool auto_update;             ///< True to auto-update from IERS
    char iers_data_path[256];     ///< Path to IERS data files
} TimeContext_t;

/* ============================================================================
 * BASIC TIME CONVERSIONS
 * ========================================================================== */

/**
 * @brief Convert calendar date to Julian Date
 * 
 * Uses NASA/JPL algorithm valid for dates from 4713 BC to 9999 AD.
 * 
 * @param year Year (negative for BC, positive for AD)
 * @param month Month 1-12
 * @param day Day 1-31
 * @param hour Hour 0-23
 * @param minute Minute 0-59
 * @param second Second 0-59.999...
 * @param jd Output Julian Date
 * @return Result_t RESULT_OK on success
 */
Result_t calendar_to_jd(int16_t year, int8_t month, int8_t day,
                       int8_t hour, int8_t minute, double second,
                       double* jd);

/**
 * @brief Convert Julian Date to calendar date
 * 
 * Inverse of calendar_to_jd().
 * 
 * @param jd Julian Date
 * @param year Output year
 * @param month Output month
 * @param day Output day
 * @param hour Output hour
 * @param minute Output minute
 * @param second Output second
 * @return Result_t RESULT_OK on success
 */
Result_t jd_to_calendar(double jd,
                       int16_t* year, int8_t* month, int8_t* day,
                       int8_t* hour, int8_t* minute, double* second);

/**
 * @brief Convert calendar date structure to Julian Date
 * 
 * @param calendar Input calendar date
 * @param jd Output Julian Date
 * @return Result_t RESULT_OK on success
 */
Result_t calendar_datetime_to_jd(const CalendarDateTime_t* calendar, double* jd);

/**
 * @brief Convert Julian Date to calendar date structure
 * 
 * @param jd Input Julian Date
 * @param calendar Output calendar date (must be allocated)
 * @return Result_t RESULT_OK on success
 */
Result_t jd_to_calendar_datetime(double jd, CalendarDateTime_t* calendar);

/**
 * @brief Convert Julian Date to Modified Julian Date
 * 
 * MJD = JD - 2400000.5
 * 
 * @param jd Julian Date
 * @return double Modified Julian Date
 */
double jd_to_mjd(double jd);

/**
 * @brief Convert Modified Julian Date to Julian Date
 * 
 * JD = MJD + 2400000.5
 * 
 * @param mjd Modified Julian Date
 * @return double Julian Date
 */
double mjd_to_jd(double mjd);

/**
 * @brief Split Julian Date into two parts for high precision
 * 
 * Prevents loss of precision when JD is large.
 * 
 * @param jd Julian Date
 * @param two_part Output two-part representation
 */
void jd_to_two_part(double jd, TwoPartJD_t* two_part);

/**
 * @brief Combine two-part Julian Date
 * 
 * @param two_part Two-part representation
 * @return double Combined Julian Date
 */
double two_part_to_jd(const TwoPartJD_t* two_part);

/**
 * @brief Convert day-of-year to calendar date
 * 
 * @param year Year
 * @param day_of_year Day of year (1-366)
 * @param month Output month (1-12)
 * @param day Output day (1-31)
 * @return Result_t RESULT_OK on success
 */
Result_t day_of_year_to_calendar(int16_t year, int16_t day_of_year,
                                int8_t* month, int8_t* day);

/**
 * @brief Convert calendar date to day-of-year
 * 
 * @param year Year
 * @param month Month (1-12)
 * @param day Day (1-31)
 * @return int16_t Day of year (1-366)
 */
int16_t calendar_to_day_of_year(int16_t year, int8_t month, int8_t day);

/* ============================================================================
 * TIME SCALE CONVERSIONS
 * ========================================================================== */

/**
 * @brief Convert time between different scales
 * 
 * Primary conversion function. Converts from any supported
 * time scale to any other.
 * 
 * @param time Input time (in input_scale)
 * @param input_scale Scale of input time
 * @param output_scale Desired output scale
 * @param output_time Output time (in output_scale)
 * @param context Time context for leap seconds and EOP (can be NULL)
 * @return Result_t RESULT_OK on success
 */
Result_t convert_time_scale(double time,
                           TimeScale_t input_scale,
                           TimeScale_t output_scale,
                           double* output_time,
                           const TimeContext_t* context);

/**
 * @brief Convert UTC to TAI (including leap seconds)
 * 
 * @param utc UTC time (Julian Date)
 * @param tai Output TAI time (Julian Date)
 * @param context Time context with leap seconds
 * @return Result_t RESULT_OK on success
 */
Result_t utc_to_tai(double utc, double* tai, const TimeContext_t* context);

/**
 * @brief Convert TAI to UTC (including leap seconds)
 * 
 * @param tai TAI time (Julian Date)
 * @param utc Output UTC time (Julian Date)
 * @param context Time context with leap seconds
 * @return Result_t RESULT_OK on success
 */
Result_t tai_to_utc(double tai, double* utc, const TimeContext_t* context);

/**
 * @brief Convert TT to TDB (approximate)
 * 
 * Uses Moyer's formula for approximate TDB-TT conversion.
 * Accurate to about 10 microseconds for Earth orbit.
 * 
 * @param tt TT time (Julian Date)
 * @param tdb Output TDB time (Julian Date)
 * @return Result_t RESULT_OK on success
 */
Result_t tt_to_tdb(double tt, double* tdb);

/**
 * @brief Convert TDB to TT (approximate)
 * 
 * Inverse of tt_to_tdb().
 * 
 * @param tdb TDB time (Julian Date)
 * @param tt Output TT time (Julian Date)
 * @return Result_t RESULT_OK on success
 */
Result_t tdb_to_tt(double tdb, double* tt);

/**
 * @brief Convert UTC to TT (common operation)
 * 
 * Combines utc_to_tai() and tai_to_tt().
 * 
 * @param utc UTC time (Julian Date)
 * @param tt Output TT time (Julian Date)
 * @param context Time context with leap seconds
 * @return Result_t RESULT_OK on success
 */
Result_t utc_to_tt(double utc, double* tt, const TimeContext_t* context);

/**
 * @brief Convert TT to UTC (common operation)
 * 
 * Combines tt_to_tai() and tai_to_utc().
 * 
 * @param tt TT time (Julian Date)
 * @param utc Output UTC time (Julian Date)
 * @param context Time context with leap seconds
 * @return Result_t RESULT_OK on success
 */
Result_t tt_to_utc(double tt, double* utc, const TimeContext_t* context);

/**
 * @brief Get current UTC time as Julian Date
 * 
 * Uses system clock to get current time.
 * 
 * @param jd_utc Output Julian Date (UTC)
 * @return Result_t RESULT_OK on success
 */
Result_t get_current_utc_jd(double* jd_utc);

/**
 * @brief Get current TT time as Julian Date
 * 
 * @param jd_tt Output Julian Date (TT)
 * @param context Time context for leap seconds (can be NULL)
 * @return Result_t RESULT_OK on success
 */
Result_t get_current_tt_jd(double* jd_tt, const TimeContext_t* context);

/* ============================================================================
 * LEAP SECOND HANDLING
 * ========================================================================== */

/**
 * @brief Initialize time context with default leap seconds
 * 
 * Loads built-in leap second table (updated through 2023).
 * For production use, should load from IERS file.
 * 
 * @param context Time context to initialize
 * @return Result_t RESULT_OK on success
 */
Result_t time_context_init(TimeContext_t* context);

/**
 * @brief Initialize time context from IERS file
 * 
 * @param context Time context to initialize
 * @param iers_file Path to IERS leap second file
 * @return Result_t RESULT_OK on success
 */
Result_t time_context_init_from_file(TimeContext_t* context, const char* iers_file);

/**
 * @brief Free resources associated with time context
 * 
 * @param context Time context to destroy
 */
void time_context_destroy(TimeContext_t* context);

/**
 * @brief Get TAI - UTC offset for a given UTC time
 * 
 * @param utc UTC time (Julian Date)
 * @param context Time context with leap seconds
 * @param delta_tai_utc Output TAI-UTC in seconds
 * @return Result_t RESULT_OK on success
 */
Result_t get_tai_minus_utc(double utc, const TimeContext_t* context, double* delta_tai_utc);

/**
 * @brief Check if a UTC time is during a leap second
 * 
 * @param utc UTC time (Julian Date)
 * @param context Time context with leap seconds
 * @param is_leap_second Output true if during leap second
 * @return Result_t RESULT_OK on success
 */
Result_t is_leap_second(double utc, const TimeContext_t* context, bool* is_leap_second);

/**
 * @brief Add leap second to the context
 * 
 * For updating leap second table dynamically.
 * 
 * @param context Time context
 * @param mjd MJD of leap second
 * @param delta Leap second value (+1 or -1)
 * @return Result_t RESULT_OK on success
 */
Result_t add_leap_second(TimeContext_t* context, double mjd, int8_t delta);

/**
 * @brief Get number of leap seconds between two UTC times
 * 
 * @param utc1 Start UTC time (Julian Date)
 * @param utc2 End UTC time (Julian Date)
 * @param context Time context with leap seconds
 * @param count Output number of leap seconds
 * @return Result_t RESULT_OK on success
 */
Result_t count_leap_seconds(double utc1, double utc2, const TimeContext_t* context, int32_t* count);

/* ============================================================================
 * TIME ARITHMETIC AND COMPARISON
 * ========================================================================== */

/**
 * @brief Add seconds to a Julian Date
 * 
 * Handles leap seconds correctly when adding to UTC.
 * 
 * @param jd Input Julian Date
 * @param seconds Seconds to add (can be negative)
 * @param scale Time scale of input/output
 * @param context Time context for leap seconds (can be NULL for non-UTC)
 * @param result_jd Output Julian Date
 * @return Result_t RESULT_OK on success
 */
Result_t jd_add_seconds(double jd, double seconds, TimeScale_t scale,
                       const TimeContext_t* context, double* result_jd);

/**
 * @brief Add days to a Julian Date
 * 
 * @param jd Input Julian Date
 * @param days Days to add (can be negative)
 * @param result_jd Output Julian Date
 * @return Result_t RESULT_OK on success
 */
Result_t jd_add_days(double jd, double days, double* result_jd);

/**
 * @brief Difference between two times in seconds
 * 
 * Accounts for leap seconds when times are in UTC.
 * 
 * @param jd1 First Julian Date
 * @param jd2 Second Julian Date
 * @param scale Time scale of both times
 * @param context Time context for leap seconds (can be NULL for non-UTC)
 * @param diff_seconds Output difference (jd2 - jd1) in seconds
 * @return Result_t RESULT_OK on success
 */
Result_t jd_difference_seconds(double jd1, double jd2, TimeScale_t scale,
                              const TimeContext_t* context, double* diff_seconds);

/**
 * @brief Compare two Julian Dates with tolerance
 * 
 * @param jd1 First Julian Date
 * @param jd2 Second Julian Date
 * @param tolerance_seconds Tolerance in seconds
 * @param scale Time scale of both times
 * @param context Time context for leap seconds (can be NULL for non-UTC)
 * @return int -1 if jd1 < jd2, 0 if equal within tolerance, 1 if jd1 > jd2
 */
int jd_compare(double jd1, double jd2, double tolerance_seconds,
              TimeScale_t scale, const TimeContext_t* context);

/**
 * @brief Calculate mean Julian Date (average of two times)
 * 
 * @param jd1 First Julian Date
 * @param jd2 Second Julian Date
 * @param mean_jd Output mean Julian Date
 * @return Result_t RESULT_OK on success
 */
Result_t jd_mean(double jd1, double jd2, double* mean_jd);

/**
 * @brief Linear interpolation between two times
 * 
 * @param jd1 Time 1 (Julian Date)
 * @param value1 Value at time 1
 * @param jd2 Time 2 (Julian Date)
 * @param value2 Value at time 2
 * @param jd_interp Interpolation time (Julian Date)
 * @param value_interp Output interpolated value
 * @return Result_t RESULT_OK on success
 */
Result_t jd_interpolate(double jd1, double value1,
                       double jd2, double value2,
                       double jd_interp, double* value_interp);

/* ============================================================================
 * EARTH ROTATION AND SIDEREAL TIME
 * ========================================================================== */

/**
 * @brief Calculate Greenwich Mean Sidereal Time (GMST)
 * 
 * Based on IAU 1982 model.
 * 
 * @param ut1 UT1 time (Julian Date)
 * @param gmst_radians Output GMST in radians (0 to 2π)
 * @return Result_t RESULT_OK on success
 */
Result_t calculate_gmst(double ut1, double* gmst_radians);

/**
 * @brief Calculate Greenwich Apparent Sidereal Time (GAST)
 * 
 * Includes nutation correction.
 * 
 * @param ut1 UT1 time (Julian Date)
 * @param gast_radians Output GAST in radians (0 to 2π)
 * @return Result_t RESULT_OK on success
 */
Result_t calculate_gast(double ut1, double* gast_radians);

/**
 * @brief Calculate Local Mean Sidereal Time (LMST)
 * 
 * @param ut1 UT1 time (Julian Date)
 * @param longitude_radians Longitude in radians (positive east)
 * @param lmst_radians Output LMST in radians (0 to 2π)
 * @return Result_t RESULT_OK on success
 */
Result_t calculate_lmst(double ut1, double longitude_radians, double* lmst_radians);

/**
 * @brief Calculate Local Apparent Sidereal Time (LAST)
 * 
 * @param ut1 UT1 time (Julian Date)
 * @param longitude_radians Longitude in radians (positive east)
 * @param last_radians Output LAST in radians (0 to 2π)
 * @return Result_t RESULT_OK on success
 */
Result_t calculate_last(double ut1, double longitude_radians, double* last_radians);

/**
 * @brief Calculate Earth rotation angle (ERA)
 * 
 * Based on IAU 2000 model.
 * 
 * @param ut1 UT1 time (Julian Date)
 * @param era_radians Output ERA in radians (0 to 2π)
 * @return Result_t RESULT_OK on success
 */
Result_t calculate_earth_rotation_angle(double ut1, double* era_radians);

/**
 * @brief Convert UT1 to UTC using Earth orientation parameters
 * 
 * @param ut1 UT1 time (Julian Date)
 * @param eop Earth orientation parameters
 * @param utc Output UTC time (Julian Date)
 * @return Result_t RESULT_OK on success
 */
Result_t ut1_to_utc(double ut1, const EarthOrientation_t* eop, double* utc);

/**
 * @brief Convert UTC to UT1 using Earth orientation parameters
 * 
 * @param utc UTC time (Julian Date)
 * @param eop Earth orientation parameters
 * @param ut1 Output UT1 time (Julian Date)
 * @return Result_t RESULT_OK on success
 */
Result_t utc_to_ut1(double utc, const EarthOrientation_t* eop, double* ut1);

/* ============================================================================
 * STRING FORMATTING AND PARSING
 * ========================================================================== */

/**
 * @brief Format Julian Date as ISO 8601 string
 * 
 * Example: "2024-01-15T12:30:45.123456Z"
 * 
 * @param jd Julian Date
 * @param scale Time scale (for suffix: Z for UTC, etc.)
 * @param context Time context for conversions (can be NULL)
 * @param buffer Output buffer
 * @param buffer_size Size of buffer (minimum 32 bytes)
 * @return Result_t RESULT_OK on success
 */
Result_t jd_to_iso8601(double jd, TimeScale_t scale,
                      const TimeContext_t* context,
                      char* buffer, size_t buffer_size);

/**
 * @brief Parse ISO 8601 string to Julian Date
 * 
 * Supports formats:
 * - "2024-01-15T12:30:45Z"
 * - "2024-01-15 12:30:45.123"
 * - "20240115T123045.123456"
 * 
 * @param iso_string ISO 8601 string
 * @param scale Time scale of string (input)
 * @param context Time context for conversions (can be NULL)
 * @param jd Output Julian Date
 * @return Result_t RESULT_OK on success
 */
Result_t iso8601_to_jd(const char* iso_string, TimeScale_t scale,
                      const TimeContext_t* context, double* jd);

/**
 * @brief Format Julian Date as TLE epoch string
 * 
 * Example: "24015.12345678" (YYDDD.dddddddd)
 * 
 * @param jd Julian Date (UTC)
 * @param context Time context for leap seconds
 * @param buffer Output buffer (minimum 16 bytes)
 * @return Result_t RESULT_OK on success
 */
Result_t jd_to_tle_epoch(double jd, const TimeContext_t* context, char* buffer);

/**
 * @brief Parse TLE epoch string to Julian Date
 * 
 * @param tle_epoch TLE epoch string (YYDDD.dddddddd)
 * @param context Time context for leap seconds
 * @param jd Output Julian Date (UTC)
 * @return Result_t RESULT_OK on success
 */
Result_t tle_epoch_to_jd(const char* tle_epoch, const TimeContext_t* context, double* jd);

/**
 * @brief Format time as human-readable string
 * 
 * Example: "2024-Jan-15 12:30:45.123 UTC"
 * 
 * @param jd Julian Date
 * @param scale Time scale
 * @param context Time context (can be NULL)
 * @param buffer Output buffer
 * @param buffer_size Size of buffer
 * @return Result_t RESULT_OK on success
 */
Result_t jd_to_human_readable(double jd, TimeScale_t scale,
                             const TimeContext_t* context,
                             char* buffer, size_t buffer_size);

/* ============================================================================
 * VALIDATION AND UTILITIES
 * ========================================================================== */

/**
 * @brief Validate calendar date
 * 
 * Checks for valid year, month, day, hour, minute, second.
 * 
 * @param year Year
 * @param month Month
 * @param day Day
 * @param hour Hour
 * @param minute Minute
 * @param second Second
 * @return Result_t RESULT_OK if valid
 */
Result_t validate_calendar_date(int16_t year, int8_t month, int8_t day,
                               int8_t hour, int8_t minute, double second);

/**
 * @brief Check if year is a leap year
 * 
 * @param year Year
 * @return true Year is a leap year
 * @return false Year is not a leap year
 */
bool is_leap_year(int16_t year);

/**
 * @brief Get number of days in month
 * 
 * @param year Year (for February leap year calculation)
 * @param month Month 1-12
 * @return int8_t Days in month (28, 29, 30, or 31)
 */
int8_t days_in_month(int16_t year, int8_t month);

/**
 * @brief Calculate day of week (0=Sunday, 6=Saturday)
 * 
 * @param jd Julian Date
 * @return int8_t Day of week (0-6)
 */
int8_t day_of_week(double jd);

/**
 * @brief Calculate Julian century from J2000
 * 
 * T = (JD - 2451545.0) / 36525.0
 * 
 * @param jd Julian Date
 * @return double Julian centuries since J2000
 */
double jd_to_julian_century(double jd);

/**
 * @brief Convert Julian century to Julian Date
 * 
 * JD = T * 36525.0 + 2451545.0
 * 
 * @param t Julian centuries since J2000
 * @return double Julian Date
 */
double julian_century_to_jd(double t);

/**
 * @brief Calculate mean anomaly of Sun (approximate)
 * 
 * Used for some time conversions.
 * 
 * @param t Julian centuries since J2000
 * @return double Mean anomaly in radians
 */
double sun_mean_anomaly(double t);

/**
 * @brief Print time debugging information
 * 
 * @param jd Julian Date
 * @param scale Time scale
 * @param label Label for output
 * @param context Time context (can be NULL)
 */
void print_time_debug(double jd, TimeScale_t scale, const char* label,
                     const TimeContext_t* context);

#endif /* TIME_UTILITIES_H */