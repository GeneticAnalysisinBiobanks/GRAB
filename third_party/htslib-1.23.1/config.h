/* config.h — minimal htslib configuration for GRAB
 *
 * Only zlib (required) is enabled.  No bz2, lzma, curl, GCS, S3, plugins.
 * This file is NOT auto-generated; it replaces the autoconf config.h.
 */

#define HAVE_LIBZ 1

/* We bundle htscodecs inside htslib (not external) */
/* #undef HAVE_EXTERNAL_LIBHTSCODECS */

/* POSIX features available on MinGW, macOS, Linux */
#ifdef _WIN32
/* MinGW: drand48 family not available */
/* #undef HAVE_DRAND48 */
#else
#  define HAVE_DRAND48 1
#endif

#ifdef _WIN32
/* MinGW-w64 */
#  define HAVE_STRUCT_TIMESPEC 1
#else
/* Linux / macOS */
#  define HAVE_FDATASYNC 1
#  define HAVE_FSYNC 1
#  define HAVE_GETPAGESIZE 1
#  define HAVE_GMTIME_R 1
#  define HAVE_MMAP 1
#endif

/* x86 SIMD auto-detection (__get_cpuid_max / __cpuid_count) */
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#  define HAVE_DECL___GET_CPUID_MAX 1
#  define HAVE_DECL___CPUID_COUNT 1
#  ifdef __SSSE3__
#    define HAVE_SSSE3 1
#  endif
#  ifdef __POPCNT__
#    define HAVE_POPCNT 1
#  endif
#  ifdef __SSE4_1__
#    define HAVE_SSE4_1 1
#  endif
#  ifdef __AVX2__
#    define HAVE_AVX2 1
#  endif
#  ifdef __AVX512F__
#    define HAVE_AVX512 1
#  endif
#endif

/* GCC / Clang built-in prefetch */
#if defined(__GNUC__) || defined(__clang__)
#  define HAVE_BUILTIN_PREFETCH 1
#endif

/* Version string placeholder */
#define PACKAGE_VERSION "1.23.1"
