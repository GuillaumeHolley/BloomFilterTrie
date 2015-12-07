#ifndef DEF_USEFUL_MACROS
#define DEF_USEFUL_MACROS

#include <sys/time.h>

#if (defined(__x86_64__) || defined(_M_X64) || defined(_WIN64) \
  || defined(__powerpc64__) || defined(__ppc64__) || defined(__PPC64__) \
  || defined(__64BIT__) || defined(_LP64) || defined(__LP64__) \
  || defined(__ia64) || defined(__itanium__) || defined(_M_IA64) )
#  define _WORDx64
#else
#  define _WORDx86
#endif

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define CEIL(a, b) (((a) / (b)) + (((a) % (b)) > 0 ? 1 : 0))

#define IS_ODD( num )   ((num) & 1)

#define IS_EVEN( num )  (!IS_ODD( (num) ))

#define ERROR( error_message ) {                                    \
    fprintf (stderr, error_message);                                \
    exit (EXIT_FAILURE);                                            \
}

#define ASSERT_NULL_PTR(_p,_f)                                      \
if( (_p) == NULL )                                                  \
{                                                                   \
     fprintf(stderr,"usage of a null pointer in function %s \n", _f);  \
     exit(EXIT_FAILURE);                                            \
}

typedef struct{
    uint64_t elem1;
    int elem2;
} Duo;

inline void sum_up_duos(Duo* d1, Duo* d2){

    ASSERT_NULL_PTR(d1, "sum_up_duos()")
    ASSERT_NULL_PTR(d2, "sum_up_duos()")

    d1->elem1 += d2->elem1;
    d1->elem2 += d2->elem2;

    return;
}

inline void time_spent(struct timeval *start_time, struct timeval *end_time, struct timeval *resulting_time){

    ASSERT_NULL_PTR(start_time, "time_spent()")
    ASSERT_NULL_PTR(end_time, "time_spent()")
    ASSERT_NULL_PTR(resulting_time, "time_spent()")

    resulting_time->tv_sec = end_time->tv_sec - start_time->tv_sec;
    resulting_time->tv_usec = end_time->tv_usec - start_time->tv_usec;

    if (resulting_time->tv_usec < 0) {
        --resulting_time->tv_sec;
        resulting_time->tv_usec += 1000000;
    }

    return;
}

#endif
