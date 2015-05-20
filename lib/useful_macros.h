#ifndef DEF_USEFUL_MACROS
#define DEF_USEFUL_MACROS

#include <sys/time.h>

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
     fprintf(stderr,"usage of a null pointer in function %s", _f);  \
     exit(EXIT_FAILURE);                                            \
}

#define FREE_PTR(_p) do{ \
	 free( (_p) ); \
     (_p) = NULL; \
}while(0)

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
