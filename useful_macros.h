#ifndef DEF_USEFUL_MACROS
#define DEF_USEFUL_MACROS

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

#endif
