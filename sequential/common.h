/* Xchem: An interface to the ERD integrals package. */

/* Copyright (C) Aftab Patel, Xing Liu  */

/* See ../COPYRIGHT and ../LISENCE */


#ifndef __COMMON_H__
#define __COMMON_H__


typedef signed long long i64;
typedef unsigned long long u64;
typedef signed int i32;
typedef unsigned int u32;
typedef signed char i8;
typedef unsigned char u8;
typedef unsigned short u16;
typedef signed short i16;


#define MIN(a, b)    ((a) < (b) ? (a) : (b))
#define MAX(a, b)    ((a) > (b) ? (a) : (b))

#define _DEBUG_LEVEL_    1   //  0 to 10,
                             //  0 is no debug print info at all,
                             //  10 is full info

#if ( _DEBUG_LEVEL_ == -1 )
#define DPRINTF( level, fmt, args... )        {}
#else
#define DPRINTF( level, fmt, args... )                              \
        do                                                          \
        {                                                           \
            if ( (unsigned)(level) <= _DEBUG_LEVEL_ )               \
            {                                                       \
                fprintf( stdout, fmt, ##args );                     \
                fflush( stdout );                                   \
            }                                                       \
        } while ( 0 )
#endif


#endif /* __COMMON_H__ */
