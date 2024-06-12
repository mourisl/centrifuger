#ifndef _LSONG_SIMPLE_VECTOR_HEADER
#define _LSONG_SIMPLE_VECTOR_HEADER

// A light version of vector, which increase the size of the array by 
// a value no more than specified if it got overflow.
// And the type of elements is basic.

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//const int maxInc = -1 ;
// size_type determines the type for represent size, e.g int, size_t, or even short. This can give the ability to control some overhead.
template <class T, class size_type = int>
class SimpleVector
{
private:
	size_type size ;
	size_type capacity ;
	int maxInc ; // The maximal value we can use to increase the capacity.
	int inc ;
	T *s ;
public:
	SimpleVector() : maxInc( -1 )
	{ 
		s = NULL ;
		size = capacity = 0 ;
		inc = 1 ;
	}
	
	SimpleVector( int mi ): maxInc( mi ) 
	{ 
		s = NULL ;
		size = capacity = 0 ;
		inc = 1 ;
	}

	SimpleVector( const SimpleVector &in )
	{
		size = in.size ;
		capacity = in.capacity ;
		if ( capacity > 0 )
		{
			//s = in.s ;
			//if ( in.s == NULL )
			//	printf( "null s. %d %d\n", in.size, in.capacity ) ;
			s = (T *)malloc( sizeof( T ) * capacity ) ;
			memcpy( s, in.s, sizeof( T ) * capacity ) ;
		}
		else 
			s = NULL ;
		inc = in.inc ;
		maxInc = in.maxInc ;
	}
	
	SimpleVector& operator=( const SimpleVector &in )
	{
		if ( this != &in )
		{
			if ( s != NULL )
				free( s ) ;
			size = in.size ;
			capacity = in.capacity ;

			if ( capacity > 0 )
			{
				//s = in.s ;
				s = (T *)malloc( sizeof( T ) * capacity ) ;
				memcpy( s, in.s, sizeof( T ) * capacity ) ;
			}
			else 
				s = NULL ;

			inc = in.inc ;
			maxInc = in.maxInc ;
		}
		return *this ;
	}

	~SimpleVector()
	{
		if ( s != NULL )
			free( s ) ;
		capacity = 0 ;
		size = 0 ;
	}
	
	void Release()
	{
		if ( s != NULL )
			free( s ) ;
		s = NULL ;
		size = capacity = 0 ;
	}

	void Reserve( int sz )
	{
		if ( s != NULL )
			free( s ) ;
		s = (T *)malloc( sizeof( T ) * sz ) ;
		size = 0 ;
		capacity = sz ;
		inc = sz ;

		if ( maxInc > 0 && inc > maxInc )
			inc = maxInc ;
	}

  size_t GetSpace()
  {
    return sizeof(T) * capacity ;
  }

	size_type PushBack( const T &in )	
	{
		if ( size == capacity )
		{
			//int tmp = capacity ;
			capacity += inc ;
			inc *= 2 ;
			if ( maxInc > 0 && inc > maxInc )
				inc = maxInc ;
			if ( size == 0 )
				s = (T *)malloc( sizeof( T ) * capacity ) ;
			else
				s = (T *)realloc( s, sizeof( T ) * capacity ) ;
			if ( s == NULL ) 
			{
				fprintf( stderr, "%s: Failed to allocate memory.\n", __func__ ) ;
				exit( 1 ) ;
			}
		}
		s[ size ] = in ;
		++size ;
		return size ;
	}

	size_type PushBack( const SimpleVector<T> &in )
	{
		int newsize = size + in.size ;
		if ( newsize > capacity )
		{
			//int tmp = capacity ;
			capacity = newsize + inc ;
			inc *= 2 ;
			if ( maxInc > 0 && inc > maxInc )
				inc = maxInc ;
			if ( size == 0 )
				s = (T *)malloc( sizeof( T ) * capacity ) ;
			else
				s = (T *)realloc( s, sizeof( T ) * capacity ) ;
			if ( s == NULL ) 
			{
				fprintf( stderr, "%s: Failed to allocate memory.\n", __func__ ) ;
				exit( 1 ) ;
			}
		}
		memcpy( s + size, in.s, sizeof( T ) * in.size ) ;
		size = newsize ;
		return size ;
	}

	T PopBack()
	{
		if ( size == 0 )
		{
			fprintf( stderr, "%s: empty array.\n", __func__ ) ;
			exit( 1 ) ;
		}
		--size ;
		return s[size] ;
	}
	
	int GetInc()
	{
		return inc ;
	}

	void SetInc( int in )
	{
		inc = in ;
	}

	void SetMaxInc( int in )
	{
		maxInc = in ;
	}

	int GetMaxInc()
	{
		return maxInc ;
	}

	size_type Size() const
	{
		return size ;
	}

	size_type Resize( size_type s ) 
	{
		size = s ;
		return size ;
	}

	size_type Capacity()
	{
		return capacity ;
	}

	T &Get( size_type i )
	{
		if ( i >= size )
		{
			fprintf( stderr, "%s: Access out of the vector.\n", __func__ ) ;
			exit( 1 ) ;
		}
		return s[i] ;
	}

	T &operator[]( size_type i ) const
	{
		/*if ( i >= size )
		{
			printf( "ERROR\n" ) ;
		}*/
		//assert( i < size ) ;
		/*if ( i >= size )
		{
			fprintf( stderr, "%s: Access out of the vector.\n", __func__ ) ;
			exit( 1 ) ;
		}*/
		return s[i] ;
	}
	
	// Return how many element left.
	size_type Remove( size_type ind )
	{
		int i ;
		if ( ind >= size )
		{
			fprintf( stderr, "%s: Access out of the vector.\n", __func__ ) ;
			exit( 1 ) ;
		}

		//if ( size == 1 )
		//	return 0 ;
		for ( i = ind ; i < size - 1 ; ++i )
			s[i] = s[i + 1] ;
		--size ;
		return size ;
	}

	// Allocate less memory. 
	size_type Shrink()
	{
		if ( size < capacity / 4 )
		{
			capacity /= 2 ;
			inc = capacity ;
			if ( inc > maxInc )
				inc = maxInc ;
			s = (T *)realloc( s, sizeof( T ) * capacity ) ;				
		}
		return capacity ;
	}

	void Clear()
	{
		size = 0 ;
	}

	void QSort( int (*compare)(const void*,const void*) )
	{
		qsort( s, size, sizeof( T ), compare ) ;
	}

	size_type BinarySearch( const T &v )
	{
		int l, r, m ;
		l = 0 ; 
		r = size - 1 ;

		while ( l <= r )
		{
			m = ( l + r ) / 2 ;
			if ( s[m] == v )
				return m ;
			else if ( s[m] < v )
				l = m + 1 ;
			else
				r = m - 1 ;

		}
		return l - 1 ; // Should be between  l - 1 and l
	}

	void Destroy()
	{
		if ( s != NULL )
			free( s ) ;
		s = NULL ;
		size = capacity = 0 ;
		inc = 1 ;
	}

	void Overwrite( const SimpleVector<T> &in )
	{
		if ( s != NULL )
			free( s ) ;
		s = NULL ;
		if ( in.s != NULL )
			s = (T *)malloc( sizeof( T ) * in.capacity ) ;
		size = in.size ;
		capacity = in.capacity ;
		inc = in.inc ;
		size_type i ;
		for ( i = 0 ; i < size ; ++i )
			s[i] = in.s[i] ;
	}

	void Reverse()
	{
		size_type i, j ;
		T tmp ;
		for ( i = 0, j = size - 1 ; i < j ; ++i, --j )
		{
			tmp = s[j] ;
			s[j] = s[i] ;
			s[i] = tmp ;
		}
	}
	
	// Expand the array by given size.
	// Does not care about the value in the new allocated space.
	size_type ExpandBy( size_type expandSize )
	{
	  size_type newSize = size + expandSize ;
		if ( newSize <= capacity )
		{
			size = newSize ;
		}
		else
		{
			//int tmp = capacity ;
			capacity = newSize + inc ;
			inc *= 2 ;
			if ( maxInc > 0 && inc > maxInc )
				inc = maxInc ;
			if ( size == 0 )
				s = (T *)malloc( sizeof( T ) * capacity ) ;
			else
				s = (T *)realloc( s, sizeof( T ) * capacity ) ;
			if ( s == NULL ) 
			{
				fprintf( stderr, "%s: Failed to allocate memory.\n", __func__ ) ;
				exit( 1 ) ;
			}
			size = newSize ;
		}
		return size ;
	}

  size_type ExpandTo( size_type newSize )
	{
		return ExpandBy( newSize - size ) ;
	}

	void ShiftRight( size_type shift )
	{
		size = ExpandBy( shift ) ;
		size_type i ;

		for ( i = size - 1 ; i >= shift ; --i )
			s[i] = s[i - shift] ;
		return ;
	}
	
	// Set the content to zero in the range
	void SetZero( size_type start, size_type len )
	{
		memset( s + start, 0, sizeof( T ) * len ) ;
	}

	T *BeginAddress()
	{
		return s ;
	}

	T *EndAddress() 
	{
		return s + size ;
	}
} ;

#endif
