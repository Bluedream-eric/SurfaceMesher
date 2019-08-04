#ifndef MATRIXT_H
#define MATRIXT_H

//== INCLUDES =================================================================

#include <OpenMesh\Core\Geometry\VectorT.hh>
#include <iostream>

//== CLASS DEFINITION =========================================================

/** The M by N values of the template Scalar type are the only data members
    of the class MatrixT<Scalar,M,N>. This guarantees 100% compatibility
    with arrays of type Scalar and size M by N, allowing us to define the
    cast operators to and from arrays and array pointers.
*/

template <typename Scalar,int M,int N> struct MatrixDataT
{
    Scalar values_[M][N];
};

#define unroll(expr)        for (int i=0; i<M; ++i) for (int j=0; j<N; ++j) expr(i,j)
// #define unroll_lower(expr)  for (int i=0; i<M; ++i) for (int j=0; j<=i; ++j) expr(i,j)  // especially for square matrix

/** \class MatrixT  MatrixT.h
    A matrix is an bi-array of <M*N> values of type <Scalar>.
    The actual data is stored in an MatrixDataT, this class just adds
    the necessary operators.
*/

template <typename Scalar, int M, int N>
class MatrixT : public MatrixDataT<Scalar, M, N>
{
private:
	typedef MatrixDataT<Scalar, M, N>                      Base;

public:
	//== class info ===========================================================

    // the type of the scalar used in this template
    typedef Scalar value_type;

    // type of this vector
    typedef MatrixT<Scalar,M,N>  matrix_type;

    // returns number of rows and column of the matrix (deprecated)
    static inline int n_rows()     { return M; }
	static inline int n_columns()  { return N; }

    //== constructors & destructor ============================================

	// destructor
	~MatrixT() {}

    // default constructor creates uninitialized values.
    inline MatrixT() {}

    // no 1 by 1 matrix, set all matrix elements as 0.0, seems still need 1 by 1 matrix
    inline MatrixT(const Scalar& _s) 
    {
#define expr(i,j)  Base::values_[i][j] = _s;
		unroll(expr);
#undef expr
    }

	// special constructor for 2 by 1 matrices
	inline MatrixT(const Scalar& _m0, const Scalar& _m1)
	{
		assert(M==2 && N==1);

		Base::values_[0][0] = _m0;
		Base::values_[1][0] = _m1;
	}

    // special constructor for 2 by 2 matrices, row majored
    inline MatrixT(const Scalar& _m0, const Scalar& _m1, const Scalar& _m2, const Scalar& _m3) 
	{
        assert((M==2 && N==2) || (M==4 && N==1));

        if (M==2 && N==2)
		{
			Base::values_[0][0] = _m0; Base::values_[0][1] = _m1;
			Base::values_[1][0] = _m2; Base::values_[1][1] = _m3;
		}
		else
		{
			Base::values_[0][0] = _m0; Base::values_[1][0] = _m1;
			Base::values_[2][0] = _m2; Base::values_[3][0] = _m3;
		}
    }

	// special constructor for 2 by 3 or 3 by 2 matrices, row majored
	MatrixT(const Scalar& _m0, const Scalar& _m1, const Scalar& _m2, const Scalar& _m3, 
		    const Scalar& _m4, const Scalar& _m5) // const and & ?
	{
		assert((M==2 && N ==3) || (M==3 && N==2));

		if (M==2 && N==3) // Scalar m[2][4], belongs to "if-else" name domain
		{
			Base::values_[0][0] = _m0; Base::values_[0][1] = _m1; Base::values_[0][2] = _m2; 
			Base::values_[1][0] = _m3; Base::values_[1][1] = _m4; Base::values_[1][2] = _m5;
		}
		else
		{
			Base::values_[0][0] = _m0; Base::values_[0][1] = _m1; 
			Base::values_[1][0] = _m2; Base::values_[1][1] = _m3; 
			Base::values_[2][0] = _m4; Base::values_[2][1] = _m5;
		}
	}

	// special constructor for 2 by 4 or 4 by 2 matrices, row majored
	MatrixT(const Scalar& _m0, const Scalar& _m1, const Scalar& _m2, const Scalar& _m3, 
		    const Scalar& _m4, const Scalar& _m5, const Scalar& _m6, const Scalar& _m7)
	{
		assert((M==2 && N==4) || (M==4 && N==2));

		if (M==2 && N==4) // Scalar m[2][4], belongs to "if-else" name domain
		{
			Base::values_[0][0] = _m0; Base::values_[0][1] = _m1; 
			Base::values_[0][2] = _m2; Base::values_[0][3] = _m3;
			Base::values_[1][0] = _m4; Base::values_[1][1] = _m5;
			Base::values_[1][2] = _m6; Base::values_[1][3] = _m7;
		}
		else
		{
			Base::values_[0][0] = _m0; Base::values_[0][1] = _m1; 
			Base::values_[1][0] = _m2; Base::values_[1][1] = _m3;
			Base::values_[2][0] = _m4; Base::values_[2][1] = _m5;
			Base::values_[3][0] = _m6; Base::values_[3][1] = _m7;
		}
	}

	// special constructor for 3 by 1 matrices, row majored
	inline MatrixT(const Scalar& _m0, const Scalar& _m1, const Scalar& _m2)
	{
		assert(M==3 && N==1);

		Base::values_[0][0] = _m0;
		Base::values_[1][0] = _m1;
		Base::values_[2][0] = _m2;
	}

	// special constructor for 3 by 3 matrices, row majored
    inline MatrixT(const Scalar& _m0, const Scalar& _m1, const Scalar& _m2, const Scalar& _m3, 
		    const Scalar& _m4, const Scalar& _m5, const Scalar& _m6, const Scalar& _m7, const Scalar& _m8) 
	{
        assert(M==3 && N==3);

		Scalar m[3][3] = {{_m0, _m1, _m2}, {_m3, _m4, _m5}, {_m6, _m7, _m8}};
#define expr(i,j)  Base::values_[i][j] = m[i][j];
		unroll(expr);
#undef expr
    }

	// special constructor for 3 by 4 or 4 by 3 matrices, row majored
    inline MatrixT(const Scalar& _m0, const Scalar& _m1, const Scalar& _m2, const Scalar& _m3, 
		    const Scalar& _m4, const Scalar& _m5, const Scalar& _m6, const Scalar& _m7, 
			const Scalar& _m8, const Scalar& _m9, const Scalar& _m10, const Scalar& _m11) 
	{
        assert((M==3 && N==4) || (M==4 && N==3));

		if (M==3 && N==4)
		{
			Base::values_[0][0] = _m0; Base::values_[0][1] = _m1; 
			Base::values_[0][2] = _m2; Base::values_[0][3] = _m3;
			Base::values_[1][0] = _m4; Base::values_[1][1] = _m5;
			Base::values_[1][2] = _m6; Base::values_[1][3] = _m7;
			Base::values_[2][0] = _m8; Base::values_[2][1] = _m9;
			Base::values_[2][2] = _m10; Base::values_[2][3] = _m11;
		}
		else
		{
			Base::values_[0][0] = _m0; Base::values_[0][1] = _m1; Base::values_[0][2] = _m2; 
			Base::values_[1][0] = _m3; Base::values_[1][1] = _m4; Base::values_[1][2] = _m5;
			Base::values_[2][0] = _m6; Base::values_[2][1] = _m7; Base::values_[2][2] = _m8;
            Base::values_[3][0] = _m0; Base::values_[3][1] = _m10; Base::values_[3][2] = _m11;
		}
    }

	// special constructor for 4 by 4 matrices, row majored
    inline MatrixT(const Scalar& _m0, const Scalar& _m1, const Scalar& _m2, const Scalar& _m3, 
		    const Scalar& _m4, const Scalar& _m5, const Scalar& _m6, const Scalar& _m7, 
			const Scalar& _m8, const Scalar& _m9, const Scalar& _m10, const Scalar& _m11,
			const Scalar& _m12, const Scalar& _m13, const Scalar& _m14, const Scalar& _m15) 
	{
        assert(M==4 && N==4);

		Scalar m[4][4] = {{_m0, _m1, _m2, _m3}, {_m4, _m5, _m6, _m7}, {_m8, _m9, _m10, _m11}, {_m12, _m13, _m14, _m15}};
#define expr(i,j)  Base::values_[i][j] = m[i][j];
		unroll(expr);
#undef expr
    }

	// copy constructor
    inline MatrixT(const matrix_type& _rhs) 
	{
        memcpy(Base::values_, _rhs.Base::values_, M*N*sizeof(Scalar));
    }	

    // construct from a value array (explicit), need to test
    explicit inline MatrixT(const Scalar _values[M][N]) 
	{
        memcpy(Base::values_, _values, M*N*sizeof(Scalar));
    }

    // copy & cast constructor (explicit)
    template<typename otherScalarType>
    explicit inline MatrixT(const MatrixT<otherScalarType,M,N>& _rhs) 
	{
        operator=(_rhs);
    }

	// assignment operator
    inline matrix_type& operator=(const matrix_type& _rhs) 
	{
        memcpy(Base::values_, _rhs.Base::values_, M*N*sizeof(Scalar));
        return *this;
    } 

	// assignment from a matrix with a different scalar type
    template<typename otherScalarType>
    inline matrix_type& operator=(const MatrixT<otherScalarType,M,N>& _rhs) 
	{
#define expr(i,j)  Base::values_[i][j] = (Scalar)_rhs[i][j];
        unroll(expr);
#undef expr
        return *this;
    }

	//== casts ================================================================

	// \attention: may be wrong! In VectorT, it's Scalar*
    // cast to Scalar array, originally return Base::values_;
	// can not be const Scalar**(), otherwise compile errer
    inline operator Scalar**() { return Base::values_; }

	 //== element access ======================================================

	// set/get the _i,_j'th element, read-write, the size_t counter-part
	inline Scalar& operator()(int _i, int _j)
	{
		assert(_i>=0 && _i<M && _j>=0 && _j<N);
		return Base::values_[_i][_j];
	}

	// conflict with the above operator ?
	// get the _i,_j'th element
	inline Scalar operator()(int _i, int _j) const
	{
		assert(_i>=0 && _i<M && _j>=0 && _j<N);
		return Base::values_[_i][_j];
	}

    // get row and column data, subscripting, read-write

	inline Scalar* operator[](const size_t _row) const
	{
        assert(_row < M); 
		return Base::values_[_row];
    }

	inline Scalar* get_row(const size_t _row) const
	{
        assert(_row < M);
		return Base::values_[_row];
	}

	inline void get_row(const size_t _row, Scalar _data[N]) const
	{
		assert(_row < M);
		for (int i = 0; i < N; ++i)
			_data[i] = Base::values_[_row][i];
	}

	inline void get_column(const size_t _column, Scalar _data[M]) const
	{
		assert(_column < N);
		for (int i = 0; i < M; ++i)
			_data[i] = Base::values_[i][_column];
	}

	//== arithmetic +,-,+=,-=,*=,/£½,^= and so on =============================

	// operator +
	inline matrix_type operator+()
	{ return (*this); }

	// operator -
	inline matrix_type operator-()
	{ return (Scalar)1.0 * (*this); }

	// matrix self-addition, with same matrix type
	// \todo_clg: matrix self-addition, with almost the same matrix type
    inline matrix_type& operator+=(const matrix_type& _rhs) 
	{
#define expr(i,j) Base::values_[i][j] += _rhs(i,j);
        unroll(expr);
#undef expr
        return *this;
    }

    // matrix difference from this
	// \todo_clg: matrix self-substraction, with almost the same matrix type
    inline matrix_type& operator-=(const matrix_type& _rhs) 
	{
#define expr(i,j) Base::values_[i][j] -= _rhs(i,j);
        unroll(expr);
#undef expr
        return *this;
    }

    // component-wise self-multiplication with scalar
    inline matrix_type& operator*=(const Scalar& _s) 
	{
#define expr(i,j) Base::values_[i][j] *= _s;
        unroll(expr);
#undef expr
        return *this;
    }

	// self-multiplication, M*N, N*N -> M*N
	// \attention: could be a redefinition with the above function in some cases!
    inline matrix_type& operator*=(const MatrixT<Scalar,N,N>& _m)
	{
        matrix_type tmp_mat = *this;

		int i,j,k;
		Scalar sum;
		for (i=0; i<M; ++i)
			for (j=0; j<N; ++j)
			{
                sum = (Scalar)0;
				for (k=0; k<N; ++k)
					sum += Base::values_[i][k] * _m(k,j);
				(*this)(i,j) = sum;
			}

		return *this;
    }

    // component-wise self-division by scalar
    // \attention: v *= (1/_s) is much faster than this  
    inline matrix_type& operator/=(const Scalar& _s) 
	{
#define expr(i,j) Base::values_[i][j] /= _s;
        unroll(expr);
#undef expr
        return *this;
    }

    // combined power and assignment operator
    inline matrix_type& operator^=(const size_t& _power)
	{
		assert(M==N);

        matrix_type tmp_mat(*this);
        for (size_t i=2; i <= _power; ++i)
        *this = *this * tmp_mat;

		return *this;
    }

	//== misc functions =======================================================
    
	// store the same value in each component (e.g. to clear all entries)
	inline void reset(const Scalar& _s = (Scalar)0.0)
	{
#define expr(i,j) Base::values_[i][j] = _s;
        unroll(expr);
#undef expr
	}

	// set all elements to zero
	inline void null()
	{
        reset((Scalar)0.0);
	}

	// set to unit matrix, it does not have to be square matrix
	inline void unit()
	{
		null();
#define expr(i,j)  Base::values_[i][j] = (i==j ? (Scalar)1.0 : (Scalar)0.0);
        unroll(expr);
#undef expr
	}

	//== utility methods ======================================================

	Scalar min() const
	{
        Scalar min_val = 1e38;
#define expr(i,j)  min_val = (Base::values_[i][j] < min_val ? Base::values_[i][j] : min_val);
        unroll(expr);
#undef expr

		return min_val;
	}

	Scalar max() const
	{
        Scalar max_val = -1e38;
#define expr(i,j)  max_val = (Base::values_[i][j] > max_val ? Base::values_[i][j] : max_val);
        unroll(expr);
#undef expr

		return max_val;
	}

	Scalar mean() const
	{
        Scalar sum = 0;
#define expr(i,j)  sum += Base::values_[i][j];
        unroll(expr);
#undef expr

		return sum / (M*N);
	}

	Scalar Norm()
	{
        Scalar square_sum = 0;
#define expr(i,j)  square_sum += Base::values_[i][j] * Base::values_[i][j];
        unroll(expr);
#undef expr

		return (Scalar)sqrt(square_sum);
	}

	// private partial pivoting method, used in gauss elimination with column majored 
	int pivot(int _ii) const
	{
        int i, k = int(_ii);
        double tmp, abs_max = -1;
        for (i=_ii; i<M; ++i)
			if ((tmp = fabs(Base::values_[i][_ii])) > abs_max && (tmp != 0.0))
            {
                abs_max = tmp;
                k = i;
            }

		if (Base::values_[k][_ii] == double(0))
            return -1;

        if (k != _ii)  // \attention: here we do not exchange rows k and _ii
			return k;

        return 0;  // k == _ii
	}

	Scalar determinant() const
	{
        assert(M==N);

		int i, j, k, index;
        Scalar piv_val = 0, det_val = Scalar(1); 
		Scalar tmp_exch[N];// try to use VectorT<Scalar, N>
        matrix_type tmp_mat(*this);
        for (i=0; i<M; ++i)
		{
			index = tmp_mat.pivot(i);
            if (index == -1)
				return 0; // singular
			if (index != 0)
			{
				det_val = -det_val;

                tmp_mat.get_row(index, tmp_exch);
			    for (j=0; j<N; ++j)  // exchange rows i and index
				    tmp_mat(index,j) = tmp_mat(i,j);
                for (j=0; j<N; ++j)
				    tmp_mat(i,j) = tmp_exch[j];
			}
			det_val *= tmp_mat(i,i);

			for (j=i+1; j<M; ++j)
			{
				piv_val = tmp_mat(j,i) / tmp_mat(i,i);
				for (k=i+1; k<M; ++k)
					tmp_mat(j,k) -= piv_val*tmp_mat(i,k);
			}
		}

		return det_val;
	}

	// different from operator!,this function operates on itself
	// use gauss-jordan elimination with column majored, ref. <<LinChengsen>>, p.85
	matrix_type& inversion()
	{
        assert(M==N);
		// assert(!is_singular());

		int i, j, k, r, piv[M];
		for (i=0; i<M; ++i)
			piv[i] = i;
		Scalar tmp_exch[N];
		// matrix_type res_mat(*this);

		for (i=0; i<M; ++i)
		{
            r = pivot(i); // piv[i]
			if (r==-1)
			{
				assert(false); // or print singular natrix message
				return *this;  // \attention: here *this could have been changed
			}
			if (r!=0)  // the same as r!=i, then exchange rows r and i
			{
				piv[i] = r;
                get_row(r, tmp_exch);
			    for (j=0; j<N; ++j)  
					Base::values_[r][j] = Base::values_[i][j];
                for (j=0; j<N; ++j)
				    Base::values_[i][j] = tmp_exch[j];
			}

			Base::values_[i][i] = 1 / Base::values_[i][i];
			for (j=0; j<M; ++j)
				if (j!=i)
					Base::values_[j][i] *= -Base::values_[i][i];
			for (j=0; j<M; ++j)
				for (k=0; k<N; ++k)
					if (j!=i && k!=i)
						Base::values_[j][k] += Base::values_[j][i] * Base::values_[i][k];
			for (j=0; j<N; ++j)
				if (j!=i)
					Base::values_[i][j] *= Base::values_[i][i];
		}

		for (i=N-1; i>=0; --i)
			if (i!=piv[i])
			{
                get_column(i, tmp_exch);
			    for (j=0; j<M; ++j)  
				    Base::values_[j][i] = Base::values_[j][piv[i]];
                for (j=0; j<M; ++j)
				    Base::values_[j][piv[i]] = tmp_exch[j];
			}

		return *this;
	}
	
	/**
	matrixT Solve (const matrixT& v) const _THROW_MATRIX_ERROR;
    matrixT Adj () _THROW_MATRIX_ERROR;
    T Cofact (size_t row, size_t col) _THROW_MATRIX_ERROR;
    T Cond () _NO_THROW;
    */

	//== type of matrices =====================================================

	bool is_square()
	{ return M==N; }

	bool is_null()
	{
#define expr(i,j)  if (Base::values_[i][j] != (Scalar)0.0) return false;
        unroll(expr);
#undef expr

		return true;
	}

	bool is_unit()
	{
#define expr(i,j)  if (i!=j) {if (Base::values_[i][j] != (Scalar)0.0) return false;} \
	else {if (Base::values_[i][j] != (Scalar)1.0) return false;}
        unroll(expr);
#undef expr

		return true;
	}

	bool is_symmetric()
	{
		if (M!=N)
			return false;

        int i, j;
		for (i = 1; i < M; ++i)
			for (j = 0; j < i; ++j)
				if (Base::values_[i][j] != Base::values_[j][i])
					return false;

		return true;
	}

	bool is_diagonal()
	{
        if (M!=N)
			return false;

#define expr(i,j)  if ((i != j) && (Base::values_[i][j] != Scalar(0))) return false;
        unroll(expr);
#undef expr

		return true;
	}

	// whether it is scalar: s*unit
	bool is_scalar()
	{
		if (!is_diagonal())
			return false;

		Scalar s = Base::values_[0][0];
		int i;
		for (i=1; i<M; ++i)
			if (Base::values_[i][i] != s)
				return false;

		return true;
	}

	bool is_upper_triangle()
	{
		if (M!=N)
			return false;

        int i, j;
		for (i = 1; i < M; ++i)
			for (j = 0; j < i; ++j)
				if (Base::values_[i][j] != Scalar(0))
					return false;

		return true;
	}

    bool is_lower_triangle()
	{
		if (M!=N)
			return false;

        int i, j;
		for (i = 0; i < M - 1; ++i)
			for (j = i + 1; j < N; ++j)
				if (Base::values_[i][j] != Scalar(0))
					return false;

		return true;
	}

	bool is_singular()
	{
        if (M!=N)
			return false;

		return (determinant() == Scalar(0));
	}
};



//== PARTIAL TEMPLATE SPECIALIZATIONS =========================================

//== FULL TEMPLATE SPECIALIZATIONS ============================================


//== GLOBAL FUNCTIONS =========================================================
// these could be friend function

// read the space-separated components of a matrix from a stream
template <typename Scalar, int M, int N>
inline std::istream&
operator>>(std::istream& _is, MatrixT<Scalar,M,N>& _m)
{
#define expr(i,j) _is >> _m(i,j);
    unroll(expr);
#undef expr

    return _is;
}

// output a matrix by printing its space-separated compontens
template <typename Scalar, int M, int N>
inline std::ostream&
operator<<(std::ostream& _os, const MatrixT<Scalar,M,N>& _m)
{
#define expr(i,j) _os << _m(i,j) << " ";
	unroll(expr);
#undef expr

    return _os;
}

// component-wise comparison, should be the same matrix type
template <typename Scalar, int M, int N>
inline bool 
operator==(const MatrixT<Scalar,M,N>& _m1, const MatrixT<Scalar,M,N>& _m2) 
{
#define expr(i,j) if(_m1(i,j) != _m2(i,j)) return false;
    unroll(expr);
#undef expr
    return true;
}

// component-wise comparison
template <typename Scalar, int M, int N>
inline bool 
operator!=(const MatrixT<Scalar,M,N>& _m1, const MatrixT<Scalar,M,N>& _m2)
{
    return !(_m1 == _m2);
}

// #undef  unroll  // \attention!

// component-wise matrix addition
template<typename Scalar,int M,int N>
inline MatrixT<Scalar,M,N> 
operator+(const MatrixT<Scalar,M,N>& _m1, const MatrixT<Scalar,M,N>& _m2)
{
    return MatrixT<Scalar,M,N>(_m1) += _m2;
}

// component-wise matrix difference
template<typename Scalar,int M,int N>
inline MatrixT<Scalar,M,N>
operator-(const MatrixT<Scalar,M,N>& _m1, const MatrixT<Scalar,M,N>& _m2)
{
    return MatrixT<Scalar,M,N>(_m1) -= _m2;
}

// scalar * matrix
template<typename Scalar,int M,int N>
inline MatrixT<Scalar,M,N> 
operator*(const Scalar& _s, const MatrixT<Scalar,M,N>& _m) 
{
    return MatrixT<Scalar,M,N>(_m) *= _s;
}

// matrix * scalar
template<typename Scalar,int M,int N>
inline MatrixT<Scalar,M,N> 
operator*(const MatrixT<Scalar,M,N>& _m, const Scalar& _s) 
{
    return MatrixT<Scalar,M,N>(_m) *= _s;
}

// matrix * matrix, general case, M*N, N*L -> M*L
// \attention: here we have used _m(k,j), so the parameter _m can not be const
template<typename Scalar,int M,int N,int L>
inline MatrixT<Scalar,M,L>
operator*(MatrixT<Scalar,M,N>& _m1, MatrixT<Scalar,N,L>& _m2 )
{
	MatrixT<Scalar,M,L> res_mat;

	int i,j,k;
	Scalar sum;
	for (i=0; i<M; ++i)
		for (j=0; j<L; ++j)
		{
            sum = (Scalar)0.0;
	    	for (k=0; k<N; ++k)
		    	sum += _m1(i,k) * _m2(k,j); 
			res_mat(i,j) = sum;
		}

	return res_mat;
}

// component-wise division by scalar
template<typename Scalar,int M,int N>
inline MatrixT<Scalar,M,N> 
operator/(const MatrixT<Scalar,M,N>& _m, const Scalar& _s)
{
    return MatrixT<Scalar,M,N>(_m) /= _s;
}

// binary power operator
template<typename Scalar,int M>
inline MatrixT<Scalar,M,M> 
operator^(const MatrixT<Scalar,M,M>& _m, const size_t& power)
{
   MatrixT<Scalar,M,M> res_mat = _m;
   res_mat ^= power;

   return res_mat;
}

// unary transpose operator
template<typename Scalar,int M,int N>
inline MatrixT<Scalar,N,M> 
operator~(const MatrixT<Scalar,M,N>& _m)
{
    MatrixT<Scalar,N,M> res_mat; // set null 
#define expr(i,j) res_mat(j,i) = _m(i,j);
    unroll(expr);
#undef expr

	return res_mat;
}

// unary inversion operator
template<typename Scalar,int M>
inline MatrixT<Scalar,M,M> 
operator!(const MatrixT<Scalar,M,M>& _m)
{
    MatrixT<Scalar,M,M> res_mat(_m);
	return res_mat.inversion();
}

// unary minus operator, seems not necessary
template<typename Scalar,int M,int N>
inline MatrixT<Scalar,M,N> 
operator-(const MatrixT<Scalar,M,N>& _m)
{
    MatrixT<Scalar,M,N> res_mat; // set null 
#define expr(i,j) res_mat(i,j) = -_m(i,j);
    unroll(expr);
#undef expr

	return res_mat;
}

#undef  unroll  // \attention

//== TYPEDEFS =================================================================

// 1 by 1 - byte signed matrix
typedef MatrixT<signed char,1,1> Mat1by1c;
// 1 by 1 - byte unsigned matrix
typedef MatrixT<unsigned char,1,1> Mat1by1uc;
// 1 by 1 - short signed matrix 
typedef MatrixT<signed short int,1,1> Mat1by1s;
// 1 by 1 - short unsigned matrix 
typedef MatrixT<unsigned short int,1,1> Mat1by1us;
// 1 by 1 - int signed matrix 
typedef MatrixT<signed int,1,1> Mat1by1i;
// 1 by 1 - int unsigned matrix 
typedef MatrixT<unsigned int,1,1> Mat1by1ui;
// 1 by 1 - float matrix 
typedef MatrixT<float,1,1> Mat1by1f;
// 1 by 1 - double matrix 
typedef MatrixT<double,1,1> Mat1by1d;

// 1 by 2 - byte signed matrix
typedef MatrixT<signed char,1,2> Mat1by2c;
// 1 by 2 - byte unsigned matrix
typedef MatrixT<unsigned char,1,2> Mat1by2uc;
// 1 by 2 - short signed matrix 
typedef MatrixT<signed short int,1,2> Mat1by2s;
// 1 by 2 - short unsigned matrix 
typedef MatrixT<unsigned short int,1,2> Mat1by2us;
// 1 by 2 - int signed matrix 
typedef MatrixT<signed int,1,2> Mat1by2i;
// 1 by 2 - int unsigned matrix 
typedef MatrixT<unsigned int,1,2> Mat1by2ui;
// 1 by 2 - float matrix 
typedef MatrixT<float,1,2> Mat1by2f;
// 1 by 2 - double matrix 
typedef MatrixT<double,1,2> Mat1by2d;

// 1 by 3 - byte signed matrix
typedef MatrixT<signed char,1,3> Mat1by3c;
// 1 by 3 - byte unsigned matrix
typedef MatrixT<unsigned char,1,3> Mat1by3uc;
// 1 by 3 - short signed matrix 
typedef MatrixT<signed short int,1,3> Mat1by3s;
// 1 by 3 - short unsigned matrix 
typedef MatrixT<unsigned short int,1,3> Mat1by3us;
// 1 by 3 - int signed matrix 
typedef MatrixT<signed int,1,3> Mat1by3i;
// 1 by 3 - int unsigned matrix 
typedef MatrixT<unsigned int,1,3> Mat1by3ui;
// 1 by 3 - float matrix 
typedef MatrixT<float,1,3> Mat1by3f;
// 1 by 3 - double matrix 
typedef MatrixT<double,1,3> Mat1by3d;

// 1 by 4 - byte signed matrix
typedef MatrixT<signed char,1,4> Mat1by4c;
// 1 by 4 - byte unsigned matrix
typedef MatrixT<unsigned char,1,4> Mat1by4uc;
// 1 by 4 - short signed matrix 
typedef MatrixT<signed short int,1,4> Mat1by4s;
// 1 by 4 - short unsigned matrix 
typedef MatrixT<unsigned short int,1,4> Mat1by4us;
// 1 by 4 - int signed matrix 
typedef MatrixT<signed int,1,4> Mat1by4i;
// 1 by 4 - int unsigned matrix 
typedef MatrixT<unsigned int,1,4> Mat1by4ui;
// 1 by 4 - float matrix 
typedef MatrixT<float,1,4> Mat1by4f;
// 1 by 4 - double matrix 
typedef MatrixT<double,1,4> Mat1by4d;

// 2 by 1 - byte signed matrix
typedef MatrixT<signed char,2,1> Mat2by1c;
// 2 by 1 - byte unsigned matrix
typedef MatrixT<unsigned char,2,1> Mat2by1uc;
// 2 by 1 - short signed matrix 
typedef MatrixT<signed short int,2,1> Mat2by1s;
// 2 by 1 - short unsigned matrix 
typedef MatrixT<unsigned short int,2,1> Mat2by1us;
// 2 by 1 - int signed matrix 
typedef MatrixT<signed int,2,1> Mat2by1i;
// 2 by 1 - int unsigned matrix 
typedef MatrixT<unsigned int,2,1> Mat2by1ui;
// 2 by 1 - float matrix 
typedef MatrixT<float,2,1> Mat2by1f;
// 2 by 1 - double matrix 
typedef MatrixT<double,2,1> Mat2by1d;

// 2 by 2 - byte signed matrix
typedef MatrixT<signed char,2,2> Mat2by2c;
// 2 by 2 - byte unsigned matrix
typedef MatrixT<unsigned char,2,2> Mat2by2uc;
// 2 by 2 - short signed matrix 
typedef MatrixT<signed short int,2,2> Mat2by2s;
// 2 by 2 - short unsigned matrix 
typedef MatrixT<unsigned short int,2,2> Mat2by2us;
// 2 by 2 - int signed matrix 
typedef MatrixT<signed int,2,2> Mat2by2i;
// 2 by 2 - int unsigned matrix 
typedef MatrixT<unsigned int,2,2> Mat2by2ui;
// 2 by 2 - float matrix 
typedef MatrixT<float,2,2> Mat2by2f;
// 2 by 2 - double matrix 
typedef MatrixT<double,2,2> Mat2by2d;

// 2 by 3 - byte signed matrix
typedef MatrixT<signed char,2,3> Mat2by3c;
// 2 by 3 - byte unsigned matrix
typedef MatrixT<unsigned char,2,3> Mat2by3uc;
// 2 by 3 - short signed matrix 
typedef MatrixT<signed short int,2,3> Mat2by3s;
// 2 by 3 - short unsigned matrix 
typedef MatrixT<unsigned short int,2,3> Mat2by3us;
// 2 by 3 - int signed matrix 
typedef MatrixT<signed int,2,3> Mat2by3i;
// 2 by 3 - int unsigned matrix 
typedef MatrixT<unsigned int,2,3> Mat2by3ui;
// 2 by 3 - float matrix 
typedef MatrixT<float,2,3> Mat2by3f;
// 2 by 3 - double matrix 
typedef MatrixT<double,2,3> Mat2by3d;

// 2 by 4 - byte signed matrix
typedef MatrixT<signed char,2,4> Mat2by4c;
// 2 by 4 - byte unsigned matrix
typedef MatrixT<unsigned char,2,4> Mat2by4uc;
// 2 by 4 - short signed matrix 
typedef MatrixT<signed short int,2,4> Mat2by4s;
// 2 by 4 - short unsigned matrix 
typedef MatrixT<unsigned short int,2,4> Mat2by4us;
// 2 by 4 - int signed matrix 
typedef MatrixT<signed int,2,4> Mat2by4i;
// 2 by 4 - int unsigned matrix 
typedef MatrixT<unsigned int,2,4> Mat2by4ui;
// 2 by 4 - float matrix 
typedef MatrixT<float,2,4> Mat2by4f;
// 2 by 4 - double matrix 
typedef MatrixT<double,2,4> Mat2by4d;

// 3 by 1 - byte signed matrix
typedef MatrixT<signed char,3,1> Mat3by1c;
// 3 by 1 - byte unsigned matrix
typedef MatrixT<unsigned char,3,1> Mat3by1uc;
// 3 by 1 - short signed matrix 
typedef MatrixT<signed short int,3,1> Mat3by1s;
// 3 by 1 - short unsigned matrix 
typedef MatrixT<unsigned short int,3,1> Mat3by1us;
// 3 by 1 - int signed matrix 
typedef MatrixT<signed int,3,1> Mat3by1i;
// 3 by 1 - int unsigned matrix 
typedef MatrixT<unsigned int,3,1> Mat3by1ui;
// 3 by 1 - float matrix 
typedef MatrixT<float,3,1> Mat3by1f;
// 3 by 1 - double matrix 
typedef MatrixT<double,3,1> Mat3by1d;

// 3 by 2 - byte signed matrix
typedef MatrixT<signed char,3,2> Mat3by2c;
// 3 by 2 - byte unsigned matrix
typedef MatrixT<unsigned char,3,2> Mat3by2uc;
// 3 by 2 - short signed matrix 
typedef MatrixT<signed short int,3,2> Mat3by2s;
// 3 by 2 - short unsigned matrix 
typedef MatrixT<unsigned short int,3,2> Mat3by2us;
// 3 by 2 - int signed matrix 
typedef MatrixT<signed int,3,2> Mat3by2i;
// 3 by 2 - int unsigned matrix 
typedef MatrixT<unsigned int,3,2> Mat3by2ui;
// 3 by 2 - float matrix 
typedef MatrixT<float,3,2> Mat3by2f;
// 3 by 2 - double matrix 
typedef MatrixT<double,3,2> Mat3by2d;

// 3 by 3 - byte signed matrix
typedef MatrixT<signed char,3,3> Mat3by3c;
// 3 by 3 - byte unsigned matrix
typedef MatrixT<unsigned char,3,3> Mat3by3uc;
// 3 by 3 - short signed matrix 
typedef MatrixT<signed short int,3,3> Mat3by3s;
// 3 by 3 - short unsigned matrix 
typedef MatrixT<unsigned short int,3,3> Mat3by3us;
// 3 by 3 - int signed matrix 
typedef MatrixT<signed int,3,3> Mat3by3i;
// 3 by 3 - int unsigned matrix 
typedef MatrixT<unsigned int,3,3> Mat3by3ui;
// 3 by 3 - float matrix 
typedef MatrixT<float,3,3> Mat3by3f;
// 3 by 3 - double matrix 
typedef MatrixT<double,3,3> Mat3by3d;

// 3 by 4 - byte signed matrix
typedef MatrixT<signed char,3,4> Mat3by4c;
// 3 by 4 - byte unsigned matrix
typedef MatrixT<unsigned char,3,4> Mat3by4uc;
// 3 by 4 - short signed matrix 
typedef MatrixT<signed short int,3,4> Mat3by4s;
// 3 by 4 - short unsigned matrix 
typedef MatrixT<unsigned short int,3,4> Mat3by4us;
// 3 by 4 - int signed matrix 
typedef MatrixT<signed int,3,4> Mat3by4i;
// 3 by 4 - int unsigned matrix 
typedef MatrixT<unsigned int,3,4> Mat3by4ui;
// 3 by 4 - float matrix 
typedef MatrixT<float,3,4> Mat3by4f;
// 3 by 4 - double matrix 
typedef MatrixT<double,3,4> Mat3by4d;

// 4 by 1 - byte signed matrix
typedef MatrixT<signed char,4,1> Mat4by1c;
// 4 by 1 - byte unsigned matrix
typedef MatrixT<unsigned char,4,1> Mat4by1uc;
// 4 by 1 - short signed matrix 
typedef MatrixT<signed short int,4,1> Mat4by1s;
// 4 by 1 - short unsigned matrix 
typedef MatrixT<unsigned short int,4,1> Mat4by1us;
// 4 by 1 - int signed matrix 
typedef MatrixT<signed int,4,1> Mat4by1i;
// 4 by 1 - int unsigned matrix 
typedef MatrixT<unsigned int,4,1> Mat4by1ui;
// 4 by 1 - float matrix 
typedef MatrixT<float,4,1> Mat4by1f;
// 4 by 1 - double matrix 
typedef MatrixT<double,4,1> Mat4by1d;

// 4 by 2 - byte signed matrix
typedef MatrixT<signed char,4,2> Mat4by2c;
// 4 by 2 - byte unsigned matrix
typedef MatrixT<unsigned char,4,2> Mat4by2uc;
// 4 by 2 - short signed matrix 
typedef MatrixT<signed short int,4,2> Mat4by2s;
// 4 by 2 - short unsigned matrix 
typedef MatrixT<unsigned short int,4,2> Mat4by2us;
// 4 by 2 - int signed matrix 
typedef MatrixT<signed int,4,2> Mat4by2i;
// 4 by 2 - int unsigned matrix 
typedef MatrixT<unsigned int,4,2> Mat4by2ui;
// 4 by 2 - float matrix 
typedef MatrixT<float,4,2> Mat4by2f;
// 4 by 2 - double matrix 
typedef MatrixT<double,4,2> Mat4by2d;

// 4 by 3 - byte signed matrix
typedef MatrixT<signed char,4,3> Mat4by3c;
// 4 by 3 - byte unsigned matrix
typedef MatrixT<unsigned char,4,3> Mat4by3uc;
// 4 by 3 - short signed matrix 
typedef MatrixT<signed short int,4,3> Mat4by3s;
// 4 by 3 - short unsigned matrix 
typedef MatrixT<unsigned short int,4,3> Mat4by3us;
// 4 by 3 - int signed matrix 
typedef MatrixT<signed int,4,3> Mat4by3i;
// 4 by 3 - int unsigned matrix 
typedef MatrixT<unsigned int,4,3> Mat4by3ui;
// 4 by 3 - float matrix 
typedef MatrixT<float,4,3> Mat4by3f;
// 4 by 3 - double matrix 
typedef MatrixT<double,4,3> Mat4by3d;

// 4 by 4 - byte signed matrix
typedef MatrixT<signed char,4,4> Mat4by4c;
// 4 by 4 - byte unsigned matrix
typedef MatrixT<unsigned char,4,4> Mat4by4uc;
// 4 by 4 - short signed matrix 
typedef MatrixT<signed short int,4,4> Mat4by4s;
// 4 by 4 - short unsigned matrix 
typedef MatrixT<unsigned short int,4,4> Mat4by4us;
// 4 by 4 - int signed matrix 
typedef MatrixT<signed int,4,4> Mat4by4i;
// 4 by 4 - int unsigned matrix 
typedef MatrixT<unsigned int,4,4> Mat4by4ui;
// 4 by 4 - float matrix 
typedef MatrixT<float,4,4> Mat4by4f;
// 4 by 4 - double matrix 
typedef MatrixT<double,4,4> Mat4by4d;

#endif  // MATRIXT_H defined 