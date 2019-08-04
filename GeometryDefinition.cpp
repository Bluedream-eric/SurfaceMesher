#include "GeometryDefinition.h"

void FergusonCurve::calculate_control_point_tangent_vector(EndCondition _ec, 
		Vec3d _v0, Vec3d _vn)
{
	int ncp = get_control_point_number();
	assert(ncp>=2);
	control_point_tangent_vector_.assign(ncp, Vec3d(0.0,0.0,0.0)); // its size() should be ncp
    
	Vec3d vec(0.0);
	double piv = 0.0;
	double *aux = new double [ncp];
	int i;
	for (i=0; i<ncp; ++i)
		aux[i] = 0.0;

	if (2==ncp)
	{
		control_point_tangent_vector_[0] = control_point_[1] - control_point_[0];
		control_point_tangent_vector_[1] = control_point_tangent_vector_[0];
	}
	else if (3==ncp)
	{
        if (SpecifiedTangent==_ec)
		{
            control_point_tangent_vector_[0] = _v0;
			control_point_tangent_vector_[2] = _vn;
			control_point_tangent_vector_[1] = 0.25*(3.0*(control_point_[2]-control_point_[0]) - _v0 - _vn);
		}
		else if (ZeroSecondDerivative==_ec)
		{
            control_point_tangent_vector_[0] = -1.25*control_point_[0]
			    + 1.50*control_point_[1] - 0.25*control_point_[2];
            control_point_tangent_vector_[1] = -0.5*control_point_[0] + 0.5*control_point_[2];
			control_point_tangent_vector_[2] = 0.25*control_point_[0]
			    - 1.50*control_point_[1] + 1.25*control_point_[2];
		}
		else
		{
			// todo_clg
			assert(false);
		}
	}
	else // ncp>3, gaussian elimination with backsubstitution
	{
		// attention: AX=B, then X=inverse(A)*B, here the matrix need to be dynamic,
		// thus the existing Matrix<M,N> class doesn't work
		// this is a tri-diagonal matrix, so it can be solved specially and don't use common methods

		// the first two rows
        vec = 3.0*(control_point_[2] - control_point_[0]);
        if (SpecifiedTangent==_ec)
		{
			control_point_tangent_vector_[0] = _v0;
			control_point_tangent_vector_[ncp-1] = _vn;
			aux[0] = 0.0;

			vec = vec - _v0;
			piv = 4.0;
		}
		else if (ZeroSecondDerivative==_ec)
		{
			control_point_tangent_vector_[0] = 1.5*(control_point_[1] - control_point_[0]);
			aux[0] = 0.5;

			vec = vec - 1.5*(control_point_[1] - control_point_[0]);
			piv = 4.0 - 0.5;
		}
		else
		{
			// todo_clg
			assert(false);
		}
        control_point_tangent_vector_[1] = vec / piv;

		// middle rows of the type [1, 4, 1]
		for (i=2; i<=ncp-3; ++i)
		{
			aux[i-1] = 1.0 / piv;

			piv = 4.0 - aux[i-1];
			assert(!is_zero(fabs(piv), 1.0e-12));
            vec = 3.0*(control_point_[i+1] - control_point_[i-1]);
			control_point_tangent_vector_[i] = (vec - control_point_tangent_vector_[i-1]) / piv;
		}
        aux[ncp-3] = 1.0 / piv;

		// the second last row i=ncp-2
		// first elimination with the second last row i=ncp-1
		vec = 3.0*(control_point_[ncp-1] - control_point_[ncp-3]);
		if (SpecifiedTangent==_ec)
		{
			vec = vec - _vn; 
			piv = 4.0;
		}
		else if (ZeroSecondDerivative==_ec)
		{
			vec = vec - 1.5*(control_point_[ncp-1] - control_point_[ncp-2]);
			piv = 4.0 - 0.5;
		}
		else
		{
			// todo_clg
			assert(false);
		}

		// second elimination with the above row to get the final result!
		piv = piv - aux[ncp-3];
        assert(!is_zero(fabs(piv), 1.0e-12));
		control_point_tangent_vector_[ncp-2] = (vec - control_point_tangent_vector_[ncp-3]) / piv;

		// backsubstitution
		for (i=ncp-3; i>=1; --i)
			control_point_tangent_vector_[i] -= aux[i]*control_point_tangent_vector_[i+1];
		
		// end values when _ec=ZeroSecondDerivative
		if (ZeroSecondDerivative==_ec)
		{
			control_point_tangent_vector_[0] = 1.5*(control_point_[1]-control_point_[0])
				- 0.5*control_point_tangent_vector_[1];
            control_point_tangent_vector_[ncp-1] = 1.5*(control_point_[ncp-1]-control_point_[ncp-2])
				- 0.5*control_point_tangent_vector_[ncp-2];
		}
	}

	delete [] aux; // never forget

	return;
}

void FergusonCurve::calculate_segment_coefficient()
{
	int ncp = get_control_point_number();
	int i;
	Vec3d a1, a2, a3;

	segment_coefficient_.assign(ncp-1, Vec5d(0.0)); // warning_clg
	for (i=0; i<ncp-1; ++i)
	{
		a1 = control_point_tangent_vector_[i];
		a2 = 3.0*(control_point_.at(i+1) - control_point_[i]) - 2.0*control_point_tangent_vector_[i]
		    - control_point_tangent_vector_[i+1];
		a3 = 2.0*(control_point_[i] - control_point_[i+1]) + control_point_tangent_vector_[i]
		    + control_point_tangent_vector_[i+1];

		segment_coefficient_.at(i)[0] = a1 | a1;
        segment_coefficient_.at(i)[1] = 4*(a1 | a2);
		segment_coefficient_.at(i)[2] = 4*(a2 | a2) + 6*(a1 | a3);
		segment_coefficient_.at(i)[3] = 12*(a2 | a3);
		segment_coefficient_.at(i)[4] = 9*(a3 | a3);
	}

	return;
}

void FergusonCurve::calculate_partial_segment_length(const Vec5d& _coef, double _u0, 
        double _u1, double& _len, double _tol)
{
	// assert(0.0<=_u0 && _u0<_u1 && _u1<=1.0);

    double f0, f1;
	f0 = _coef[0]+_u0*(_coef[1]+_u0*(_coef[2]+_u0*(_coef[3]+_u0*_coef[4])));
	f1 = _coef[0]+_u1*(_coef[1]+_u1*(_coef[2]+_u1*(_coef[3]+_u1*_coef[4])));
	assert(f0>=0.0 && f1>=0.0);

	if (f0<_tol*_tol || f1<_tol*_tol) 
	{
        f0 = sqrt(f0);
		f1 = sqrt(f1);
        _len = 0.5*(f0+f1)*(_u1-_u0);
		return;
	}

    f0 = sqrt(f0);
	f1 = sqrt(f1);

	const int nit = 20;
	const double u01=_u1-_u0;
	double eps, del, x, f, new_sum; 
	double orig_sum=f0+f1, orig_len=-1.0e30;
	int i, j, nn;

	// use iterative composite simpson formula to calculate the arc length integral
	for (i=1, nn=1; i<=nit; ++i, nn*=2)
	{
        del = u01 / nn;
        new_sum = 0.0;
		x = _u0 + 0.5*del;
		for (j=1; j<=nn; ++j) // this loop can be put into one "for" statement
		{
			f = _coef[0] + x*(_coef[1] + x*(_coef[2] + x*(_coef[3] + x*_coef[4])));
			new_sum += sqrt(f);
			x += del;
		}
		_len = (orig_sum + 4*new_sum) * u01 / (6.0*nn);

		eps = _tol * fabs(orig_len);
		if (fabs(_len-orig_len) < eps)
			return;

		orig_sum += 2 * new_sum;
		orig_len = _len;
	}

	cout << "calculate_partial_segment_length: number of iterations exceeded" << endl;
	assert(false);

	return;
}

void FergusonCurve::calculate_segment_length(double _tol)
{
    int ncp = get_control_point_number();
	double len;
	int i;
	segment_length_.assign(ncp-1, 0.0);
    
	for (i=0; i<ncp-1; ++i)
	{
        calculate_partial_segment_length(segment_coefficient_[i], 0.0, 1.0, len, _tol);
        segment_length_[i] = len;
	}

    for (i=1; i<ncp-1; ++i) // accumulating
        segment_length_[i] += segment_length_[i-1];

	return;
}

    
void FergusonCurve::calculate_parameter_given_arc_length(const Vec5d& _coef, double _seg_len, double _u0, 
		double& _u1, double _tar_len, double _tol)
{
    // this is calculated by means of a Newton-Raphson iteration
	const int nit = 20;
	double len, len_err, f, tol2 = 1.0e-4;
	int i;

	_u1 = _tar_len / _seg_len;
	assert(0.0<=_u1 && _u1<=1.0);
	if (_u1<_tol)
		_u1 = 0.0;
	else if (1.0-_u1<_tol)
		_u1 = 1.0;
	else
	{
		for (i=0; i<nit; ++i)
		{
			_u1 = _u1<1.0 ? _u1 : 1.0;
            calculate_partial_segment_length(_coef, _u0, _u1, len, _tol);
			len_err = _tar_len - len;
            if (fabs(len_err/_tar_len)<tol2)
				return;
            f = sqrt(_coef[0]+_u1*(_coef[1]+_u1*(_coef[2]+_u1*(_coef[3]+_u1*_coef[4]))));
			_u1 += len_err / f;
		}
		assert(false); // shouldn't get here
	}

	return;
}

void FergusonCurve::evaluate_ferguson_curve_segment(DerivativeNumber _dn, Vec3d _v0[2], 
		Vec3d _v1[2], double _u, Vec3d _res[3]) const
{
    Vec3d a[4];
	a[0] = _v0[0];
	a[1] = _v0[1];
	a[2] = 3.0*(_v1[0] - _v0[0]) - 2.0*_v0[1] - _v1[1];
	a[3] = 2.0*(_v0[0] - _v1[0]) + _v0[1] + _v1[1];

	_res[0] = a[0]+_u*(a[1]+_u*(a[2]+_u*a[3]));
	if (Zero!=_dn)
	{
        _res[1] = a[1]+_u*(2.0*a[2]+_u*3.0*a[3]);
		if (Two==_dn)
			_res[2] = 2.0*a[2]+_u*6.0*a[3];
	}

	return;
}

void FergusonCurve::initialize()
{
	calculate_control_point_tangent_vector(ZeroSecondDerivative, 
		Vec3d(0.0,0.0,0.0), Vec3d(0.0,0.0,0.0));
	calculate_segment_coefficient();
	calculate_segment_length(1.e-5);

	return;
}

void FergusonCurve::evaluate(double _t, Vec3d& _c) const
{
	int ncp = get_control_point_number();
	assert(_t>=0 && _t<=ncp-1);
	int seg_num = (int)_t;
	double u = _t - seg_num;
	Vec3d a0, a1, a2, a3;

	a0 = control_point_[seg_num];
	a1 = control_point_tangent_vector_[seg_num];
	a2 = 3.0*(control_point_[seg_num+1] - control_point_[seg_num]) - 2.0*control_point_tangent_vector_[seg_num]
	    - control_point_tangent_vector_[seg_num+1];
	a3 = 2.0*(control_point_[seg_num] - control_point_[seg_num+1]) + control_point_tangent_vector_[seg_num]
	    + control_point_tangent_vector_[seg_num+1];

    _c = a0 + u*(a1 + u*(a2 + u*a3));

	return;
}

void FergusonCurve::calculate_first_derivative(double _t, Vec3d& _ct) const
{
	int ncp = get_control_point_number();
	assert(_t>=0 && _t<=ncp-1);
	int seg_num = (int)_t;
	double u = _t - seg_num;
	Vec3d a1, a2, a3;

	a1 = control_point_tangent_vector_[seg_num];
	a2 = 3.0*(control_point_[seg_num+1] - control_point_[seg_num]) - 2.0*control_point_tangent_vector_[seg_num]
	    - control_point_tangent_vector_[seg_num+1];
	a3 = 2.0*(control_point_[seg_num] - control_point_[seg_num+1]) + control_point_tangent_vector_[seg_num]
	    + control_point_tangent_vector_[seg_num+1];

    _ct = a1 + u*(2.0*a2 + u*3*a3);

	return;
}

void FergusonCurve::calculate_second_derivative(double _t, Vec3d& _ctt) const
{
	int ncp = get_control_point_number();
	assert(_t>=0 && _t<=ncp-1);
	int seg_num = (int)_t;
	double u = _t - seg_num;
	Vec3d a2, a3;

	a2 = 3.0*(control_point_[seg_num+1] - control_point_[seg_num]) - 2.0*control_point_tangent_vector_[seg_num]
	    - control_point_tangent_vector_[seg_num+1];
	a3 = 2.0*(control_point_[seg_num] - control_point_[seg_num+1]) + control_point_tangent_vector_[seg_num]
	    + control_point_tangent_vector_[seg_num+1];

    _ctt = 2.0*a2 + 6.0*u*a3;

	return;
}

void FergusonCurve::calculate_projection(Vec3d _c, double& _t) const
{
	// todo_clg
}

void FergusonCurve::calculate_arc_length(double _t1, double _t2, double& _length) const
{
	// todo_clg
	// may calculate segment by segment
}


void FergusonSurface::calculate_ferguson_curve_control_point_tangent_vector(std::vector<Vec3d>& _control_point,
		std::vector<Vec3d>& _tangent_vector, EndCondition _ec, Vec3d _v0, Vec3d _vn)
{
	int ncp = (int)_control_point.size();
	assert(ncp>=2);
	_tangent_vector.assign(ncp, Vec3d(0.0,0.0,0.0)); // its size() should be ncp
    
	Vec3d vec(0.0);
	double piv = 0.0;
	double *aux = new double [ncp];
	int i;
	for (i=0; i<ncp; ++i)
		aux[i] = 0.0;

	if (2==ncp)
	{
		_tangent_vector[0] = _control_point[1] - _control_point[0];
		_tangent_vector[1] = _tangent_vector[0];
	}
	else if (3==ncp)
	{
        if (SpecifiedTangent==_ec)
		{
            _tangent_vector[0] = _v0;
			_tangent_vector[2] = _vn;
			_tangent_vector[1] = 0.25*(3.0*(_control_point[2]-_control_point[0]) - _v0 - _vn);
		}
		else if (ZeroSecondDerivative==_ec)
		{
            _tangent_vector[0] = -1.25*_control_point[0]
			    + 1.50*_control_point[1] - 0.25*_control_point[2];
            _tangent_vector[1] = -0.5*_control_point[0] + 0.5*_control_point[2];
			_tangent_vector[2] = 0.25*_control_point[0]
			    - 1.50*_control_point[1] + 1.25*_control_point[2];
		}
		else
		{
			// todo_clg
			assert(false);
		}
	}
	else // ncp>3, gaussian elimination with backsubstitution
	{
		// attention: AX=B, then X=inverse(A)*B, here the matrix need to be dynamic,
		// thus the existing Matrix<M,N> class doesn't work
		// this is a tri-diagonal matrix, so it can be solved specially and don't use common methods

		// the first two rows
        vec = 3.0*(_control_point[2] - _control_point[0]);
        if (SpecifiedTangent==_ec)
		{
			_tangent_vector[0] = _v0;
			_tangent_vector[ncp-1] = _vn;
			aux[0] = 0.0;

			vec = vec - _v0;
			piv = 4.0;
		}
		else if (ZeroSecondDerivative==_ec)
		{
			_tangent_vector[0] = 1.5*(_control_point[1] - _control_point[0]);
			aux[0] = 0.5;

			vec = vec - 1.5*(_control_point[1] - _control_point[0]);
			piv = 4.0 - 0.5;
		}
		else
		{
			// todo_clg
			assert(false);
		}
        _tangent_vector[1] = vec / piv;

		// middle rows of the type [1, 4, 1]
		for (i=2; i<=ncp-3; ++i)
		{
			aux[i-1] = 1.0 / piv;

			piv = 4.0 - aux[i-1];
			assert(!is_zero(fabs(piv), 1.0e-12));
            vec = 3.0*(_control_point[i+1] - _control_point[i-1]);
			_tangent_vector[i] = (vec - _tangent_vector[i-1]) / piv;
		}
        aux[ncp-3] = 1.0 / piv;

		// the second last row i=ncp-2
		// first elimination with the second last row i=ncp-1
		vec = 3.0*(_control_point[ncp-1] - _control_point[ncp-3]);
		if (SpecifiedTangent==_ec)
		{
			vec = vec - _vn; 
			piv = 4.0;
		}
		else if (ZeroSecondDerivative==_ec)
		{
			vec = vec - 1.5*(_control_point[ncp-1]-_control_point[ncp-2]); // originally is vec - _vn; wrong!
			piv = 4.0 - 0.5;
		}
		else
		{
			// todo_clg
			assert(false);
		}

		// second elimination with the above row to get the final result!
		piv = piv - aux[ncp-3];
        assert(!is_zero(fabs(piv), 1.0e-12));
		_tangent_vector[ncp-2] = (vec - _tangent_vector[ncp-3]) / piv;

		// backsubstitution
		for (i=ncp-3; i>=1; --i)
			_tangent_vector[i] -= aux[i]*_tangent_vector[i+1];
		
		// end values when _ec=ZeroSecondDerivative
		if (ZeroSecondDerivative==_ec)
		{
			_tangent_vector[0] = 1.5*(_control_point[1]-_control_point[0])
				- 0.5*_tangent_vector[1];
            _tangent_vector[ncp-1] = 1.5*(_control_point[ncp-1]-_control_point[ncp-2])
				- 0.5*_tangent_vector[ncp-2];
		}
	}

	delete [] aux; // never forget

	return;
}
void FergusonSurface::calculate_control_point_derivative()
{
	std::vector<Vec3d> control_point, tangent_vector;
	int i, j, k;
	Vec3d t0, t1;

	// allocate memory
	ru_.assign(nu_*nv_, Vec3d(0.0));
	rv_.assign(nu_*nv_, Vec3d(0.0));
	ruv_.assign(nu_*nv_, Vec3d(0.0));

    // interpolate ru in the u direction
    for (i=0; i<nv_; ++i)
	{
		for (j=i*nu_; j<(i+1)*nu_; ++j)
		    control_point.push_back(control_point_[j]);		
        calculate_ferguson_curve_control_point_tangent_vector(control_point, tangent_vector);
        for (j=0; j<nu_; ++j)
		    ru_[i*nu_+j] = tangent_vector[j];

		control_point.clear();
		tangent_vector.clear();
	}

    // interpolate rv in the v direction
    for (i=0; i<nu_; ++i)
	{
        for (j=0; j<nv_; ++j)
		{
            k = j*nu_ + i;
			control_point.push_back(control_point_[k]);	
		}	
        calculate_ferguson_curve_control_point_tangent_vector(control_point, tangent_vector);
		for (j=0; j<nv_; ++j)
		{
            k = j*nu_ + i;
			rv_[k] = tangent_vector[j];
		}

		control_point.clear();
		tangent_vector.clear();
	}

	// interpolate ruv

	/** first interpolate to get the mixed partial derivative vectors for the 
	* first and the last v-direction curve, and then use these vectors as the
	* specified tangent vector to get the mixed partial derivative vectors for 
	* all the u-direction curve, i.e. for all the control points.
	*/
    
    // interpolate in the v direction for the first and the last curves
	for (i=0; i<nv_; ++i)
		control_point.push_back(ru_[i*nu_]);
    calculate_ferguson_curve_control_point_tangent_vector(control_point, tangent_vector);
    for (i=0; i<nv_; ++i)
		ruv_[i*nu_] = tangent_vector[i];
	control_point.clear();
	tangent_vector.clear();

    for (i=0; i<nv_; ++i)
		control_point.push_back(ru_[(i+1)*nu_-1]);
    calculate_ferguson_curve_control_point_tangent_vector(control_point, tangent_vector);
    for (i=0; i<nv_; ++i)
		ruv_[(i+1)*nu_-1] = tangent_vector[i];
	control_point.clear();
	tangent_vector.clear();

	// interpolate for all the u direction curves
    for (i=0; i<nv_; ++i)
	{
		for (j=0; j<nu_; ++j)
		{
			k = i*nu_ + j;
			control_point.push_back(rv_[k]);
		}
        t0 = ruv_[i*nu_];
		t1 = ruv_[(i+1)*nu_-1];
        calculate_ferguson_curve_control_point_tangent_vector(control_point,
			tangent_vector, SpecifiedTangent, t0, t1);
        for (j=0; j<nu_; ++j)
		{
			k = i*nu_ + j;
			ruv_[k] = tangent_vector[k];
		}
	}

	return;
}

void FergusonSurface::initialize()
{
	calculate_control_point_derivative();

	return;
}

void FergusonSurface::evaluate(Vec2d _uv, Vec3d& _s) const
{
    const int nu = get_nu(), nv = get_nv();

	// for debug
	// cout << "In function evaluate(): " << "nu: " << nu << " nv: " << nv << endl;
	// cout << "_uv: " << _uv[0] << " " << _uv[1] << endl;

	assert(0.0<=_uv[0] && _uv[0]<=nu-1 && 0.0<=_uv[1] && _uv[1]<=nv-1);
	int i;

    Vec2i iuv = Vec2i(int(_uv[0]), int(_uv[1]));
	if (nu-1==iuv[0]) 
        iuv[0] -= 1;
	if (nv-1==iuv[1])
		iuv[1] -= 1;
	Vec2d uv = Vec2d(_uv[0]-iuv[0], _uv[1]-iuv[1]);

	// the control point information of the containing ferguson patch
	Vec3d r[4], ru[4], rv[4], ruv[4]; 
	int cp_index[4];
	cp_index[0] = nu*iuv[1] + iuv[0];
	cp_index[1] = cp_index[0] + 1;
	cp_index[2] = cp_index[0] + nu;
	cp_index[3] = cp_index[2] + 1;
	for (i=0; i<4; ++i)
	{
        r[i]   = control_point_[cp_index[i]];
		ru[i]  = ru_[cp_index[i]];
		rv[i]  = rv_[cp_index[i]];
		ruv[i] = ruv_[cp_index[i]];
	}

	Mat1by4d mat_u = ~Mat4by1d(1, uv[0], uv[0]*uv[0], uv[0]*uv[0]*uv[0]);
	Mat4by1d mat_v(1, uv[1], uv[1]*uv[1], uv[1]*uv[1]*uv[1]);
	Mat4by4d a, m(1, 0, 0, 0, 0, 0, 1, 0, -3, 3, -2, -1, 2, -2, 1, 1);
	Mat1by1d res;
    
    for (i=0; i<3; ++i)
	{
		// a = Mat4by4d(r[0][i], r[1][i], rv[0][i], rv[1][i], r[2][i], r[3][i], rv[2][i], rv[3][i], 
			// ru[0][i], ru[1][i], ruv[0][i], ruv[1][i], ru[2][i], ru[3][i], ruv[2][i], ruv[3][i]);
		a = Mat4by4d(r[0][i], r[2][i], rv[0][i], rv[2][i], r[1][i], r[3][i], rv[1][i], rv[3][i], 
			ru[0][i], ru[2][i], ruv[0][i], ruv[2][i], ru[1][i], ru[3][i], ruv[1][i], ruv[3][i]);
		res = mat_u * m * a * (~m) * mat_v;

		_s[i] = res(0, 0);
	}

	return;
}

// need to be incorporated in evaluate() not only for efficiency, but also for convenience
void FergusonSurface::calculate_first_derivative(Vec2d _uv, Vec3d& _su, Vec3d& _sv) const
{
	const int nu = get_nu(), nv = get_nv();
	assert(0.0<=_uv[0] && _uv[0]<=nu-1 && 0.0<=_uv[1] && _uv[1]<=nv-1);

    Vec2i iuv = Vec2i(int(_uv[0]), int(_uv[1]));
	if (nu-1==iuv[0]) 
        iuv[0] -= 1;
	if (nv-1==iuv[1])
		iuv[1] -= 1;
	Vec2d uv = Vec2d(_uv[0]-iuv[0], _uv[1]-iuv[1]);

	// the control point information of the containing ferguson patch
	Vec3d r[4], ru[4], rv[4], ruv[4]; 
	int i;
	int cp_index[4];
	cp_index[0] = nu*iuv[1] + iuv[0];
	cp_index[1] = cp_index[0] + 1;
	cp_index[2] = cp_index[0] + nu;
	cp_index[3] = cp_index[2] + 1;
	for (i=0; i<4; ++i)
	{
        r[i]   = control_point_[cp_index[i]];
		ru[i]  = ru_[cp_index[i]];
		rv[i]  = rv_[cp_index[i]];
		ruv[i] = ruv_[cp_index[i]];
	}

	
	Mat4by4d a, m(1, 0, 0, 0, 0, 0, 1, 0, -3, 3, -2, -1, 2, -2, 1, 1);
	Mat1by1d res;
    
	Mat1by4d mat_u = ~Mat4by1d(0., 1, 2*uv[0], 3*uv[0]*uv[0]);
	Mat4by1d mat_v(1, uv[1], uv[1]*uv[1], uv[1]*uv[1]*uv[1]);
    for (i=0; i<3; ++i)
	{
		// a = Mat4by4d(r[0][i], r[1][i], rv[0][i], rv[1][i], r[2][i], r[3][i], rv[2][i], rv[3][i], 
			// ru[0][i], ru[1][i], ruv[0][i], ruv[1][i], ru[2][i], ru[3][i], ruv[2][i], ruv[3][i]);
		a = Mat4by4d(r[0][i], r[2][i], rv[0][i], rv[2][i], r[1][i], r[3][i], rv[1][i], rv[3][i], 
			ru[0][i], ru[2][i], ruv[0][i], ruv[2][i], ru[1][i], ru[3][i], ruv[1][i], ruv[3][i]);
		res = mat_u * m * a * (~m) * mat_v;

		_su[i] = res(0, 0);
	}

    mat_u = ~Mat4by1d(1, uv[0], uv[0]*uv[0], uv[0]*uv[0]*uv[0]);
	mat_v = Mat4by1d(0, 1, 2*uv[1], 3*uv[1]*uv[1]);
	for (i=0; i<3; ++i)
	{
		// a = Mat4by4d(r[0][i], r[1][i], rv[0][i], rv[1][i], r[2][i], r[3][i], rv[2][i], rv[3][i], 
			// ru[0][i], ru[1][i], ruv[0][i], ruv[1][i], ru[2][i], ru[3][i], ruv[2][i], ruv[3][i]);
		a = Mat4by4d(r[0][i], r[2][i], rv[0][i], rv[2][i], r[1][i], r[3][i], rv[1][i], rv[3][i], 
			ru[0][i], ru[2][i], ruv[0][i], ruv[2][i], ru[1][i], ru[3][i], ruv[1][i], ruv[3][i]);
		res = mat_u * m * a * (~m) * mat_v;

		_sv[i] = res(0, 0);
	}

	return;
}

// need to be incorporated in evaluate() not only for efficiency, but also for convenience
void FergusonSurface::calculate_second_derivative(const Vec2d& _uv, Vec3d& _suu, Vec3d& _suv, 
		Vec3d& _svv) const
{

	return;
}

// need to be incorporated in evaluate() not only for efficiency, but also for convenience
void FergusonSurface::calculate_first_and_second_derivative(Vec2d _uv, Vec3d& _su, 
		Vec3d& _sv, Vec3d& _suu, Vec3d& _suv, Vec3d& _svv) const
{

	return;
}

void FergusonSurface::calculate_pricipal_curvature_and_direction(Vec2d _uv, 
        double& _kmin, double& _kmax, Vec3d& _dmin, Vec3d& _dmax) const
{

	return;
}

void FergusonSurface::calculate_projection(Vec3d _s, Vec2d& _uv) const
{

	return;
}


bool Geometry::read(const std::string& _filename)
{
	std::ifstream ifs(_filename.c_str());
	assert(ifs);

	std::string extname = get_file_extension(_filename);
	bool res = false;
	if ("stl"==extname || "STL"==extname)
	{
		res = read_from_stl(ifs);
		model_type_ = "stl";
	}
	else if ("fli"==extname || "FLI"==extname)
	{
		res = read_from_fli(ifs);
		model_type_ = "fli";
	}

	// todo_clg: may need to add other geometric file formats
	
	return res;
}

std::string Geometry::get_file_extension(const std::string& _filename) const
{
	std::string extname;
	std::string::size_type idx = _filename.rfind('.');
	assert(std::string::npos != idx);
	extname = _filename.substr(idx+1);

	return extname;
}

bool Geometry::read_from_stl(std::ifstream& _ifs)
{
	int i, j, k;
	Vec3d pt;
	std::vector<Vec3d> vec_pt;
	char cdummy[256];
	_ifs.getline(cdummy, 256);

	// bbox_.min() = (Vec3d)(1.e38, 1.e38, 1.e38);
	// bbox_.max() = (Vec3d)(-1.e38, -1.e38, -1.e38);

	while (1)
	{
		_ifs.getline(cdummy, 256);
		if (strncmp(cdummy, "end", 3) == 0) // end of file
			break;
		_ifs.getline(cdummy, 256);
		for (i=0; i<3; ++i)
		{
		    _ifs >> cdummy >> pt[0] >> pt[1] >> pt[2];
			vec_pt.push_back(pt);
		}
		_ifs.getline(cdummy, 256);
		_ifs.getline(cdummy, 256);
		_ifs.getline(cdummy, 256);
	}
	_ifs.close();

	// fill in the class Geometry
	bool existed = false;
	for (i=0; i<(int)vec_pt.size()/3; ++i)
	{
		Surface* psf = new LinearTriangleSurface(vec_pt[3*i], vec_pt[3*i+1], vec_pt[3*i+2]);
		surface_geometry_.push_back(psf);

		int* pcl = new int[3];  // curve loop
		for (j=0; j<=2; ++j)
		{
			existed = false;
			for (k=0; k<(int)curve_geometry_.size(); ++k)  // search to see if existed
				if ((((LineCurve*)curve_geometry_[k])->get_first_point() == vec_pt[3*i+j]  // use tolerance may be better
				    && ((LineCurve*)curve_geometry_[k])->get_second_point() == vec_pt[3*i+(j+1)%3])
					|| (((LineCurve*)curve_geometry_[k])->get_second_point() == vec_pt[3*i+j] 
				    && ((LineCurve*)curve_geometry_[k])->get_first_point() == vec_pt[3*i+(j+1)%3]))
				{
                    existed = true;
					break;
				}

			if (existed)
				pcl[j] = k;
			else
			{
				Curve* pcv = new LineCurve(vec_pt[3*i+j], vec_pt[3*i+(j+1)%3]);
				curve_geometry_.push_back(pcv);
				pcl[j] = (int)curve_geometry_.size() - 1;
			}
		}

        SurfaceRegion* psr = new SurfaceRegion((int)surface_geometry_.size(), 3, pcl); // todo_clg: originally "(int)surface_geometry_.size()-1" 20071231
		surface_region_.push_back(psr); // pass value
	}

	for (i=0; i<(int)curve_geometry_.size(); ++i)
		curve_segment_.push_back(i); // seems not needed

	return true;
}

bool Geometry::read_from_fli(std::ifstream& _ifs)
{
    int ncv, nsf, ncs, nsr;
	int n, nu, nv, surface_id, n_curve;
	Curve *pcv;
	Surface *psf;
	SurfaceRegion *psr;
	int *curve_loop;
	int i,j;
	int idummy;
    Vec3d pt;
	std::vector<Vec3d> vec_pt;
	char cdummy[256];

	_ifs.getline(cdummy, 256); // "1.- Geometry Definition" label
	_ifs >> ncv >> nsf;
	_ifs.getline(cdummy, 256);
	_ifs.getline(cdummy, 256); // "Curves" label

	// the curve number/index starts with 0 here, while in the fli file it's 1
	for (i=0; i<ncv; ++i)
	{
		_ifs.getline(cdummy, 256); // curve number, intersection type: always = 1
		_ifs >> n; // number of control points defining the ferguson curve
		_ifs.getline(cdummy, 256);

		for (j=0; j<n; ++j)
		{
            _ifs >> pt[0] >> pt[1] >> pt[2];
			vec_pt.push_back(pt);
			_ifs.getline(cdummy, 256);
		}

		pcv = new FergusonCurve(vec_pt); // value passed
		vec_pt.clear();
		curve_geometry_.push_back(pcv);
	}

	_ifs.getline(cdummy, 256); // "Support_Surfaces" label
    // the curve number/index starts with 0 here, while in the fli file it's 1
	for (i=0; i<nsf; ++i)
	{
        _ifs.getline(cdummy, 256);
		_ifs >> nu >> nv;
		_ifs.getline(cdummy, 256);
		for (j=0; j<nu*nv; ++j)
		{
            _ifs >> pt[0] >> pt[1] >> pt[2];
			vec_pt.push_back(pt);
			_ifs.getline(cdummy, 256);
		}

		psf = new FergusonSurface(nu, nv, vec_pt);
		vec_pt.clear();
		surface_geometry_.push_back(psf);
	}

	_ifs.getline(cdummy, 256); // "2.- Mesh Generation" label
	_ifs >> ncs >> nsr; 
    _ifs.getline(cdummy, 256);
	_ifs.getline(cdummy, 256); // "Segments in curves" label
	// the curve segment number/index starts with 0 here, while in the fli file it's 1
	for (i=0; i<ncs; ++i)
	{
		_ifs.getline(cdummy, 256);
		curve_segment_.push_back(i);
	}

    _ifs.getline(cdummy, 256); // "Regions on surfaces" label
    // the surface region number/index starts with 0 here, while in the fli file it's 1
	for (i=0; i<nsr; ++i)
	{
		_ifs >> idummy >> surface_id >> idummy;
        _ifs.getline(cdummy, 256);
		_ifs >> n_curve;
		_ifs.getline(cdummy, 256);
        curve_loop = new int [n_curve];
		for (j=0; j<n_curve; ++j)
			_ifs >> curve_loop[j];
		_ifs.getline(cdummy, 256);

		psr = new SurfaceRegion(surface_id, n_curve, curve_loop);
		surface_region_.push_back(psr);
	}


	// initialize the geometric entities once forever

    // first the curves
    for (i=0; i<ncv; ++i)
		curve_geometry_.at(i)->initialize(); // type conversion is necessary
	// second the surfaces
    for (i=0; i<nsf; ++i)
		surface_geometry_.at(i)->initialize(); // type conversion is necessary

    return true;
}