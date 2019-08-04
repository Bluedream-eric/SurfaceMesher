#ifndef GEOMETRYDEFINITION_H
#define GEOMETRYDEFINITION_H

//== INCLUDES =================================================================
#include <OpenMesh\Core\Geometry\MathDefs.hh>
#include <OpenMesh\Core\Geometry\VectorT.hh>
#include "MatrixT.h"
#include "BBox.h"

#include <assert.h>
#include <vector>
#include <fstream>
#include <string>

using namespace std;
using namespace OpenMesh;


//== CURVE CLASS DEFINITION ===================================================

/** \class Curve  GeometryDefinition.h

    Base class for different kinds of curves. It mainly contains some
	common interfaces for curves.
*/

// mutually exclusive, should not be put in class Curve
enum CurveType  
{
	CurveLine = 0  ,
    CurveFerguson  ,
	CurveNurbs     ,
	CurveTrimmed   
};

class Curve
{
public:
	
	//== constructors =========================================================

	Curve() : curve_type_(CurveLine), parametric_interval_(Vec2d(0.0, 1.0))    {}  

	Curve(CurveType _curve_type, Vec2d _parametric_interval = Vec2d(0.0, 1.0)) 
		: curve_type_(_curve_type), parametric_interval_(_parametric_interval) {}

	CurveType get_curve_type() { return curve_type_; }

	void set_curve_type(CurveType _curve_type) { curve_type_ = _curve_type; }

	Vec2d get_parametric_interval() { return parametric_interval_; }

	void set_parametric_interval(Vec2d _parametric_interval)
	{ parametric_interval_ = _parametric_interval; }

	//== abstract/common interfaces ===========================================

	virtual void initialize() = 0;//                                                {}

	virtual void evaluate(double _t, Vec3d& _c) const = 0;//                       {} // pure virtual function?

	virtual void calculate_first_derivative(double _t, Vec3d& _ct) const     {}

	// needed when the desired element sizes depend on the radii of curvature
	virtual void calculate_second_derivative(double _t, Vec3d& _ctt) const   {}

	virtual void calculate_projection(Vec3d _, double& _t) const             {}

	virtual void calculate_arc_length(double _t1, double _t2, double& _length) const  {}

private:
    CurveType curve_type_;
	Vec2d parametric_interval_;
};

/** \class LineCurve  GeometryDefinition.h

    Class for Line Curve. 
*/

class LineCurve : public Curve
{
public:

    //== constructors & basic functions =======================================

	LineCurve()  // pt_[0](Vec3d(0.0, 0.0, 0.0)), pt_[1](Vec3d(0.0, 0.0, 0.0)), compile error
	{ 
		set_curve_type(CurveLine); 
		set_parametric_interval(Vec2d(0.0, 1.0)); 
	    pt_[0] = Vec3d(0.0, 0.0, 0.0);
		pt_[1] = Vec3d(0.0, 0.0, 0.0);
	}

	LineCurve(Vec3d _pt[2])
    { 
		set_curve_type(CurveLine); 
		set_parametric_interval(Vec2d(0.0, 1.0)); 
	    pt_[0] = _pt[0];
		pt_[1] = _pt[1];
	}

	LineCurve(Vec3d _pt0, Vec3d _pt1)
	{ 
		set_curve_type(CurveLine); 
		set_parametric_interval(Vec2d(0.0, 1.0)); 
	    pt_[0] = _pt0;
		pt_[1] = _pt1;
	}

	void  set_first_point(const Vec3d& _pt)  {pt_[0] = _pt; }
	void  set_second_point(const Vec3d& _pt) {pt_[1] = _pt; }
	Vec3d get_first_point() const  { return pt_[0]; }
	Vec3d get_second_point() const { return pt_[1]; }

	//== functions ============================================================

	bool is_degenerate() const
	{ 
		double length = (pt_[1] - pt_[0]).norm();

		return is_zero(length); // tolerance: 1.0e-5
	}

	void initialize()
	{
		return;
	}

	void evaluate(double _t, Vec3d& _c) const
	{
        assert(0.0<=_t && _t<=1.0 ); // may not be necessary
		_c = pt_[0] + _t*(pt_[1] - pt_[0]);
	}

	void calculate_first_derivative(double _t, Vec3d& _ct) const
	{ _ct = pt_[1] - pt_[0]; } // warning_clg: what if vertical

	void calculate_second_derivative(double _t, Vec3d& _ctt) const
	{ _ctt = Vec3d(0.0, 0.0, 0.0); }

	void calculate_projection(Vec3d _c, double& _t) const
	{
        // assumption: _c is on this line and the line is not degenerate

        Vec3d line = pt_[1] - pt_[0];
		Vec3d vector = _c - pt_[0];
		_t = vector.norm() / line.norm(); // _t should be in [0.0, 1.0]
	}

	void calculate_arc_length(double _t1, double _t2, double& _length) const
	{
		// assumption: parameters are valid

		assert(0.0<=_t1 && _t1<=1.0);
		assert(0.0<=_t2 && _t2<=1.0);

		Vec3d line = pt_[1] - pt_[0];
		_length = fabs(_t2 - _t1) * line.norm();
	}

private:
	Vec3d pt_[2];
};

/** \class FergusonCurve  GeometryDefinition.h
    
    Class for Ferguson Curve
*/

class FergusonCurve : public Curve
{
public:

	//== constructors =========================================================

	FergusonCurve() 
	{ 
		set_curve_type(CurveFerguson); 
	    set_parametric_interval(Vec2d(0.0, 1.0));	
		control_point_.clear();
	}

	FergusonCurve(std::vector<Vec3d> _control_point) 
	{ 
		set_curve_type(CurveFerguson); 
        control_point_ = _control_point; // ok?
		int num = (int)_control_point.size();
		set_parametric_interval(Vec2d(0.0, num-1));
	}

	void set_control_point(std::vector<Vec3d> _control_point)
	{ control_point_ = _control_point; }

	int get_control_point_number() const { return (int)control_point_.size(); }

	std::vector<Vec3d>& get_control_point() { return control_point_; }

	std::vector<Vec3d>& get_control_point_tangent_vector() { return control_point_tangent_vector_; }

	std::vector<double>& get_segment_length() { return segment_length_; }

	std::vector<Vec5d>& get_segment_coefficient() { return segment_coefficient_; }

	//== main external functions ==============================================

	/** \attention: the following functions depend on the parameterization, i.e., the 
	    parametric curve expression. In BLSURF, it's provided by external modules.

	    \remark: here we may adopt the normative parameterization, i.e., for each curve segment,  
	    the parameter domain is [0, 1] and the entire composite curve is [0, n], where n is 
	    the number of curve segments.
    */
	void initialize();
	
    void evaluate(double _t, Vec3d& _c) const;

	void calculate_first_derivative(double _t, Vec3d& _ct) const;

	void calculate_second_derivative(double _t, Vec3d& _ctt) const;

	void calculate_projection(Vec3d _c, double& _t) const;

	void calculate_arc_length(double _t1, double _t2, double& _length) const;

    //== other external functions =============================================

	enum EndCondition 
	{
		ZeroSecondDerivative = 0,  // free form end
		Parabolic,                 // parabolic form end
		SpecifiedTangent           // clipping form end
	};
	void calculate_control_point_tangent_vector(EndCondition _bc=ZeroSecondDerivative, 
		Vec3d _v0=Vec3d(0.0,0.0,0.0), Vec3d _vn=Vec3d(0.0,0.0,0.0));  // _t1, _t2 for specified tangent

	// gets the coefficients of the polynomial |r'|**2
	// maybe it should be private cause it only serves for curve length calculation
	void calculate_segment_coefficient(); 

	// calculate the length segment by segment
	void calculate_segment_length(double _tol=1.e-6);
	void calculate_partial_segment_length(const Vec5d& _coef, double _u0, double _u1, 
		double& _len, double _tol=1.e-6);

	// calculates the position _u1 of a point on a single ferguson curve segment 
	// such that the length of the cubic segment _u0,_u1 is _tar_len                 
    void calculate_parameter_given_arc_length(const Vec5d& coef, double _seg_len, double _u0, 
		double& _u1, double _tar_len, double _tol=1.0e-6);

	enum DerivativeNumber
	{
        Zero = 0,
		One   ,
		Two
	};
	void evaluate_ferguson_curve_segment(DerivativeNumber _dn, Vec3d _v0[2], Vec3d _v1[2], double _u, 
		Vec3d _res[3]) const;

private:
	std::vector<Vec3d> control_point_;  // ordered control points

	// the following are kind of intermediate results but important and convenient, 
	// which need to be computed. This class should be checked thoroughly

	// the tangent vectors are determineed by the type of end condition, and 
	// the entire ferguson curve is totally determined by the control points and their tangent vector
	std::vector<Vec3d> control_point_tangent_vector_;

	// the following are used for calculating the length of the curve, 
	// should they be put in SurfaceMesher class?

	std::vector<Vec5d> segment_coefficient_;
	// here we use accumulating length just like psue did 
	std::vector<double> segment_length_; 
};

/** \class NurbsCurve  GeometryDefinition.h

    Class for NURBS curve
*/

class NurbsCurve : public Curve
{
public:
	// todo_clg

	//== constructors =========================================================

	NurbsCurve() 
	{ set_curve_type(CurveNurbs); set_parametric_interval(Vec2d(0.0, 1.0)); }

private:
	// todo_clg
};

/** \class TrimmedCurve  GeometryDefinition.h

    Class for trimmed curve
*/

class TrimmedCurve : public Curve
{
	// todo_clg
};


//== SURFACE CLASS DEFINITION =================================================

/** \class Surface  GeometryDefinition.h
    
	Base class for different surface types. It mainly contains some
	common interfaces for surfaces.
*/

enum SurfaceType  // mutually exclusive, be aware of the name domain.
{
	SurfaceLinearTriangle = 0  ,  // quadrilateral or other polygon
	SurfaceQuadraticTriangle   ,  // quadratic quadrilateral
	SurfaceFerguson            ,  
    SurfaceNurbs               ,
	SurfaceTrimmed 
};

class Surface 
{
public:

	//== constructors =========================================================

	Surface() : surface_type_(SurfaceLinearTriangle) {}  // default is plane surface

	Surface(SurfaceType _surface_type) : surface_type_(_surface_type) {}

	void set_surface_type(SurfaceType _surface_type) { surface_type_ = _surface_type; }
	SurfaceType get_surface_type() { return surface_type_; }

    //== abstract interfaces ==================================================

	// \remark: most of the following functions depend on the particular parameterization.
	// evaluate the influence of this surface parameterization.

	virtual void initialize()                                                {}

	virtual void evaluate(Vec2d _uv, Vec3d& _s) const                        {}

	// equivalent to get the coefficients of the first fundamental form, E, F, G.
	virtual void calculate_first_derivative(Vec2d _uv, Vec3d& _su, Vec3d& _sv) const  {}

	// needed when the desired element sizes depend on the radii of curvature
	// equivalent to get the coefficients of the second fundamental form, L, M, N.
	virtual void calculate_second_derivative(const Vec2d& _uv, Vec3d& _suu, 
		Vec3d& _suv, Vec3d& _svv) const                                      {}  

	// calculate the first and second derivatives in one function for efficiency consideration
	virtual void calculate_first_and_second_derivative(Vec2d _uv, Vec3d& _su, 
		Vec3d& _sv, Vec3d& _suu, Vec3d& _suv, Vec3d& _svv) const             {}

	// needed for curvature element-size control, either isotropic or anisotropic.
	// calculated through the first and second derivatives.
    virtual void calculate_pricipal_curvature_and_direction(Vec2d _uv, double& _kmin, 
		double& _kmax, Vec3d& _dmin, Vec3d& _dmax) const                     {}

	// at least _s is close to the surface enough, can be solved by the 
	// inverse analytic function or an iterative procedure.
	virtual void calculate_projection(Vec3d _s, Vec2d& _uv) const            {}

private:
    SurfaceType surface_type_;

	// shouldn't put here, it's definition of support surface, not for mesh generation
	// std::vector<Curve*> trimming_curve_; 
};

class LinearTriangleSurface : public Surface
{
public:

	//== constructors & basic functions =======================================

	LinearTriangleSurface() 
	{
		set_surface_type(SurfaceLinearTriangle);
		for (int i = 0; i < 3; i++)
			pt_[i] = Vec3d(0.0, 0.0, 0.0);

		coef_matrix_.reset();
		metric_tensor_.reset();
	}

	LinearTriangleSurface(Vec3d _pt0, Vec3d _pt1, Vec3d _pt2)
	{
		set_surface_type(SurfaceLinearTriangle);
		pt_[0] = _pt0; pt_[1] = _pt1; pt_[2] = _pt2;

		calc_coef_matrix();

		Vec3d ru = pt_[1] - pt_[0], rv = pt_[2] - pt_[0];
        metric_tensor_ = Mat2by2d(ru | ru, ru | rv, ru | rv, rv | rv);
	}

	LinearTriangleSurface(Vec3d _pt[3])
	{ LinearTriangleSurface(_pt[0], _pt[1], _pt[2]); }

	~LinearTriangleSurface()
	{
	}


	void set_point(const Vec3d& _pt0, const Vec3d& _pt1, const Vec3d& _pt2)
	{
		pt_[0] = _pt0; pt_[1] = _pt1; pt_[2] = _pt2;
		calc_coef_matrix();
	}

	void set_point(Vec3d _pt[3])
	{
		for (int i=0; i<3; ++i)
			pt_[i] = _pt[i];

		calc_coef_matrix();

		Vec3d ru = pt_[1] - pt_[0], rv = pt_[2] - pt_[0];
        metric_tensor_ = Mat2by2d(ru | ru, ru | rv, ru | rv, rv | rv);
	}

	Vec3d* get_point()
	{ return pt_; }

	void get_point(Vec3d& _pt0, Vec3d& _pt1, Vec3d& _pt2) const
	{ 
		_pt0 = pt_[0]; _pt1 = pt_[1]; _pt2 = pt_[2];
	}
	
	void get_point(Vec3d _pt[3]) const
	{
        for (int i=0; i<3; ++i)
			_pt[i] = pt_[i];
	}

	/**
	const Mat2by2d& get_metric_tensor() const // compile error, redefinition
	{
		return metric_tensor_;
	}
	*/

	Mat2by2d get_metric_tensor() const
	{
		return metric_tensor_;
	}

	//== functions ============================================================

	/** \remark: the particular parameterization (parametric domain) we adopt here
	    is the triangle with three nodes (0.0, 0.0), (1.0, 0.0) and (0.0, 1.0), 
		the linear mapping from parametric triangle to physical triangle can be 
		easily determined.
		we can also use a different parameterization with equilateral parametric 
		triangle with nodes (0.0, 0.0), (1.0, 0.0) and (1/2, sqrt(3)/2).
	*/ 

    void evaluate(const Vec2d& _uv, Vec3d& _s) const
	{
		Vec3d a = pt_[1] - pt_[0], b = pt_[2] - pt_[0], c = pt_[0];
		_s = _uv[0]*a + _uv[1]*b +c;
	}

	void calculate_first_derivative(Vec2d _uv, Vec3d& _su, Vec3d& _sv) const // shouldn't be Vec2d&
	{
        Vec3d a = pt_[1] - pt_[0], b = pt_[2] - pt_[0];
        _su = a; _sv = b;
	}

	void calculate_second_derivative(const Vec2d& _uv, Vec3d& _suu, Vec3d& _suv, 
		    Vec3d& _svv) const
	{
		_suu = _suv = _svv = Vec3d(0.0, 0.0, 0.0);
	}

	void calculate_first_and_second_derivative(Vec2d _uv, Vec3d& _su, 
		    Vec3d& _sv, Vec3d& _suu, Vec3d& _suv, Vec3d& _svv) const             
	{
        Vec3d a = pt_[1] - pt_[0], b = pt_[2] - pt_[0];
        _su = a; _sv = b;
		_suu = _suv = _svv = Vec3d(0.0, 0.0, 0.0);
	}

	void calculate_pricipal_curvature_and_direction(const Vec2d& _uv, double& _kmin, 
		    double& _kmax, Vec3d& _dmin, Vec3d& _dmax) const
	{
        // actually not needed. _dmin and _dmax all the same.
		_kmin = _kmax = 0.0;
	}

	void calculate_projection(const Vec3d& _s, Vec2d& _uv)
	{
		Mat3by1d s(_s[0], _s[1], _s[2]), uv1; // conversion operator from vector to matrix
		uv1 = coef_matrix_ * s;
		_uv[0] = uv1(0, 0);
		_uv[1] = uv1(1, 0);
	}

private:
	void calc_coef_matrix()
	{
        double mat[3][3] = {{(pt_[1]-pt_[0])[0], (pt_[2]-pt_[0])[0], pt_[0][0]}, 
		    {(pt_[1]-pt_[0])[1], (pt_[2]-pt_[0])[1], pt_[0][1]}, 
		    {(pt_[1]-pt_[0])[2], (pt_[2]-pt_[0])[2], pt_[0][2]}};
		for (int i=0; i<3; ++i)
			for (int j=0; j<3; ++j)
				coef_matrix_(i,j) = mat[i][j];
		coef_matrix_.inversion();
	}

private:
	Vec3d pt_[3];

	Mat3by3d coef_matrix_;  // auxiliary inverse matrix of (pt_[1]-pt_[0],pt_[2]-pt_[0],pt_[0]) 

	Mat2by2d metric_tensor_;
};

class FergusonSurface : public Surface
{
public:

	/**
	typedef struct ControlPoint
	{
        Vec3d r, ru, rv, ruv;
	}ControlPoint;
	*/

	//== constructors & basic functions =======================================

	FergusonSurface()
	{
		set_surface_type(SurfaceFerguson);
		nu_ = nv_ = 0;
		control_point_.clear();
	}

	FergusonSurface(int _nu, int _nv, std::vector<Vec3d> _control_point)
		: nu_(_nu), nv_(_nv)
	{
		assert(_control_point.size() == _nu*_nv);
		set_surface_type(SurfaceFerguson);
        control_point_ = _control_point; // warning_clg: could be probably wrong!
	}

	~FergusonSurface()
	{
	}

	void set(int _nu, int _nv, std::vector<Vec3d> _control_point)
	{
		assert(_control_point.size()==_nu*_nv);

		nu_ = _nu;
		nv_ = _nv;
		control_point_ = _control_point;
	}

	int get_nu() const { return nu_; }

	int get_nv() const { return nv_; }

	std::vector<Vec3d>& get_control_point()  { return control_point_; }

	//== functions ============================================================

	/** \remark: the particular parameterization (parametric domain) we adopt here
	    is the quadrilateral with four nodes (0.0, 0.0), (nu_-1, 0.0), (nu_-1, nv_-1) 
		and (0.0, nv_-1), the function value, the first and the second derivatives and 
		so on can be easily determined by the following functions.
	*/ 

	// calculate necessary information
	void initialize();

	void evaluate(Vec2d _uv, Vec3d& _s) const;

	void calculate_first_derivative(Vec2d _uv, Vec3d& _su, Vec3d& _sv) const;

	void calculate_second_derivative(const Vec2d& _uv, Vec3d& _suu, Vec3d& _suv, 
		Vec3d& _svv) const;  
	
	void calculate_first_and_second_derivative(Vec2d _uv, Vec3d& _su, 
		Vec3d& _sv, Vec3d& _suu, Vec3d& _suv, Vec3d& _svv) const;

    void calculate_pricipal_curvature_and_direction(Vec2d _uv, double& _kmin, 
		double& _kmax, Vec3d& _dmin, Vec3d& _dmax) const;

	void calculate_projection(Vec3d _s, Vec2d& _uv) const;

private:

	//== inner called functions for initialize() ==============================

	// according to control points r, calculate ru, rv, ruv, to get the ferguson surface
	void calculate_control_point_derivative();

	enum EndCondition 
	{
		ZeroSecondDerivative = 0,  // free form end
		Parabolic,                 // parabolic form end
		SpecifiedTangent           // clipping form end
	};
	void calculate_ferguson_curve_control_point_tangent_vector(std::vector<Vec3d>& _control_point,
		std::vector<Vec3d>& _tangent_vector, EndCondition _bc=ZeroSecondDerivative, 
		Vec3d _v0=Vec3d(0.0,0.0,0.0), Vec3d _vn=Vec3d(0.0,0.0,0.0));  // _t1, _t2 for specified tangent

private:
	std::vector<Vec3d> control_point_;
	int nu_, nv_; // parametric domain: [0,nu_-1]*[0,nv_-1]

	// necessary information for the surface, need to be computed
	// ru, rv, ruv for control points
	std::vector<Vec3d> ru_;
	std::vector<Vec3d> rv_;
	std::vector<Vec3d> ruv_;
};


//== GEOMETRY CLASS DEFINITION =================================================

/** \class Geometry  GeometryDefinition.h
    
	class for a geometry model for surface mesh generation.
*/

typedef struct SurfaceRegion
{
    int surface_id;   // support surface id
    int curve_number; // number of curve segments in the loop

	// curve segment ids of the loop, 
	// may not use std::vector for it's used in std::vector<SurfaceRegion>
	int* curve_loop;  

	SurfaceRegion()
	{
		surface_id = 0;
		curve_number = 0;
		curve_loop = NULL;
	}

	SurfaceRegion(int _surface_id, int _curve_number, int* _curve_loop)
	{
		surface_id = _surface_id;
		curve_number = _curve_number;
		curve_loop = _curve_loop;
	}

	~SurfaceRegion()  // free memory/resource
	{
		if (curve_loop!=NULL)
			delete curve_loop;
	}
} SurfaceRegion;

class Geometry
{
public:

	//== constructors & destructors & basic functions =========================

	Geometry()                                
	{
		curve_geometry_.clear();
		surface_geometry_.clear();
		curve_segment_.clear();
		surface_region_.clear(); 

		bbox_.reset();
		model_type_.clear();
	}

	~Geometry()  // free memory/resource
	{
		int i; 

		for (i=0; i<(int)curve_geometry_.size(); ++i)
			delete curve_geometry_[i];
		curve_geometry_.clear();

		for (i=0; i<(int)surface_geometry_.size(); ++i)
			delete surface_geometry_[i];
		surface_geometry_.clear();

		for (i=0; i<(int)surface_region_.size(); ++i)
			delete surface_region_[i];
		surface_region_.clear();
	}

    int get_curve_geometry_number() const
	{ return (int)curve_geometry_.size(); }

	int get_surface_geometry_number() const
	{ return (int)surface_geometry_.size(); }

	int get_curve_segment_number() const
	{ return (int)curve_segment_.size(); }

	int get_surface_region_number() const
	{ return (int)surface_region_.size(); }

	/**
	int get_trimming_curve_number(int _surface_region_id)
	{
		assert(_surface_region_id>=0 && _surface_region_id<get_surface_region_number());
		return surface_region_[_surface_region_id].curve_number;
	}
	*/

	std::vector<Curve*>& get_curve_geometry()  // not const 
	{ return curve_geometry_; }

	std::vector<Surface*>& get_surface_geometry()
	{ return surface_geometry_; }

	std::vector<int>& get_curve_segment()
	{ return curve_segment_; }

	std::vector<SurfaceRegion*>& get_surface_region()
	{ return surface_region_; }

	void append_curve_geometry(Curve* _pcg)  // pass value
	{
		assert(_pcg!=NULL);
		curve_geometry_.push_back(_pcg); 
	}

	void append_surface_geometry(Surface* _psg)
	{
		assert(_psg!=NULL);
		surface_geometry_.push_back(_psg); 
	}

	void append_curve_segment(int _cs)
	{
		assert(_cs>=0);
		curve_segment_.push_back(_cs); 
	}

	void append_surface_region(SurfaceRegion* _psr)
	{
		assert(_psr!=NULL);
		surface_region_.push_back(_psr); 
	}

	const std::string& get_model_type() const
	{ return model_type_; }

	//== I/O with given file formats ==========================================

	bool read(const std::string& _filename);

private:
	std::string get_file_extension(const std::string& _filename) const;
	bool read_from_stl(std::ifstream& _ifs);
	bool read_from_fli(std::ifstream& _ifs);

	BBox bbox_; // haven't been used
	std::string model_type_;

private:
	// geometry definition
	std::vector<Curve*>    curve_geometry_;    // curve, id <=> array index
	std::vector<Surface*>  surface_geometry_;  // support surface

	// mesh generation information
	std::vector<int>       curve_segment_;     // reference to curve
	std::vector<SurfaceRegion*> surface_region_;
};

#endif  // GEOMETRYDEFINITION_H defined