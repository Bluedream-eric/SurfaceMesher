#include "MeshSizeSpecification.h"
#include <fstream>

using namespace std;

bool MeshSizeSpecification::read_from_ba3(std::string& _filename)
{
	std::ifstream ba3_file(_filename.c_str());
	assert(ba3_file);

	int np, ne, nps, nls, nts;
	int i, j;
	double ddummy;
	char cdummy[1024];
	// Node nd;
	// Elem el;
	PointSource    pt_src;
	LineSource     ln_src;
	TriangleSource tri_src;

	ba3_file.getline(cdummy, 1024);
	ba3_file >> np >> ne >> nps >> nls >> nts;
	ba3_file.getline(cdummy, 1024); // never forget!
	cout<<"    #PtBkg #ElBkg #PtSr #LnSr #TrSr: "<<np<<" "<<ne<<" "<<nps<<" "<<nls<<" "<<nts<<endl; 
	
	// filter with bkg
	for (i=0; i<np; i++)
	{
		ba3_file.getline(cdummy, 1024);
		ba3_file >> ddummy >> ddummy >> ddummy >> constant_spacing_;
		ba3_file.getline(cdummy, 1024);
		ba3_file.getline(cdummy, 1024);
		ba3_file.getline(cdummy, 1024);
	}
	for (i=0; i<ne; i++)
		ba3_file.getline(cdummy, 1024);

	ba3_file.getline(cdummy, 1024);  // *points
	for (i=0; i<nps; i++)
	{
		ba3_file.getline(cdummy, 1024);  // index, sometimes it has other words
		ba3_file >> pt_src.center[0] >> pt_src.center[1] >> pt_src.center[2]
		         >> pt_src.spacing >> pt_src.inner_radius >> pt_src.outer_radius;
		ba3_file.getline(cdummy, 1024);
		point_source_.push_back(pt_src);
	}

	ba3_file.getline(cdummy, 1024);  // *lines
	for (i=0; i<nls; i++)
	{
		ba3_file.getline(cdummy, 1024);  // index
		for (j=0; j<2; j++)
		{
			ba3_file >> ln_src.point_source[j].center[0] >> ln_src.point_source[j].center[1] 
			    >> ln_src.point_source[j].center[2] >> ln_src.point_source[j].spacing
				>> ln_src.point_source[j].inner_radius >> ln_src.point_source[j].outer_radius;
			ba3_file.getline(cdummy, 1024); // never forget!
		}
		line_source_.push_back(ln_src);
	}

	ba3_file.getline(cdummy, 1024);  // * triangles
	for (i=0; i<nts; i++)
	{
		ba3_file.getline(cdummy, 1024); // index
		for (j=0; j<3; j++)
		{
			ba3_file >> tri_src.point_source[j].center[0] >> tri_src.point_source[j].center[1]
			    >> tri_src.point_source[j].center[2] >> tri_src.point_source[j].spacing
				>>tri_src.point_source[j].inner_radius >> tri_src.point_source[j].outer_radius;
			ba3_file.getline(cdummy, 1024);
		}
		triangle_source_.push_back(tri_src);
	}

	ba3_file.close();

	smooth_mesh_gradation();
	calculate_min_spacing();

	return true;
}

void MeshSizeSpecification::smooth_mesh_gradation()
{
	int i, j;
	const double factor = 0.35, tolg = 1.e-4;
	double fac, f;

    // ignore to modify the background mesh temporally

    // Modifies the source distribution
	cout << "SurfaceMesher > Smooth mesh gradation, ignore to modify the background mesh temporally" << endl;
	cout << "SurfaceMesher > Modifying source distribution" << endl;
    
	fac = factor / log(2.0);

	for (i=0; i<(int)point_source_.size(); ++i)
	{
		point_source_[i].outer_radius = point_source_[i].outer_radius > (point_source_[i].inner_radius + tolg) ?
			point_source_[i].outer_radius : (point_source_[i].inner_radius + tolg);
		f = point_source_[i].spacing / (point_source_[i].outer_radius - point_source_[i].inner_radius);
		if (f>fac)
		{
			point_source_[i].outer_radius = point_source_[i].inner_radius + point_source_[i].spacing / fac;
			cout << "Point source " << i+1 << " has been modified." << endl;
		}
	}

    for (i=0; i<(int)line_source_.size(); ++i)
	{
		for (j=0; j<2; ++j)
		{
			line_source_[i].point_source[j].outer_radius = line_source_[i].point_source[j].outer_radius >
				(line_source_[i].point_source[j].inner_radius + tolg) ? line_source_[i].point_source[j].outer_radius :
                (line_source_[i].point_source[j].inner_radius + tolg);
			f = line_source_[i].point_source[j].spacing / (line_source_[i].point_source[j].outer_radius - line_source_[i].point_source[j].inner_radius);
			if (f>fac)
			{
                line_source_[i].point_source[j].outer_radius = line_source_[i].point_source[j].inner_radius + line_source_[i].point_source[j].spacing / fac;
				cout << "Line source " << i+1 << " has been modified." << endl;
			}
		}
	}

	for (i=0; i<(int)triangle_source_.size(); ++i)
	{
		for (j=0; j<3; ++j)
		{
            triangle_source_[i].point_source[j].outer_radius = triangle_source_[i].point_source[j].outer_radius >
				(triangle_source_[i].point_source[j].inner_radius + tolg) ? triangle_source_[i].point_source[j].outer_radius :
                (triangle_source_[i].point_source[j].inner_radius + tolg);
			f = triangle_source_[i].point_source[j].spacing / (triangle_source_[i].point_source[j].outer_radius - triangle_source_[i].point_source[j].inner_radius);
			if (f>fac)
			{
                triangle_source_[i].point_source[j].outer_radius = triangle_source_[i].point_source[j].inner_radius + triangle_source_[i].point_source[j].spacing / fac;
				cout << "Triangle source " << i+1 << " has been modified." << endl;
			}
		}
	}

	return;
}

void MeshSizeSpecification::calculate_min_spacing()
{
    int i, j;
	double tmp_spa;
	min_spacing_ = constant_spacing_;

    for (i=0; i<(int)point_source_.size(); ++i)
	{
		tmp_spa = point_source_[i].spacing;
		if (tmp_spa<min_spacing_)
			min_spacing_ = tmp_spa;
	}

	for (i=0; i<(int)line_source_.size(); ++i)
	{
		for (j=0; j<2; ++j)
		{
			tmp_spa = line_source_[i].point_source[j].spacing;
		    if (tmp_spa<min_spacing_)
			    min_spacing_ = tmp_spa;
		}
	}

	for (i=0; i<(int)triangle_source_.size(); ++i)
	{
		for (j=0; j<3; ++j)
		{
			tmp_spa = triangle_source_[i].point_source[j].spacing;
		    if (tmp_spa<min_spacing_)
			    min_spacing_ = tmp_spa;
		}
	}

	assert(min_spacing_>1.0e-6);

	return;
}

double MeshSizeSpecification::get_final_spacing(Vec3d &_pt)
{
    double tmp_spa, res_spa = 1.e30;
    tmp_spa = get_spacing_from_mesh_source(_pt);
	if (tmp_spa < res_spa)
		res_spa = tmp_spa;
	if (constant_spacing_ < res_spa)
		res_spa = constant_spacing_;

	return res_spa;
}

double MeshSizeSpecification::get_spacing_from_mesh_source(Vec3d& _pt)
{
    int i;
	double tmp_spa, res_spa = 1.e30;
	
	for (i=0; i<(int)point_source_.size(); ++i)
	{
		tmp_spa = get_spacing_from_point_source(point_source_[i], _pt);
		if (tmp_spa<res_spa)
			res_spa = tmp_spa;
	}

	for (i=0; i<(int)line_source_.size(); ++i)
	{
		tmp_spa = get_spacing_from_line_source(line_source_[i], _pt);
		if (tmp_spa<res_spa)
			res_spa = tmp_spa;
	}

	for (i=0; i<(int)triangle_source_.size(); ++i)
	{
		tmp_spa = get_spacing_from_triangle_source(triangle_source_[i], _pt);
		if (tmp_spa<res_spa)
			res_spa = tmp_spa;
	}

	return res_spa;
}

double MeshSizeSpecification::get_spacing_from_point_source(PointSource& _pnt_src, Vec3d& _pt) const
{
    double diff, ae;
	const double BIG = 50.0;
	const double CLG2 = log(2.0);
	Vec3d dist_vec = _pnt_src.center - _pt;
	double dist = dist_vec.norm();
	if (dist<=_pnt_src.inner_radius)
		return _pnt_src.spacing;
	else
	{
		dist = dist - _pnt_src.inner_radius;
		diff = _pnt_src.outer_radius - _pnt_src.inner_radius;
		assert(diff>0.0);
		diff = CLG2 / diff;
		ae = dist * diff < BIG ? dist * diff : BIG;

		return _pnt_src.spacing * exp(ae);
	}
}
	
double MeshSizeSpecification::get_spacing_from_line_source(LineSource& _ln_src, Vec3d& _pt) const
{
    double proj_len = 0.0;  // signed project length
	double w1, w2;
	// REAL tolg = 1e-5;
	Vec3d v01, line = _ln_src.point_source[1].center - _ln_src.point_source[0].center;
	PointSource pnt_src2;

	double dist = line.norm(); 
	assert(dist>1.0e-6);  // avoid coincident points and degenerate line 
	v01 = line / dist;

	proj_len = (_pt - _ln_src.point_source[0].center) | v01;
	if (proj_len<=0.0)
		return get_spacing_from_point_source(_ln_src.point_source[0], _pt);
	else if (proj_len>=dist) 
		return get_spacing_from_point_source(_ln_src.point_source[1], _pt);
	else
	{
		w2 = proj_len / dist;
		w1 = 1.0 - w2;

		pnt_src2.center = w1 * _ln_src.point_source[0].center + w2 * _ln_src.point_source[1].center;
		pnt_src2.inner_radius = w1*_ln_src.point_source[0].inner_radius + w2*_ln_src.point_source[1].inner_radius;
		pnt_src2.outer_radius = w1*_ln_src.point_source[0].outer_radius + w2*_ln_src.point_source[1].outer_radius;
		pnt_src2.spacing = w1*_ln_src.point_source[0].spacing + w2*_ln_src.point_source[1].spacing;

		return get_spacing_from_point_source(pnt_src2, _pt);
	}
}

double MeshSizeSpecification::get_spacing_from_triangle_source(TriangleSource& _tri_src, Vec3d& _pt) const
{
	double w0, w1, w2, dist;
	// REAL tolg = 1e-5;
	const double EPS_ZERO_SQ = 1.e-12;
	PointSource pnt_src3;
	LineSource  ln_src3;
	Vec3d v01, v12, v20;
	Vec3d ntri, n01, n12; 
	Vec3d diff_vec;  // _pt - _tri_src.point_source[1]

	v01 = _tri_src.point_source[1].center - _tri_src.point_source[0].center;
	dist = v01.sqrnorm(); 
	assert(dist>EPS_ZERO_SQ);

	v12 = _tri_src.point_source[2].center - _tri_src.point_source[1].center;
	dist = v12.sqrnorm(); 
	assert(dist>EPS_ZERO_SQ);

	v20 = _tri_src.point_source[0].center - _tri_src.point_source[2].center;
	dist = v20.sqrnorm(); 
	assert(dist>EPS_ZERO_SQ);

	ntri = v01 % v12;
	n01 = ntri % v01;
	n12 = ntri % v12;
	diff_vec = _pt - _tri_src.point_source[1].center;
	
	// compute area weights, using the point projected to the triangle 
    w0 = (diff_vec | n12) / (-v01 | n12); // the ratio of two heights of two triangles
	w2 = (diff_vec | n01) / ( v12 | n01); // detailed in the notebook
	w1 = 1.0 - w0 - w2;	
	
	if (w0<=0.0)
	{
		ln_src3.point_source[0] = _tri_src.point_source[1];
		ln_src3.point_source[1] = _tri_src.point_source[2];

		return get_spacing_from_line_source(ln_src3, _pt);
	}
	else if (w1 <= 0.0)
	{
		ln_src3.point_source[0] = _tri_src.point_source[2];
		ln_src3.point_source[1] = _tri_src.point_source[0];

		return get_spacing_from_line_source(ln_src3, _pt);
	}
	else if (w2<=0.0)
	{
        ln_src3.point_source[0] = _tri_src.point_source[0];
		ln_src3.point_source[1] = _tri_src.point_source[1];

		return get_spacing_from_line_source(ln_src3, _pt);
	}
	else
	{
		pnt_src3.center = w0*_tri_src.point_source[0].center + w1*_tri_src.point_source[1].center
			+ w2*_tri_src.point_source[2].center;
		pnt_src3.inner_radius = w0*_tri_src.point_source[0].inner_radius 
			+ w1*_tri_src.point_source[1].inner_radius + w2*_tri_src.point_source[2].inner_radius;
		pnt_src3.outer_radius = w0*_tri_src.point_source[0].outer_radius
			+ w1*_tri_src.point_source[1].outer_radius + w2*_tri_src.point_source[2].outer_radius;
		pnt_src3.spacing = w0*_tri_src.point_source[0].spacing
			+ w1*_tri_src.point_source[1].spacing + w2*_tri_src.point_source[2].spacing;

		return get_spacing_from_point_source(pnt_src3, _pt);
	}
}