#include "DTIso2D.h"  // for DTIso2D
#include "SurfaceMesher.h"
#include <list>
#include <deque>

using namespace std;
// using namespace OpenMesh;

OpenMesh::EPropHandleT<double> riemannian_length;
OpenMesh::EPropHandleT<std::list<MeshPoint3D> *> discretized_point;   
OpenMesh::VPropHandleT<double> iso3d_size;
OpenMesh::VPropHandleT<int> global_id;
// OpenMesh::EPropHandleT<double> quality_ratio;

void SurfaceMesher::riemannian_length_3d_line(Vec3d& _pt1, Vec3d& _pt2, double& _rl) const
{
	_rl = (_pt1-_pt2).norm();

	double h1, h2, h3, h;
	h1 = mesh_size_spec_->get_final_spacing(_pt1);
	h2 = mesh_size_spec_->get_final_spacing((_pt1+_pt2)/2.0);
	h3 = mesh_size_spec_->get_final_spacing(_pt2);
    h = (1.0/h1 + 4.0/h2 + 1.0/h3) / 6;

	_rl *= h;

	return;
}

// _mse has been allocated memory outside this function
void SurfaceMesher::create_binary_tree(MidsplitedEdge *_mse) const
{
	const double threshold_length = 0.1; // _tol!
	if (_mse->riemannian_length <= threshold_length)  // recursive termination condition
		return;

	_mse->left = new MidsplitedEdge;
	_mse->right = new MidsplitedEdge;

	_mse->left->point[0] = _mse->point[0];
	_mse->left->point[1] = (_mse->point[0]+_mse->point[1])/2.0;
	_mse->left->left = _mse->left->right = NULL;
	_mse->left->parent = _mse;
	riemannian_length_3d_line(_mse->left->point[0], _mse->left->point[1], _mse->left->riemannian_length);

    _mse->right->point[0] = (_mse->point[0]+_mse->point[1])/2.0;
	_mse->right->point[1] = _mse->point[1];
	_mse->right->left = _mse->right->right = NULL;
	_mse->right->parent = _mse;
	riemannian_length_3d_line(_mse->right->point[0], _mse->right->point[1], _mse->right->riemannian_length);

	create_binary_tree(_mse->left);
	create_binary_tree(_mse->right);
}

void SurfaceMesher::traverse_to_get_midsplited_point(MidsplitedEdge *_mse, 
	    std::vector<MidsplitedPoint> &_midsplited_point)
{
    assert(NULL!=_mse);

	// leaf node, recursive termination condition
	if (NULL==_mse->left && NULL==_mse->right) 
	{
        MidsplitedPoint msp;
	    msp.point = _mse->point[1];
		msp.riemannian_length = _mse->riemannian_length;
		_midsplited_point.push_back(msp);

		return;
	}
	else
	{
		traverse_to_get_midsplited_point(_mse->left, _midsplited_point);
		traverse_to_get_midsplited_point(_mse->right, _midsplited_point);
	}

	// free memory!!! do it later!

	return;
}

void SurfaceMesher::discretize_line_curve(LineCurve *_lc, double _tol)
{
	// create binary tree to generate intermediate auxiliary points
    MidsplitedEdge *root = new MidsplitedEdge;
	root->point[0] = _lc->get_first_point();
	root->point[1] = _lc->get_second_point();
	root->left = root->right = NULL;
	root->parent = NULL;
	riemannian_length_3d_line(root->point[0], root->point[1], root->riemannian_length);
    create_binary_tree(root);

	// store the intermediate auxiliary points as a vector
	std::vector<MidsplitedPoint> midsplited_point;
	MidsplitedPoint msp;
	msp.point = root->point[0];
	msp.riemannian_length = 0.0;
	midsplited_point.push_back(msp);
	// may need to output the results to test the correctness
	traverse_to_get_midsplited_point(root, midsplited_point); 

	cout << "    Number of midsplited points: " << midsplited_point.size() << endl;

    // determine the total riemannian length, ideal edge number, ideal riemannian length
	double ideal_riemannian_length = 0.0, total_riemannian_length = 0.0;
	int ideal_edge_number = 1;
	int i = 0;
	for (i=0; i<(int)midsplited_point.size(); ++i)
	{
		total_riemannian_length += midsplited_point[i].riemannian_length;
		midsplited_point[i].isotropic_size = mesh_size_spec_->get_final_spacing(midsplited_point[i].point);
		if (i > 0)
			midsplited_point[i].riemannian_length += midsplited_point[i-1].riemannian_length;
	}
    ideal_edge_number = (int)total_riemannian_length;
	if (total_riemannian_length - ideal_edge_number >= 0.5)
		ideal_edge_number++;
	if (0==ideal_edge_number)
		ideal_edge_number = 1;
	ideal_riemannian_length = total_riemannian_length / (1.0 * ideal_edge_number);

	cout << "    Number of final discretized points: " << ideal_edge_number+1 << endl;
    
	// generate final discretized points
	std::vector<MeshPoint3D> *vec_mp3d = new std::vector<MeshPoint3D>(ideal_edge_number+1);
	vec_mp3d->at(0).point = midsplited_point[0].point;
	vec_mp3d->at(0).isotropic_size = midsplited_point[0].isotropic_size;
	vec_mp3d->at(ideal_edge_number).point = midsplited_point[midsplited_point.size()-1].point;
	vec_mp3d->at(ideal_edge_number).isotropic_size = midsplited_point[midsplited_point.size()-1].isotropic_size;

	int j, k = 0, m = 0;
	double interval, c, f0, fc, temp_ideal_riemannian_length = 0.0, temp_mid_riemannian_length = 0.0;
	const double threshold_length = 1.e-6; // _tol!
    for (i=1; i<=ideal_edge_number-1; ++i)
	{
		for (j=k; j<(int)midsplited_point.size()-1; ++j)
		{
			if (midsplited_point[j].riemannian_length < i * ideal_riemannian_length 
				&& i * ideal_riemannian_length <= midsplited_point[j+1].riemannian_length)
			{
				k = j;
				break;
			}
		}
		assert(midsplited_point.size()-1 != j); // find the index

		// get the final discretized point by recursion
		// attension and todo_clg: the efficiency should be quite low because of the 
		// iteration, it need not to be such precise, it should be treated more roughly
		if (fabs(i*ideal_riemannian_length-midsplited_point[j].riemannian_length) < threshold_length)
		{
			vec_mp3d->at(i).point = midsplited_point[j].point;
	        vec_mp3d->at(i).isotropic_size = midsplited_point[j].isotropic_size;
			continue;
		}
		else if (fabs(midsplited_point[j+1].riemannian_length-i*ideal_riemannian_length) < threshold_length)
		{
			vec_mp3d->at(i).point = midsplited_point[j+1].point;
	        vec_mp3d->at(i).isotropic_size = midsplited_point[j+1].isotropic_size;
			continue;
		}

		temp_ideal_riemannian_length = i * ideal_riemannian_length - midsplited_point[j].riemannian_length;
        c = interval = 0.5;
		f0 = midsplited_point[j+1].riemannian_length - midsplited_point[j].riemannian_length 
			- 0.5*(midsplited_point[j].isotropic_size - midsplited_point[j+1].isotropic_size);
        fc = f0*c+0.5*c*c*(midsplited_point[j].isotropic_size - midsplited_point[j+1].isotropic_size);
		m = 1;
		while(fabs(temp_ideal_riemannian_length-fc) > threshold_length || m <= 20)
		{
            interval *= 0.5;
			if (temp_ideal_riemannian_length < fc)
			    c -= interval; 
			else
				c += interval;
            fc = f0*c+0.5*c*c*(midsplited_point[j].isotropic_size - midsplited_point[j+1].isotropic_size);
			m++;
		}

		if (fabs(temp_ideal_riemannian_length-fc) > threshold_length && m > 20)
		{
			cout << "discretize_line_curve(): iteration exceeded!" << endl;
			return; // exit(0)? return false?
		}

		vec_mp3d->at(i).point = midsplited_point[j].point + c*(midsplited_point[j+1].point - midsplited_point[j].point);
	    vec_mp3d->at(i).isotropic_size = mesh_size_spec_->get_final_spacing(vec_mp3d->at(i).point);
	}

	all_discretized_curve_point_.push_back(vec_mp3d); // final discretized curve points result

	return;
}

void SurfaceMesher::sample_ferguson_curve_spacing(FergusonCurve *_fc, std::vector<double>& sam_spa_)
{
	assert(_fc!=NULL);

    const int cpn = _fc->get_control_point_number();
	const double cv_len = _fc->get_segment_length().at(cpn-2);
	const double min_spa = mesh_size_spec_->get_min_spacing();
	const int spn = (int)(cv_len/min_spa + 0.5) > 2 ? (int)(cv_len/min_spa + 0.5) + 1 : 3;
	const double sam_len = cv_len / (double)(spn-1);

    Vec3d vec, v0[2], v1[2], res[3];
	double u1, spa, arc_len, seg_len, tar_len;
	int i, ks, ks0;
	Vec5d coef;

	sam_spa_.assign(spn, 0.0); // actually allocate memory
    
	// first sampling point
	vec = _fc->get_control_point().at(0);
	spa = mesh_size_spec_->get_final_spacing(vec);
	sam_spa_[0] = spa;

    /** middle sampling points: first, it finds the segment of cubic
    * spline ik which contains the sampling point. st is the target
    * length, sis is the lenght of the segment and u2=st/sis is a 
    * first guess of the parameter coordinate of the sampling point. 
	*/

    ks0 = 0;
	for (i=1; i<spn-1; ++i)
	{
        arc_len = i*sam_len;
		for (ks=ks0; ks<cpn-1; ++ks)
			if (arc_len<_fc->get_segment_length().at(ks))
				break;
		assert(ks<cpn-1);

        ks0 = ks; // at the ks'th segment
		seg_len = _fc->get_segment_length().at(ks);
		tar_len = arc_len;
		if (ks>0)
		{
			seg_len -= _fc->get_segment_length().at(ks-1);
			tar_len -= _fc->get_segment_length().at(ks-1);
		}
		coef = _fc->get_segment_coefficient().at(ks);
		v0[0] = _fc->get_control_point().at(ks);
		v0[1] = _fc->get_control_point_tangent_vector().at(ks);
		v1[0] = _fc->get_control_point().at(ks+1);
		v1[1] = _fc->get_control_point_tangent_vector().at(ks+1);

		_fc->calculate_parameter_given_arc_length(coef, seg_len, 0.0, u1, tar_len, 1.0e-6);
		_fc->evaluate_ferguson_curve_segment(FergusonCurve::Zero, v0, v1, u1, res);
		spa = mesh_size_spec_->get_final_spacing(res[0]);
        sam_spa_[i] = spa;
	}

	// last sampling point
    vec = _fc->get_control_point().at(cpn-1);
	spa = mesh_size_spec_->get_final_spacing(vec);
	sam_spa_[spn-1] = spa;

	return;
}

void SurfaceMesher::calculate_discretized_arc_length_for_ferguson_curve(const std::vector<double>& _sam_spa, 
		std::vector<double>& _arc_len, const double& _sam_len, double _tol)
{
	const int spn = (int)_sam_spa.size();
	double spa0, spa1, ave_spa, tmp_arc_len, seg_elem_num, all_elem_num;
	int final_elem_num;
	double ave_elem_num, elem_num, tmp_elem_num;
	std::vector<double> sam_elem_num;
	sam_elem_num.assign(spn, 0.0);
	int i, ks0, ks;

	// integrates the distribution of spacings and calculates 
    // the number of elements to be generated
	all_elem_num = 0.0;
    for (i=0; i<spn-1; ++i)
	{
        spa0 = _sam_spa[i];
		spa1 = _sam_spa[i+1];
        ave_spa = (spa1 - spa0) / _sam_len;  
        if (ave_spa<_tol) // warning_clg: ave_spa could be less than 0.0
			seg_elem_num = 2.0*_sam_len/(spa0+spa1);
		else
            seg_elem_num = log(spa1/spa0)/ave_spa; // warning_clg: formula's meaning?
			// seg_elem_num = log(spa1/spa0)*_sam_len/(spa1-spa0);
		all_elem_num += seg_elem_num;
        sam_elem_num[i] = all_elem_num; // accumulating
	}
	// in psue, "if (itype .eq. 1)", actually it's always true
	final_elem_num = (int)(all_elem_num+0.5)>2 ? (int)(all_elem_num+0.5) : 2;
	_arc_len.assign(final_elem_num, 0.0); // allocate memory

    ave_elem_num = all_elem_num / (double)final_elem_num;

    // calculate the segment length of the final curve discretization
    ks0 = 0;
    for (i=0; i<final_elem_num-1; ++i)
	{
        elem_num = (i+1)*ave_elem_num;
        for (ks=ks0; ks<spn-1; ++ks)
            if (elem_num<sam_elem_num[ks])
			{
				spa0 = _sam_spa[ks];
				spa1 = _sam_spa[ks+1];
                ave_spa = (spa1 - spa0) / _sam_len;
                tmp_arc_len = ks * _sam_len;
				break;
			}
		assert(ks<spn-1);

		ks0 = ks;
		if (0==ks)
			tmp_elem_num = elem_num; // actually elem_num - 0
		else
			tmp_elem_num = elem_num - sam_elem_num[ks-1];
		if (ave_spa<_tol)
			tmp_arc_len += tmp_elem_num*0.5*(spa0+spa1);
		else
			tmp_arc_len += spa0*(exp(ave_spa*tmp_elem_num)-1.0)/ave_spa; // warning_clg: formula's meaning?

		_arc_len[i] = tmp_arc_len;
	}
	_arc_len[final_elem_num-1] = _sam_len*(spn-1); // entire arc length of the ferguson curve

	return;
}

void SurfaceMesher::generate_final_discrete_point_for_ferguson_curve(FergusonCurve *_fc, 
		const std::vector<double>& _arc_len, std::vector<MeshPoint3D> *_vec_mp3d, double _tol)
{
    assert(NULL!=_fc && NULL!=_vec_mp3d);

	const int dpn = (int)_arc_len.size() + 1;
	const int cpn = _fc->get_control_point_number();
    int i, ks0, ks;
	Vec3d pt, v0[2], v1[2], res[3];
	Vec5d coef;
	double u1, spa, arc_len, seg_len, tar_len;

	// the first point
	pt = _fc->get_control_point().at(0);
	spa = mesh_size_spec_->get_final_spacing(pt);
	_vec_mp3d->at(0).point = pt;
	_vec_mp3d->at(0).isotropic_size = spa;

	// the middle points
	ks0 = 0;
    for (i=1; i<dpn-1; ++i)
	{
        arc_len = _arc_len[i-1];
        for (ks=ks0; ks<cpn-2; ++ks)
			if (arc_len<_fc->get_segment_length().at(ks))
                break;

		ks0 = ks;
        seg_len = _fc->get_segment_length().at(ks);
		tar_len = arc_len;
		if (ks>0)
		{
			seg_len -= _fc->get_segment_length().at(ks-1);
			tar_len -= _fc->get_segment_length().at(ks-1);
		}
		coef = _fc->get_segment_coefficient().at(ks);
		v0[0] = _fc->get_control_point().at(ks);
		v0[1] = _fc->get_control_point_tangent_vector().at(ks);
		v1[0] = _fc->get_control_point().at(ks+1);
		v1[1] = _fc->get_control_point_tangent_vector().at(ks+1);

		_fc->calculate_parameter_given_arc_length(coef, seg_len, 0.0, u1, tar_len, 1.0e-6);
		_fc->evaluate_ferguson_curve_segment(FergusonCurve::Zero, v0, v1, u1, res);
		spa = mesh_size_spec_->get_final_spacing(res[0]);
		_vec_mp3d->at(i).point = res[0];
	    _vec_mp3d->at(i).isotropic_size = spa;
	}

	// the last point
    pt = _fc->get_control_point().at(cpn-1);
	spa = mesh_size_spec_->get_final_spacing(pt);
	_vec_mp3d->at(dpn-1).point = pt;
	_vec_mp3d->at(dpn-1).isotropic_size = spa;

	return;
}

void SurfaceMesher::discretize_ferguson_curve_with_psue(FergusonCurve *_fc, double _tol)
{
	assert(NULL!=_fc);

	const int cpn = _fc->get_control_point_number();
	const double cv_len = _fc->get_segment_length().at(cpn-2); // out of range
	const double min_spa = mesh_size_spec_->get_min_spacing();
	const int spn = (int)(cv_len/min_spa + 0.5) > 2 ? (int)(cv_len/min_spa + 0.5) + 1 : 3;
	const double sam_len = cv_len / (double)(spn-1);

	// sample to generate intermediate points and their spacing
	std::vector<double> sam_spa;	
	sample_ferguson_curve_spacing(_fc, sam_spa);
	assert(sam_spa.size()==spn);

	// calculate the arc length for each finally discretized point
	std::vector<double> arc_len;
    calculate_discretized_arc_length_for_ferguson_curve(sam_spa, 
		arc_len, sam_len, 1.0e-5);
	const int dpn = (int)arc_len.size() + 1;

	// generate the final discrete points on the curve
	std::vector<MeshPoint3D> *vec_mp3d = new std::vector<MeshPoint3D>(dpn);
    generate_final_discrete_point_for_ferguson_curve(_fc, arc_len, vec_mp3d, _tol);

	all_discretized_curve_point_.push_back(vec_mp3d);

	cout << "    No. of control points: " << cpn << "  length: " << cv_len << endl;
	cout << "    No. sampling points: " << spn << "  interval: " << sam_len << endl;
	cout << "    No. of generated points: " << dpn << endl;
	
	return;
}

void SurfaceMesher::discretize_ferguson_curve_with_linearization(FergusonCurve *_fc, double _tol)
{

	return;
}

void SurfaceMesher::discretize_one_curve(Curve *_cv, double _tol)
{
	CurveType curve_type = _cv->get_curve_type();
	double arc_length = 0.0;

	switch (curve_type)
	{
	case CurveLine: 
		discretize_line_curve((LineCurve*)_cv, _tol);
		break;
	case CurveFerguson:
		discretize_ferguson_curve_with_psue((FergusonCurve*)_cv, _tol);
		break;
	case CurveNurbs:
		// todo_clg
		break;
	default:
		// todo_clg
		break;
	}

	return;
}

void SurfaceMesher::discretize_all_curve()
{
	// std::vector<Curve*> &curve_geomety = geometry_->get_curve_geometry(); then curve_geometry is invalid?!
	int i;
	for (i=0; i<(int)geometry_->get_curve_geometry().size(); ++i)
	{
		cout << "SurfaceMesher -> discretizing the " << i+1 << "th curve ..." << endl;
		discretize_one_curve(geometry_->get_curve_geometry().at(i), 1.e-6);
	}

	return;
}

void SurfaceMesher::project_discrete_bounary_point_for_linear_triangle_surface_region(SurfaceRegion *_sr, double _tol)
{
    assert(NULL!=_sr && _tol>0.0);
	assert(3==_sr->curve_number);

	int surface_id = _sr->surface_id;
	LinearTriangleSurface *plts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	Vec3d *point = plts->get_point();
	int i, j, icv, k, sz;
	bool normal_ordered;
	Vec3d pt1, pt2;

	// order the boundary curve loop
	// for edge(point[i], point[(i+1)%3])
    for (i=0; i<3; ++i)
	{
		// find the right curve with right order
        for (j=i; j<3; ++j)
		{
			icv = _sr->curve_loop[j];
			sz = (int)all_discretized_curve_point_[icv]->size();
			if (is_equal(point[i], all_discretized_curve_point_[icv]->at(0).point, _tol*_tol)
				&& is_equal(point[(i+1)%3], all_discretized_curve_point_[icv]->at(sz-1).point, _tol*_tol))
			{
				k = j;
				normal_ordered = true;
				break;
			}
			else if (is_equal(point[i], all_discretized_curve_point_[icv]->at(sz-1).point, _tol*_tol)
				&& is_equal(point[(i+1)%3], all_discretized_curve_point_[icv]->at(0).point, _tol*_tol))
			{
                k = j;
				normal_ordered = false;
				break;
			}
		}
		assert(3!=j);

        // adjust the curve index
		// because the minus wouldn't affect the "0"th curve, here we increase temporarily the 
		// curve index by 1, which will be recovered back in this function later.
		if (k==i)
		{
			_sr->curve_loop[i]++; 
			if (!normal_ordered)
				_sr->curve_loop[k] *= -1;
		}
		else
		{
			int itmp = _sr->curve_loop[i];
			_sr->curve_loop[i] = _sr->curve_loop[k];
			_sr->curve_loop[k] = itmp;

			_sr->curve_loop[i]++; 
			if (!normal_ordered)
				_sr->curve_loop[k] *= -1;
		}
	}

	// for debugging 
	for (i=0; i<3; ++i)
	{
        icv = _sr->curve_loop[i];
		sz = (int)all_discretized_curve_point_[abs(icv)-1]->size();
		if (icv>0)
			pt1 = all_discretized_curve_point_[abs(icv)-1]->at(sz-1).point;
		else
			pt1 = all_discretized_curve_point_[abs(icv)-1]->at(0).point;

		icv = _sr->curve_loop[(i+1)%3];
        sz = (int)all_discretized_curve_point_[abs(icv)-1]->size();
		if (icv>0)
			pt2 = all_discretized_curve_point_[abs(icv)-1]->at(0).point;
		else
			pt2 = all_discretized_curve_point_[abs(icv)-1]->at(sz-1).point;

		assert(is_equal(pt1, pt2, _tol*_tol));
	}

	// compute the parametric ordered boundary points
	std::vector<MeshPoint2D> *vec_mp2d = new std::vector<MeshPoint2D>;
	MeshPoint2D mp2d;
	for (i=0; i<3; ++i)
	{
		icv = _sr->curve_loop[i];
		sz = (int)all_discretized_curve_point_[abs(icv)-1]->size();
		if (icv>0)
		{
            for (j=0; j<sz-1; ++j)
			{
				// this segment of code is the same as that of below, can they share as one?
				mp2d.isotropic_size_3d = all_discretized_curve_point_[abs(icv)-1]->at(j).isotropic_size;
				pt1 = all_discretized_curve_point_[abs(icv)-1]->at(j).point;
				plts->calculate_projection(pt1, mp2d.point);
				vec_mp2d->push_back(mp2d); // pass value
			}
		}
		else
		{
			for (j=sz-1; j>0; --j)
			{
				mp2d.isotropic_size_3d = all_discretized_curve_point_[abs(icv)-1]->at(j).isotropic_size;
				pt1 = all_discretized_curve_point_[abs(icv)-1]->at(j).point;
				plts->calculate_projection(pt1, mp2d.point);
				vec_mp2d->push_back(mp2d);
			}
		}
	}

	// only one loop, need to test, need to clear after generating mesh in this surface region
	one_initial_front_param_point_.push_back(vec_mp2d); 

	// todo_clg (20080114-0957): here we don't recover the curve index of loop
	// which will be used when getting the global mesh (global indices of vertices)
	// for (i=0; i<3; ++i)
		// _sr->curve_loop[i] = abs(_sr->curve_loop[i]) - 1;

	// output the result only for debugging, it will be deleted eventually!

	return;
}

void SurfaceMesher::boundary_value(SurfaceRegion *_sr, Vec2d _box[2], Vec2d _up[4], Vec2d _us[4], 
		Vec2i _kp[4], Vec2d _xd[4], Vec2d _bl[4])
{
    int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	int nu = fs->get_nu(), nv = fs->get_nv();

	int i;
	int cor_index[4];
	cor_index[0] = 0;
	cor_index[1] = nu - 1;
	cor_index[2] = nu*nv - 1;
	cor_index[3] = (nv-1)*nu;

	_up[0]  = Vec2d(0.0, 0.0);
	_up[1]  = Vec2d(nu-1, 0.0);
	_up[2]  = Vec2d(nu-1, nv-1);
	_up[3]  = Vec2d(0.0, nv-1);
	_box[0] = _up[0];
	_box[1] = _up[2];
	for (i=0; i<4; ++i)
		_us[i] = 0.5*(_up[i] + _up[(i+1)%4]);
	for (i=0; i<4; ++i)
	{
		_kp[i][0] = cor_index[i];
		_kp[i][1] = cor_index[(i+1)%4];
	}
	_xd[0] = Vec2d(1.0, 0.0);
	_xd[1] = Vec2d(0.0, 1.0);
	_xd[2] = Vec2d(-1.0, 0.0);
	_xd[3] = Vec2d(0.0, -1.0);
	_bl[0] = Vec2d(0.0, nu-1);
	_bl[1] = Vec2d(0.0, nv-1);
	_bl[2] = Vec2d(0.0, nu-1);
	_bl[3] = Vec2d(0.0, nv-1);

	return;
}

void SurfaceMesher::localm(SurfaceRegion *_sr, Vec3d _cur_pt, Vec2d _uv0, Vec2d _dir, 
		double _ax, double _bx, double _tol, double &_xmin, double &_fmin)
{
    int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	int nu = fs->get_nu(), nv = fs->get_nv();

    const int max_iter = 100; 
	const double golden_sec = 0.3819660, zero_eps = 1.0e-10; // the same as t in PSUE
    int iter;
    double a=_ax, b=_bx, v, w, x, u, fv, fw, fx, fu;
    double xm, d=0.0, e=0.0, etemp;
	double p, q, r; 
	double tol1, tol2;
	Vec2d uv;
	Vec3d eva_pt;

	// what exactly is the relationship between (a, b) and (v, w), are they equivalent?
	// this is the main point of this algorithm

	// initialization
	v = w = x = a + golden_sec*(b-a);
    uv = _uv0 + x*_dir;
    fs->evaluate(uv, eva_pt);
	fv = fw = fx = (_cur_pt-eva_pt).sqrnorm();

	// main loop
	for (iter=0; iter<max_iter; ++iter)
	{
        // ending configuration
		xm = 0.5*(a+b);
		tol1 = _tol*fabs(x) + zero_eps;
		tol2 = 2*tol1;
		if (fabs(x-xm)<=(tol2-0.5*(b-a))) // done
		{
			_xmin = x;
			_fmin = fx;
			return;
		}

		// calculate the moving step d to get the new evaluating point u, 
		// either use parabolic step or golden section step
		if (fabs(e)<=tol1) // use golden section step
		{
			e = x > xm ? a - x : b - x;
			d = golden_sec * e;
		}
		else
		{
			r = (x-w)*(fx-fv);
			q = (x-v)*(fx-fw);
			p = (x-v)*q - (x-w)*r;
			q = 2.0*(q-r);
			if (q>0.0)
				p = -p; 
            q = fabs(q); // then p/q is the parabolic step
			etemp = e;
			e = d;

			if (fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
			{
				// invalid parabolic step, use golden section step
				e = x > xm ? a - x : b - x;
			    d = golden_sec * e;
			}
			else // take the parabolic step
			{
                d = p / q;
				u = x + d;
				if (u-a<tol2 || b-u<tol2)
				{
					if (x<xm)
                        d = tol1;
					else
						d = -tol1;
				}
			}
		}

		// get the final new evaluating point (u, fu)
		if (fabs(d)>=tol1)
			u = x + d;
		else
		{
			if (d>0)
				u = x + tol1;
			else
				u = x - tol1;
		}
		uv = _uv0 + u*_dir;
		fs->evaluate(uv, eva_pt);
	    fu = (_cur_pt-eva_pt).sqrnorm();

		// now decide what to do with our function evaluation
		if (fu<=fx)
		{
			if (u>=x)
				a = x;
			else
				b = x;
			v = w; 
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		}
		else
		{
			if (u<x)
				a = u;
			else
				b = u;

			if (fu<=fw || w==x)
			{
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu<=fv || v==x || v==w)
			{
				v = u;
				fv = fu;
			}
		}

	} // end of main loop

	// shouldn't get here
	cout << "localm: maximum iteration exceeds!" << endl;
	_xmin = x;
	_fmin = fx;
	assert(false);

	return;
}

void SurfaceMesher::fguess(SurfaceRegion *_sr, Vec3d _cur_pt, Vec2d &_uk)
{
    int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	int nu = fs->get_nu(), nv = fs->get_nv();
	std::vector<Vec3d> &cp = fs->get_control_point();

	double drm, dr2;
	int iu, iv, index[4];
	Vec3d rp;

    // initial guess: parametric coordinates (u,v)
    // of the closest support point
	drm = 1.e+30;
    for (iu=0; iu<nu-1; ++iu)
		for (iv=0; iv<nv-1; ++iv)
		{
			index[0] = iv*nu + iu;
			index[1] = index[0] + 1;
			index[2] = index[0] + nu;
			index[3] = index[2] + 1;
            rp = 0.25*(cp[index[0]]+cp[index[1]]+cp[index[2]]+cp[index[3]]);
			dr2 = (rp - _cur_pt).sqrnorm();
			if (dr2<drm)
			{
				drm = dr2;
				_uk = Vec2d(iu+0.5, iv+0.5);
			}
		}

	return;
}

void SurfaceMesher::locu2(FergusonSurface *_fs, Vec3d _cur_pt, Vec2d &_uk, Vec3d &_eva_pt, double &_dis, int &_flag)
{
    const int nu = _fs->get_nu(), nv = _fs->get_nv();
	const int MAX_ITER = 1000;
	const int MBND = 5;
	const double eps3d = 1.e-10, eps2d = 1.e-8, eps = 1.e-10; // todo_clg: be aware of eps3d, it's not the same as PSUE
	Vec2d xuvt[(2*MBND+1)*(2*MBND+1)]; // ns=3, originally wrong: Vec2d xuvt[7]!
	const int nbnd = 3;
	int iu, iv, iter, iter1, ntr, it;
	Vec2d uv, uvr, uvk, uvc;
	double du, dis1, disk;

    _flag = 0; // 0: success; 1: iterations exceeded; -1: failure

    _dis = 1.e+36;
	for (iu=0; iu<nu; ++iu)
	{
		uv[0] = iu;
		for (iv=0; iv<nv; ++iv)
		{
            uv[1] = iv;
			_fs->evaluate(uv, _eva_pt);
			dis1 = (_eva_pt - _cur_pt).sqrnorm();
			if (dis1<_dis)
			{
				uvr = uv;
				_dis = dis1;

				if (_dis<eps3d)
				{
					_uk = uvr;
					return;
				}
			}
		}
	}

	uv = Vec2d(uv[0]+eps, uv[1]+eps); // originally wrong: "Vec2d(uv[0]-1.0+eps, uv[1]-1.0+eps)"

	du = 1.0;
	for (iter=0; iter<MAX_ITER; ++iter)
	{
        du *= 0.5;
		// nbnd = 3;

		for (iter1=0; true; ++iter1)
		{
			uvk = uvr;
			disk = _dis;
			ntr = 0;
			for (uvc[0]=uvr[0]-du*nbnd, iu=0; iu<2*nbnd+1; ++iu, uvc[0]+=du)
			    for (uvc[1]=uvr[1]-du*nbnd, iv=0; iv<2*nbnd+1; ++iv, uvc[1]+=du)
				    if (uvc[0]>-eps && uvc[1]>-eps && uvc[0]<uv[0] && uvc[1]<uv[1])
					{
						xuvt[ntr] = uvc;
						++ntr;
					}
			
            for (it=0; it<ntr; ++it)
			{
                _fs->evaluate(xuvt[it], _eva_pt);
				dis1 = (_eva_pt - _cur_pt).sqrnorm();
			    if (dis1<_dis)
			    {
				    uvr = xuvt[it];
				    _dis = dis1;

				    if (_dis<eps3d)
				    {
					    _uk = uvr;
					    return;
				    }
			    }
			}

			if (uvr[0]>uvk[0]-du*nbnd+eps && uvr[0]<uvk[0]+du*nbnd-eps
				&& uvr[1]>uvk[1]-du*nbnd+eps && uvr[1]<uvk[1]+du*nbnd-eps)
			    break;
		}

		if (_dis>disk)
		{
			_flag = -1;
			return;
		}
		else
			disk = _dis;

		if (du<=eps2d)
			break;
	}

	if (iter>=MAX_ITER)
	{
		_flag = 1;
		return;
	}

	_uk = uvr;
	
	return;
}

bool SurfaceMesher::feasib(Vec2d _uk, Vec2d _dr, Vec2d _box[2], double _tol)
{
	Vec2d dn = _dr / _dr.norm();
	Vec2d up = _uk + _tol*dn; // New position

    // still in the parameter plane ?
	if (up[0]>=_box[0][0] && up[0]<=_box[1][0] && up[1]>=_box[0][1] && up[1]<=_box[1][1])
		return true;
	else
		return false;
}

bool SurfaceMesher::newdr(Vec2d _uk, Vec2d _gk, Vec2d &_dr, Vec2d _box[2], double _tol)
{
	const double eps = 1.e-5;
	bool flag = true;

	// Checks the active constraints
	double dub = fabs(_uk[0] - _box[0][0]), 
		   dut = fabs(_uk[0] - _box[1][0]),
		   dvb = fabs(_uk[1] - _box[0][1]),
		   dvt = fabs(_uk[1] - _box[1][1]);

	/** Projects the descent direction ( opposite to the gradient )
    * along the active boundary and verifies whether the new 
    * direction is acceptable
	*/
	if (dub<_tol || dut<_tol)
		if (fabs(_gk[0])>eps)
		{
			if (-_gk[0]>0.0)
				_dr[0] = 1.0;
			else
				_dr[0] = -1.0;
			_dr[1] = 0.0;

			flag = feasib(_uk, _dr, _box, _tol);
			if (flag)
				return true;
		}

	if (dvb<_tol || dvt<_tol)
		if (fabs(_gk[1])>eps)
		{
			_dr[0] = 0.0;
			if (-_gk[1]>0.0)
				_dr[1] = 1.0;
			else
				_dr[1] = -1.0;

			flag = feasib(_uk, _dr, _box, _tol);
			if (flag) 
                return true;
		}

	return false;
}

void SurfaceMesher::sort(double _a[], int _sz)
{
	int i, j, k;
	double tmp;

	for (i=0; i<_sz-1; ++i)
	{
		k = i;
		for (j=i+1; j<_sz; ++j)
			if (_a[j]<_a[k])
				k = j;

		tmp = _a[k];
		_a[k] = _a[i];
		_a[i] = tmp;
	}

	return;
}

void SurfaceMesher::project_one_loop_for_ferguson_surface_region(SurfaceRegion *_sr, 
		std::vector<int> &_loop, double _tol)
{
    assert(NULL!=_sr && _tol>0.0);

	/**
    Tolerances for the iterative process
    EPS1 ...... distance in 3D coordinates
    EPS2 ...... distance in (u,v) coordinates
    EPS3 ...... modulus of search direction
    EPS4 ...... parallel directions ( 6 deg )
	*/
	const int nit = 100;
	const double tol = 1.0e-5;
    const double tolg = 1.0e-4; // todo_clg: not the same as PSUE
	const double eps1 = tolg, eps = eps1*eps1, // todo_clg: eps = 1.0e-7 would cause invalid elements
		eps2 = 1.0e-4, eps3 = 1.0e-10, eps4 = 0.1;

    int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	int nu = fs->get_nu(), nv = fs->get_nv();

    std::vector<MeshPoint2D> *vec_mp2d = new std::vector<MeshPoint2D>;
	MeshPoint2D mp2d;
	Vec3d cur_pt, eva_pt, diff_vec, su, sv;
	Vec2d uk, gk, dr, gk1;
	int i, j, k, it, icv, sz, kb, istck, iflag;
	double xmin, fmin, sqr_len, len, dist, vsc, dmd, vmd, bet;
	bool flag;
	double u1, v1, w[4], ax, bx;

	// calculate some boundary values, the same as subroutine bval in PSUE
	// xd, bl, integer kp(2,4), for what?
	Vec2d box[2], us[4], up[4], xd[4], bl[4]; 
	Vec2i kp[4];
    boundary_value(_sr, box, up, us, kp, xd, bl);

	for (i=0; i<_loop.size(); ++i)
	{
		icv = abs(_loop[i]) - 1;
		sz = all_discretized_curve_point_[icv]->size();
		if (_loop[i]>0)
		{
            for (j=0; j<sz-1; ++j)
		    {
				mp2d.isotropic_size_3d = all_discretized_curve_point_[icv]->at(j).isotropic_size;
				cur_pt = all_discretized_curve_point_[icv]->at(j).point;
				
				// First, checks for singular points
                for (k=0; k<4; ++k)
				{
					fs->evaluate(us[k], eva_pt);
					sqr_len = (eva_pt-cur_pt).sqrnorm();
                    if (sqr_len<eps)
					{
						mp2d.point = us[k];
						len = sqrt(sqr_len);
						goto finished0;
					}
				}

				// Finds the minimum along the boundary
                sqr_len = 1.0e+30;
                kb = -1;
                for (k=0; k<4; ++k)
				{
					dist = (fs->get_control_point().at(kp[k][0]) - fs->get_control_point().at(kp[k][1])).sqrnorm();
					if (dist>=eps) // not degenerated, todo_clg: be aware of eps
					{
                        localm(_sr, cur_pt, up[k], xd[k], bl[k][0], bl[k][1], tol, xmin, fmin);
                        
						// record the closest edge
						if (fmin<sqr_len)
						{
							sqr_len = fmin;
							kb = k;
                            uk = up[k] + xmin*xd[k];
						}

						// find the final projected point, ending configuration
						if (sqr_len<eps)
						{
							mp2d.point = uk;
							len = sqrt(sqr_len);
							goto finished0;
						}
					}
				}

				/** If the minimum has not been found on the boundary, the point 
                * with minimum distance is used as the initial guess in the 
                * iterative process. But first we check that the first search 
                * direction does not coincide with the current one on the 
                * boundary. If it does we select an alternative first guess.
                */
                
				fs->evaluate(uk, eva_pt);
				fs->calculate_first_derivative(uk, su, sv);
				diff_vec = eva_pt - cur_pt;
				sqr_len = diff_vec.sqrnorm();
				if (sqr_len<eps)
				{
					mp2d.point = uk;
					len = sqrt(sqr_len);
					goto finished0;
				}
				gk = 2.0*Vec2d(diff_vec|su, diff_vec|sv);
				vmd = gk.norm();
				if (vmd<eps3)
				{
					mp2d.point = uk;
					len = sqrt(sqr_len);
					goto finished0;
				}
				vmd = 1. / vmd;
				dr = -gk / gk.norm(); 

				// If the direction is not valid uses a different initial 
                // guess: the closest patch center point
				vsc = fabs(dr[0]*xd[kb][1]-dr[1]*xd[kb][0]);
                if (vsc<eps4)
				{
                    fguess(_sr, cur_pt, uk);
					fs->evaluate(uk, eva_pt);
					fs->calculate_first_derivative(uk, su, sv);
					diff_vec = eva_pt - cur_pt;
					sqr_len = diff_vec.sqrnorm();
					if (sqr_len<eps)
					{
						mp2d.point = uk;
						len = sqrt(sqr_len);
						goto finished0;
					}
					gk = 2.0*Vec2d(diff_vec|su, diff_vec|sv);
					dr = -gk / gk.norm();
				}

                // Here tries a more robust but slower algorithm 
                // when at a boundary point the line search has taken the
                // iterative procedure out-of-bounds
                istck = 0;
tryagain0:
				if (0!=istck)
				{
					cout << "Target point .. " << "i: " << i << " _loop[i]: " << _loop[i] << " j: " << j << endl;
					cout << "Coordinates ... " << cur_pt[0] << "  " << cur_pt[1] << "  " << cur_pt[2] << endl;
					cout << "Last .......... " << eva_pt[0] << "  " << eva_pt[1] << "  " << eva_pt[2] << endl;
					cout << "Local ......... " << uk[0] << "  " << uk[1] << endl;
					cout << "Distance ...... " << sqrt(sqr_len) << endl;

                    locu2(fs, cur_pt, uk, eva_pt, dist, iflag);                    
                    
                    if (0==iflag)
					    cout << "LOCUV > Successful second attempt" << endl;
					else if (1==iflag)
						cout << "LOCUV > 2nd attempt: number of iteration exceeded" << endl;
					else if (-1==iflag)
                        cout << "LOCUV > Second attempt failed" << endl;
 
                    len = sqrt(dist);
					cout << "LOCUV > Values after second attempt" << endl;
					cout << "Target Point .. " << "i: " << i << " _loop[i]: " << _loop[i] << " j: " << j << endl;
					cout << "Coordinates ... " << cur_pt[0] << "  " << cur_pt[1] << "  " << cur_pt[2] << endl;
					cout << "Last .......... " << eva_pt[0] << "  " << eva_pt[1] << "  " << eva_pt[2] << endl;
					cout << "Local ......... " << uk[0] << "  " << uk[1] << endl;
					cout << "Distance ...... " << len << endl;

					if (len>tolg)
					{
						cout << "LOCUV >        !!! *** WARNING *** !!!" << endl;
						cout << "The distance between the target point and" << endl;
						cout << "and the last point in the iteration is: " << len << endl;
						cout << "This is BIGGER than the tolerance used by" << endl;
						cout << "SURFACE for coincident points which is: " << endl;
						cout << "The run continues but you are advised to " << endl;
						cout << "check your input geometry data for errors" << endl;
					}

					mp2d.point = uk;
					goto finished0;
				}
				// end of tryagain0

                // iteration for finding the local coordinates
                for (it=0; it<nit; ++it)
				{
					// ending configuration
					dmd = dr.norm();
					if (dmd<eps3) 
					{
						mp2d.point = uk;
						goto finished0;
					}
					dmd = 1. / dmd;
					dr.normalize();

					// is the current direction valid?
                    flag = feasib(uk, dr, box, tol);
					if (!flag)
					{
                        flag = newdr(uk, gk, dr, box, tol);
                        if (!flag)
						{
							istck += 1;
							cout << "LOCUV > got stuck ... trying again" << endl;
							goto tryagain0;
						}
					} 

                    // Max. values of the line parameter in the (du,dv) direction
					u1 = fabs(dr[0])>eps3 ? fabs(dr[0]) : eps3;
					if (dr[0]<0.0)
						u1 = -u1;
					u1 = 1.0 / u1;
					v1 = fabs(dr[1])>eps3 ? fabs(dr[1]) : eps3;
					if (dr[1]<0.0)
						v1 = -v1;
					v1 = 1.0 / v1;

					w[0] = (box[1][0]-uk[0])*u1;
                    w[1] = (box[0][0]-uk[0])*u1;
					w[2] = (box[1][1]-uk[1])*v1;
					w[3] = (box[0][1]-uk[1])*v1;
                    sort(w, 4);
                    ax = 0.0;
					bx = w[2];
                    if (fabs(bx)<tol)
                        bx = w[3];

					// Brent's procedure localm for line search
					localm(_sr, cur_pt, uk, dr, ax, bx, tol, xmin, fmin);
                    uk += xmin * dr;
					fs->evaluate(uk, eva_pt);
					fs->calculate_first_derivative(uk, su, sv);
					diff_vec = eva_pt - cur_pt;
                    sqr_len = diff_vec.sqrnorm();
					if (sqr_len<eps)
					{
						mp2d.point = uk;
						len = sqrt(sqr_len);
						goto finished0;
					}
					gk1 = 2.0*Vec2d(diff_vec|su, diff_vec|sv);
					if (1==it%2)
						dr = -gk1;
					else
					{
						bet = (gk1-gk) | gk1;
						bet *= vmd;
						dr = -gk1 + bet*dr;
					}
					gk = gk1;
					vmd = gk.norm();
					if (vmd<eps)
					{
						mp2d.point = uk;
						len = sqrt(sqr_len);
						goto finished0;
					}
					vmd = 1.0 / vmd;

				} // end of "for (it=0; i<nit; ++it)"

				len = sqrt(sqr_len);
				if (len>tolg)
				{
					cout << "LOCUV > Number of iterations exceeded ... trying again" << endl;
					istck++;
					goto tryagain0;
				}


finished0: 
				vec_mp2d->push_back(mp2d);
		    }
		}
		else
		{
            // for (j=0; j<sz-1; ++j)
			for (j=sz-1; j>0; --j)
		    {
				mp2d.isotropic_size_3d = all_discretized_curve_point_[icv]->at(j).isotropic_size;
				cur_pt = all_discretized_curve_point_[icv]->at(j).point;
				
				// First, checks for singular points
                for (k=0; k<4; ++k)
				{
					fs->evaluate(us[k], eva_pt);
					sqr_len = (eva_pt-cur_pt).sqrnorm();
                    if (sqr_len<eps)
					{
						mp2d.point = us[k];
						len = sqrt(sqr_len);
						goto finished1;
					}
				}

				// Finds the minimum along the boundary
                sqr_len = 1.0e+30;
                kb = -1;
                for (k=0; k<4; ++k)
				{
					dist = (fs->get_control_point().at(kp[k][0]) - fs->get_control_point().at(kp[k][1])).sqrnorm();
					if (dist>=eps) // not degenerated, todo_clg: be aware of eps
					{
                        localm(_sr, cur_pt, up[k], xd[k], bl[k][0], bl[k][1], tol, xmin, fmin);
                        
						// record the closest edge
						if (fmin<sqr_len)
						{
							sqr_len = fmin;
							kb = k;
                            uk = up[k] + xmin*xd[k];
						}

						// find the final projected point, ending configuration
						if (sqr_len<eps)
						{
							mp2d.point = uk;
							len = sqrt(sqr_len);
							goto finished1;
						}
					}
				}

				/** If the minimum has not been found on the boundary, the point 
                * with minimum distance is used as the initial guess in the 
                * iterative process. But first we check that the first search 
                * direction does not coincide with the current one on the 
                * boundary. If it does we select an alternative first guess.
                */
                
				fs->evaluate(uk, eva_pt);
				fs->calculate_first_derivative(uk, su, sv);
				diff_vec = eva_pt - cur_pt;
				sqr_len = diff_vec.sqrnorm();
				if (sqr_len<eps)
				{
					mp2d.point = uk;
					len = sqrt(sqr_len);
					goto finished1;
				}
				gk = 2.0*Vec2d(diff_vec|su, diff_vec|sv);
				vmd = gk.norm();
				if (vmd<eps3)
				{
					mp2d.point = uk;
					len = sqrt(sqr_len);
					goto finished1;
				}
				vmd = 1. / vmd;
				dr = -gk / gk.norm(); 

				// If the direction is not valid uses a different initial 
                // guess: the closest patch center point
				vsc = fabs(dr[0]*xd[kb][1]-dr[1]*xd[kb][0]);
                if (vsc<eps4)
				{
                    fguess(_sr, cur_pt, uk);
					fs->evaluate(uk, eva_pt);
					fs->calculate_first_derivative(uk, su, sv);
					diff_vec = eva_pt - cur_pt;
					sqr_len = diff_vec.sqrnorm();
					if (sqr_len<eps)
					{
						mp2d.point = uk;
						len = sqrt(sqr_len);
						goto finished1;
					}
					gk = 2.0*Vec2d(diff_vec|su, diff_vec|sv);
					dr = -gk / gk.norm();
				}

                // Here tries a more robust but slower algorithm 
                // when at a boundary point the line search has taken the
                // iterative procedure out-of-bounds
                istck = 0;
tryagain1:
				if (0!=istck)
				{
                    cout << "Target point .. " << "i: " << i << " _loop[i]: " << _loop[i] << " j: " << j << endl;
					cout << "Coordinates ... " << cur_pt[0] << "  " << cur_pt[1] << "  " << cur_pt[2] << endl;
					cout << "Last .......... " << eva_pt[0] << "  " << eva_pt[1] << "  " << eva_pt[2] << endl;
					cout << "Local ......... " << uk[0] << "  " << uk[1] << endl;
					cout << "Distance ...... " << len << endl;

                    locu2(fs, cur_pt, uk, eva_pt, dist, iflag);                    
                    
                    if (0==iflag)
					    cout << "LOCUV > Successful second attempt" << endl;
					else if (1==iflag)
						cout << "LOCUV > 2nd attempt: number of iteration exceeded" << endl;
					else if (-1==iflag)
                        cout << "LOCUV > Second attempt failed" << endl;
 
                    len = sqrt(dist);
					cout << "LOCUV > Values after second attempt" << endl;
					cout << "Target Point .. " << "i: " << i << " _loop[i]: " << _loop[i] << " j: " << j << endl;
					cout << "Coordinates ... " << cur_pt[0] << "  " << cur_pt[1] << "  " << cur_pt[2] << endl;
					cout << "Last .......... " << eva_pt[0] << "  " << eva_pt[1] << "  " << eva_pt[2] << endl;
					cout << "Local ......... " << uk[0] << "  " << uk[1] << endl;
					cout << "Distance ...... " << len << endl;

					if (len>tolg)
					{
						cout << "LOCUV >        !!! *** WARNING *** !!!" << endl;
						cout << "The distance between the target point and" << endl;
						cout << "and the last point in the iteration is: " << len << endl;
						cout << "This is BIGGER than the tolerance used by" << endl;
						cout << "SURFACE for coincident points which is: " << endl;
						cout << "The run continues but you are advised to " << endl;
						cout << "check your input geometry data for errors" << endl;
					}

					mp2d.point = uk;
					goto finished1;
				}
				// end of tryagain1

                // iteration for finding the local coordinates
                for (it=0; i<nit; ++it)
				{
					// ending configuration
					dmd = dr.norm();
					if (dmd<eps3) 
					{
						mp2d.point = uk;
						goto finished1;
					}
					dmd = 1. / dmd;
					dr.normalize();

					// is the current direction valid?
                    flag = feasib(uk, dr, box, tol);
					if (!flag)
					{
                        flag = newdr(uk, gk, dr, box, tol);
                        if (!flag)
						{
							istck += 1;
							cout << "LOCUV > got stuck ... trying again" << endl;
							goto tryagain1;
						}
					} 

                    // Max. values of the line parameter in the (du,dv) direction
					u1 = fabs(dr[0])>eps3 ? fabs(dr[0]) : eps3;
					if (dr[0]<0.0)
						u1 = -u1;
					u1 = 1.0 / u1;
					v1 = fabs(dr[1])>eps3 ? fabs(dr[1]) : eps3;
					if (dr[1]<0.0)
						v1 = -v1;
					v1 = 1.0 / v1;

					w[0] = (box[1][0]-uk[0])*u1;
                    w[1] = (box[0][0]-uk[0])*u1;
					w[2] = (box[1][1]-uk[1])*v1;
					w[3] = (box[0][1]-uk[1])*v1;
                    sort(w, 4);
                    ax = 0.0;
					bx = w[2];
                    if (fabs(bx)<tol)
                        bx = w[3];

					// Brent's procedure localm for line search
					localm(_sr, cur_pt, uk, dr, ax, bx, tol, xmin, fmin);
                    uk += xmin * dr;
					fs->evaluate(uk, eva_pt);
					fs->calculate_first_derivative(uk, su, sv);
					diff_vec = eva_pt - cur_pt;
                    sqr_len = diff_vec.sqrnorm();
					if (sqr_len<eps)
					{
						mp2d.point = uk;
						len = sqrt(sqr_len);
						goto finished1;
					}
					gk1 = 2.0*Vec2d(diff_vec|su, diff_vec|sv);
					if (0==it%2)
						dr = -gk1;
					else
					{
						bet = (gk1-gk) | gk1;
						bet *= vmd;
						dr = -gk1 + bet*dr;
					}
					gk = gk1;
					vmd = gk.norm();
					if (vmd<eps)
					{
						mp2d.point = uk;
						len = sqrt(sqr_len);
						goto finished1;
					}
					vmd = 1.0 / vmd;

				} // end of "for (it=0; i<nit; ++it)"

				len = sqrt(sqr_len);
				if (len>tolg)
				{
					cout << "LOCUV > Number of iterations exceeded ... trying again" << endl;
					istck++;
					goto tryagain1;
				}


finished1: 
				vec_mp2d->push_back(mp2d);
			}
		}
	}

	one_initial_front_param_point_.push_back(vec_mp2d);
}

void SurfaceMesher::orientate_initial_parametric_front(SurfaceRegion *_sr)
{
	int surface_id = _sr->surface_id, curve_number = _sr->curve_number, *curve_loop = _sr->curve_loop;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);

	const int nu = fs->get_nu(), nv = fs->get_nv();
    const double EPS1 = 1.E-6, ARMX = (nu-1+EPS1)*(nv-1+EPS1);
	std::vector<MeshPoint2D> *vec_mp2d = NULL;
	std::vector<double> area;
	std::vector<int> *loop = NULL;
	double ar, tmp_ar;
	int i, j, k, nseg, np, nloop = one_initial_front_param_point_.size();

	// compute signed area of each loop
	for (i=0; i<nloop; ++i)
	{
		ar = 0;
        vec_mp2d = one_initial_front_param_point_[i];
		np = vec_mp2d->size();
		for (j=0; j<np-1; ++j)
		{
			ar += (vec_mp2d->at(j).point[0]+vec_mp2d->at(j+1).point[0])
				*(vec_mp2d->at(j+1).point[1]-vec_mp2d->at(j).point[1]);
		}
		ar += (vec_mp2d->at(np-1).point[0]+vec_mp2d->at(0).point[0])
			*(vec_mp2d->at(0).point[1]-vec_mp2d->at(np-1).point[1]);
		ar *= 0.5;
		
		if (fabs(ar)>ARMX)
		{
			cout << "SurfaceMesher-err > Area error, Loop of segments has an invalid area." << endl;
			assert(false);
		}
		
		area.push_back(ar);
	}

    // Sorts the loops according to the area
	// Use selective sort, from biggest to smallest
	for (i=0; i<nloop-1; ++i)
	{
        k = i;
		for (j=i+1; j<nloop; ++j)
			if (fabs(area[j])>fabs(area[k]))
				k = j;
		if (k!=i)
		{
            tmp_ar = area[i];
			area[i] = area[k];
			area[k] = tmp_ar;

			loop = loops[i];
			loops[i] = loops[k];
			loops[k] = loop;

			vec_mp2d = one_initial_front_param_point_[i];
			one_initial_front_param_point_[i] = one_initial_front_param_point_[k];
			one_initial_front_param_point_[k] = vec_mp2d;
		}
	}

	// Modifies the orientation of the loops and the sides in the front
    // according to the sign of the computed areas.
	MeshPoint2D tmp_mp2d;
    for (i=0; i<nloop; ++i)
	{
		if (0==i)
			ar = area[i];
		else
			ar = -area[i];

		if (ar<0) // opposite orientation
		{
			loop = loops[i];
			nseg = loop->size();
			for (j=0; j<=(int)((nseg-1)/2); ++j)
			{
				k = loop->at(j);
				loop->at(j) = -loop->at(nseg-1-j);
				loop->at(nseg-1-j) = -k;
			}

			vec_mp2d = one_initial_front_param_point_[i];
			np = vec_mp2d->size();
			for (j=1; j<=(int)((np-1)/2); ++j) // attention: j doesn't start from 0
			{
				tmp_mp2d.isotropic_size_3d = vec_mp2d->at(j).isotropic_size_3d;
				tmp_mp2d.metric_tensor = vec_mp2d->at(j).metric_tensor;
				tmp_mp2d.point = vec_mp2d->at(j).point;
				vec_mp2d->at(j).isotropic_size_3d = vec_mp2d->at(np-j).isotropic_size_3d;
				vec_mp2d->at(j).metric_tensor = vec_mp2d->at(np-j).metric_tensor;
				vec_mp2d->at(j).point = vec_mp2d->at(np-j).point;
                vec_mp2d->at(np-j).isotropic_size_3d = tmp_mp2d.isotropic_size_3d;
				vec_mp2d->at(np-j).metric_tensor = tmp_mp2d.metric_tensor;
                vec_mp2d->at(np-j).point = tmp_mp2d.point;
			}
		}
	}

	k = 0;
	for (i=0; i<loops.size(); ++i)
		for (j=0; j<loops[i]->size(); ++j)
			curve_loop[k++] = loops[i]->at(j);

	return;
}

void SurfaceMesher::project_discrete_bounary_point_for_ferguson_surface_region(SurfaceRegion *_sr, double _tol)
{
    assert(NULL!=_sr && _tol>0.0);

	int surface_id = _sr->surface_id, curve_number = _sr->curve_number, *curve_loop = _sr->curve_loop;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	int i, j, k, icv, sz, itemp;
	bool normal_ordered;
	Vec3d pt0, pt1, pt2, pt3;

	// std::vector< std::vector<int>* > loops; make it as class member
	for (i=0; i<loops.size(); ++i)
	{
		delete loops[i];
	}
	loops.clear();

	std::vector<int> *loop = NULL;

	// for (i=0; i<one_initial_front_param_point_.size(); ++i)
	// {
	//     one_initial_front_param_point_[i]->clear();
	//	   delete one_initial_front_param_point_[i];
	// }
	// one_initial_front_param_point_.clear();

    // order the boundary curve loop
	// first increase the curve index by 1 for minus does not affect 0, which will be recovered later
	// in this function. the sign of the curve index will be used to determine the order of the curve 
	// in this loop. we do it only because the "0" curve. it should be improved in the future.

	// the above is wrong, curve_loop[0]==1, thus no need to have the following segment of code
    // for (i=0; i<curve_number; ++i)
		// ++curve_loop[i];

	loop = new std::vector<int>;
    for (i=0; i<curve_number; ++i)
	{
		if (0==loop->size())
		{
			loop->push_back(curve_loop[i]);
			sz = all_discretized_curve_point_.size(); // for debug
			pt0 = all_discretized_curve_point_[curve_loop[i]-1]->front().point;

			pt1 = all_discretized_curve_point_[curve_loop[i]-1]->back().point;
			if (is_equal(pt0, pt1, 1.0e-12)) // one loop only made up of one curve 
			{
				loops.push_back(loop);
				loop = new std::vector<int>;
			}

			continue;
		}

		icv = loop->at(loop->size()-1);
		if (icv>0)
			pt1 = all_discretized_curve_point_[abs(icv)-1]->back().point;
		else
			pt1 = all_discretized_curve_point_[abs(icv)-1]->front().point;

		for (j=i; j<curve_number; ++j)
		{
            pt2 = all_discretized_curve_point_[curve_loop[j]-1]->front().point;
		    pt3 = all_discretized_curve_point_[curve_loop[j]-1]->back().point;
		    if (is_equal(pt1, pt2, 1.0e-5)) // todoclg: should not be 1.0e-8
		    {
			    loop->push_back(curve_loop[j]);
				if (i!=j)
				{
                    itemp = curve_loop[i];
					curve_loop[i] = curve_loop[j];
					curve_loop[j] = itemp;
				}
                if (is_equal(pt0, pt3, 1.0e-5))
			    {
				    loops.push_back(loop);
				    loop = new std::vector<int>;
			    }

				break;
		    }
		    else if (is_equal(pt1, pt3, 1.0e-5)) // todoclg: should not be 1.0e-8
		    {
			    loop->push_back(-curve_loop[j]);
				if (i!=j)
				{
                    itemp = curve_loop[i];
					curve_loop[i] = curve_loop[j];
					curve_loop[j] = itemp;
				}
			    if (is_equal(pt0, pt2, 1.0e-5))
			    {
				    loops.push_back(loop);
				    loop = new std::vector<int>;
			    }

				break;
		    }
		    else
			    continue;
		}
		assert(j!=curve_number);
	}

	// project the discretized points of each loop
	for (i=0; i<loops.size(); ++i)
        project_one_loop_for_ferguson_surface_region(_sr, *loops[i], _tol);
	
	orientate_initial_parametric_front(_sr);

	return;
}

void SurfaceMesher::project_discrete_boundary_point(SurfaceRegion *_sr, double _tol)
{
    assert(NULL!=_sr && _tol>0.0);

	SurfaceType surface_type = geometry_->get_surface_geometry().at(_sr->surface_id-1)->get_surface_type(); // originally "_sr->surface_id-1"
	switch(surface_type)
	{
    case SurfaceLinearTriangle:
		project_discrete_bounary_point_for_linear_triangle_surface_region(_sr, _tol);
		break;
	case SurfaceQuadraticTriangle:
		// todo_clg
		break;
	case SurfaceFerguson:
		project_discrete_bounary_point_for_ferguson_surface_region(_sr, _tol);
		break;
	case SurfaceNurbs:
		// todo_clg
		break;
	default:
		;
	}

	return;
}

void SurfaceMesher::generate_boundary_mesh_for_linear_triangle_surface_region(SurfaceRegion *_sr, 
	    SurfaceMesh *_sm, double _tol)
{
    assert(NULL!=_sr && NULL!=_sm && _tol>0.0);

	// int surface_id = _sr->surface_id;
	// LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id);

	// generate the boundary mesh
	int i;
	DTIso2D generator;
	assert(generator.crtEnv());
	generator.initBouInfo(one_initial_front_param_point_);
	generator.bndPntInst();
	generator.recoverBnds();
	generator.clrOuterEles();
	generator.recDistPnts();

	// for (i=0; i<4; ++i)
	// 	generator.innerPntInst();
	// generator.smooth();

	generator.rmvEmpNods();
	generator.rmvEmpEles();	
	generator.updateBndParent();

	generator.output();
	generator.writePl2("boundary_param_mesh.pl2");
	generator.writeOff("boundary_param_mesh.off");

	// clear one_initial_front_param_point_
	for (i=0; i<one_initial_front_param_point_.size(); ++i)
	{
		one_initial_front_param_point_[i]->clear();
        delete one_initial_front_param_point_[i];
	}
	one_initial_front_param_point_.clear();

	// write into _sm;

	return;
}



// functions below are kind of messed up
void SurfaceMesher::riemannian_length_2d_line(SurfaceRegion *_sr, Vec3d& _pt1, Vec3d& _pt2, double &_rl)
{
    int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	Vec3d ru, rv, pt;
	lts->calculate_first_derivative(Vec2d(0.0, 0.0), ru, rv); // can get 2d metric tensor from ru, rv
	Mat2by2d metric_tensor(ru | ru, ru | rv, ru | rv, rv | rv);
	Mat2by1d vec(_pt2[1]-_pt1[1], _pt2[0]-_pt1[0]);
	Mat1by1d cons = ~vec*metric_tensor *vec;

    _rl = sqrt(cons(0,0));

	double h1, h2, h3, h;
	Vec2d uv(_pt1[0], _pt1[1]);
	lts->evaluate(uv, pt);
	h1 = mesh_size_spec_->get_final_spacing(pt);
	uv[0] = 0.5*(_pt1[0]+_pt2[0]);
	uv[1] = 0.5*(_pt1[1]+_pt2[1]);
    lts->evaluate(uv, pt);
	h2 = mesh_size_spec_->get_final_spacing(pt);
	uv[0] = _pt2[0];
	uv[1] = _pt2[1];
    lts->evaluate(uv, pt); 
	h3 = mesh_size_spec_->get_final_spacing(pt);
    h = (1.0/h1 + 4.0/h2 + 1.0/h3) / 6;

	_rl *= h;

	return;
}

void SurfaceMesher::create_binary_tree_2(SurfaceRegion *_sr, MidsplitedEdge *_mse)
{
	const double threshold_length = 0.1; // _tol!
	if (_mse->riemannian_length <= threshold_length)  // recursive termination condition
		return;

	_mse->left = new MidsplitedEdge;
	_mse->right = new MidsplitedEdge;

	_mse->left->point[0] = _mse->point[0];
	_mse->left->point[1] = (_mse->point[0]+_mse->point[1])/2.0;
	_mse->left->left = _mse->left->right = NULL;
	_mse->left->parent = _mse;
	riemannian_length_2d_line(_sr, _mse->left->point[0], _mse->left->point[1], _mse->left->riemannian_length);

    _mse->right->point[0] = (_mse->point[0]+_mse->point[1])/2.0;
	_mse->right->point[1] = _mse->point[1];
	_mse->right->left = _mse->right->right = NULL;
	_mse->right->parent = _mse;
	riemannian_length_2d_line(_sr, _mse->right->point[0], _mse->right->point[1], _mse->right->riemannian_length);

	create_binary_tree_2(_sr, _mse->left);
	create_binary_tree_2(_sr, _mse->right);
}

void SurfaceMesher::traverse_to_get_midsplited_point_2(MidsplitedEdge *_mse, 
	    std::vector<MidsplitedPoint> &_midsplited_point)
{
    assert(NULL!=_mse);

	// leaf node, recursive termination condition
	if (NULL==_mse->left && NULL==_mse->right) 
	{
        MidsplitedPoint msp;
	    msp.point = _mse->point[1];
		msp.riemannian_length = _mse->riemannian_length;
		_midsplited_point.push_back(msp);

		return;
	}
	else
	{
		traverse_to_get_midsplited_point_2(_mse->left, _midsplited_point);
		traverse_to_get_midsplited_point_2(_mse->right, _midsplited_point);
	}

	// free memory!!! do it later!

	return;
}

void SurfaceMesher::collapse_2d_short_edges(SurfaceRegion *_sr, MyMesh &_mesh, double _tol)
{
    int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	const double length_threshold = 0.5;

	HalfedgeHandle heh0, heh1;
	VertexHandle   vh0, vh1;
	EdgeHandle     eh;
	double total_riemannian_length = 0.0;
	int i;
	MeshPoint3D mp3d;

	std::vector<MidsplitedPoint> midsplited_point;
	MidsplitedPoint msp;
	std::list<MeshPoint3D> *lst_mp3d = NULL;

	MyMesh::EdgeIter e_it;
	for (e_it=_mesh.edges_begin(); e_it!=_mesh.edges_end(); ++e_it)
	{
		eh = e_it.handle();
		if (_mesh.is_boundary(eh) || _mesh.status(eh).deleted())
			continue;
		
        heh0   = _mesh.halfedge_handle(e_it.handle(), 0);
	    heh1   = _mesh.halfedge_handle(e_it.handle(), 1);
        vh0    = _mesh.to_vertex_handle(heh0);
	    vh1    = _mesh.to_vertex_handle(heh1);

		if (_mesh.face_handle(heh0).is_valid()&&_mesh.status(_mesh.face_handle(heh0)).deleted())
			continue;
		if (_mesh.face_handle(heh1).is_valid()&&_mesh.status(_mesh.face_handle(heh1)).deleted())
			continue;

        // create binary tree to generate intermediate auxiliary points
        MidsplitedEdge *root = new MidsplitedEdge; // where to free memory?!
	    root->point[0] = _mesh.point(vh0);
	    root->point[1] = _mesh.point(vh1);
	    root->left = root->right = NULL;
	    root->parent = NULL;
	    riemannian_length_2d_line(_sr, root->point[0], root->point[1], root->riemannian_length);
        create_binary_tree_2(_sr, root); 

	    msp.point = root->point[0];
	    msp.riemannian_length = 0.0;
	    midsplited_point.push_back(msp);
	    traverse_to_get_midsplited_point_2(root, midsplited_point);

		total_riemannian_length = 0.0;
	    for (i=0; i<(int)midsplited_point.size(); ++i)
		    total_riemannian_length += midsplited_point[i].riemannian_length;

        if (total_riemannian_length<length_threshold)
		{
            if (_mesh.is_boundary(vh0)&&_mesh.is_collapse_ok(heh0))
				_mesh.collapse(heh0);
			else if (_mesh.is_boundary(vh1)&&_mesh.is_collapse_ok(heh1))
				_mesh.collapse(heh1);
			else if (_mesh.is_collapse_ok(heh0))
				_mesh.collapse(heh0);
			else if (_mesh.is_collapse_ok(heh1))
				_mesh.collapse(heh1);
				
		}
	}

	return;
}

void SurfaceMesher::discretize_2d_mesh_inner_edges(SurfaceRegion *_sr, MyMesh &_mesh, double _tol)
{
	int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	Vec3d ru, rv;
	lts->calculate_first_derivative(Vec2d(0.0, 0.0), ru, rv); // can get 2d metric tensor from ru, rv
	Mat2by2d metric_tensor(ru | ru, ru | rv, ru | rv, rv | rv);

	HalfedgeHandle heh0, heh1;
	VertexHandle   vh0, vh1;
	EdgeHandle     eh; // for debug
	Vec3d pt(0.0, 0.0, 0.0);
	Vec2d uv(0.0, 0.0);
	double total_riemannian_length = 0.0, ideal_riemannian_length = 0.0;
	int i, ideal_edge_number;
	MeshPoint3D mp3d;

	const int NP_THRE = 10;
	// int id_mid_pt;
	static int times = 0;
	double edge_riemannian_length = 0.0;
	static SurfaceRegion *sr_old = NULL;
	if (0==times)
		sr_old = _sr;
	if (sr_old!=_sr)
	{
		sr_old = _sr;
		times = 0;
	}

	std::vector<MidsplitedPoint> midsplited_point;
	MidsplitedPoint msp;
	std::list<MeshPoint3D> *lst_mp3d = NULL;

	MyMesh::EdgeIter e_it(_mesh.edges_begin()),
		             e_end(_mesh.edges_end());

	// compute the entire 2d riemmanian length and the discretized points on inner mesh edges
    for (; e_it!=e_end; ++e_it)
	{
		eh = e_it.handle();
        heh0   = _mesh.halfedge_handle(e_it.handle(), 0);
	    heh1   = _mesh.halfedge_handle(e_it.handle(), 1);
        vh0    = _mesh.to_vertex_handle(heh0);
	    vh1    = _mesh.to_vertex_handle(heh1);

		edge_riemannian_length = _mesh.property(riemannian_length, e_it);

		if (_mesh.is_boundary(e_it))
		{
            _mesh.property(discretized_point, e_it) = new std::list<MeshPoint3D>;
		    lst_mp3d = _mesh.property(discretized_point, e_it);
            mp3d.point = _mesh.point(vh0);
			mp3d.isotropic_size = _mesh.property(iso3d_size, vh0);
			lst_mp3d->push_back(mp3d);
			mp3d.point = _mesh.point(vh1);
			mp3d.isotropic_size = _mesh.property(iso3d_size, vh1);
			lst_mp3d->push_back(mp3d);

			continue;
		}

		/**
		if (_mesh.property(riemannian_length, e_it)>0.0 && _mesh.property(riemannian_length, e_it)<1.5)
		{
            _mesh.property(discretized_point, e_it) = new std::list<MeshPoint3D>;
		    lst_mp3d = _mesh.property(discretized_point, e_it);
            mp3d.point = _mesh.point(vh0);
			mp3d.isotropic_size = _mesh.property(iso3d_size, vh0);
			lst_mp3d->push_back(mp3d);
			mp3d.point = _mesh.point(vh1);
			mp3d.isotropic_size = _mesh.property(iso3d_size, vh1);
			lst_mp3d->push_back(mp3d);

			continue;
		}
		*/

        // create binary tree to generate intermediate auxiliary points
        MidsplitedEdge *root = new MidsplitedEdge; // where to free memory?!
	    root->point[0] = _mesh.point(vh0);
	    root->point[1] = _mesh.point(vh1);
	    root->left = root->right = NULL;
	    root->parent = NULL;
	    riemannian_length_2d_line(_sr, root->point[0], root->point[1], root->riemannian_length);
        create_binary_tree_2(_sr, root); // By 20071210-2209, it had been "create_binary_tree(root);"

	    // store the intermediate auxiliary points as a vector
	    // std::vector<MidsplitedPoint> midsplited_point;
	    // MidsplitedPoint msp;
	    msp.point = root->point[0];
	    msp.riemannian_length = 0.0;
	    midsplited_point.push_back(msp);
	    // may need to output the results to test the correctness
	    traverse_to_get_midsplited_point_2(root, midsplited_point);


		// determine the total riemannian length, ideal edge number, ideal riemannian length
        ideal_riemannian_length = 0.0;
		total_riemannian_length = 0.0;
	    ideal_edge_number = 1;
	    for (i=0; i<(int)midsplited_point.size(); ++i)
	    {
		    total_riemannian_length += midsplited_point[i].riemannian_length;
			uv[0] = midsplited_point[i].point[0];
			uv[1] = midsplited_point[i].point[1];
			lts->evaluate(uv, pt);
		    midsplited_point[i].isotropic_size = mesh_size_spec_->get_final_spacing(pt);
		    if (i > 0)
			    midsplited_point[i].riemannian_length += midsplited_point[i-1].riemannian_length;
	    }
        ideal_edge_number = (int)total_riemannian_length;
	    if (total_riemannian_length - ideal_edge_number >= 0.5)
	     	ideal_edge_number++;
	    if (0==ideal_edge_number)
		    ideal_edge_number = 1;
	    ideal_riemannian_length = total_riemannian_length / (1.0 * ideal_edge_number);

		_mesh.property(riemannian_length, e_it) = total_riemannian_length; // but this edge will be discretized later, 

	    cout << "    Number of final discretized points: " << ideal_edge_number+1 << endl;

	    // generate final discretized points
		// actually it's MeshPoint2D! because of OpenMesh, we have to use 3D instead of 2D!
	    // std::list<MeshPoint3D> *lst_mp3d = NULL; // = new std::list<MeshPoint3D>; do not have "(ideal_edge_number+1)"!
		_mesh.property(discretized_point, e_it) = new std::list<MeshPoint3D>;
		lst_mp3d = _mesh.property(discretized_point, e_it);
		mp3d.point = midsplited_point[0].point;
		mp3d.isotropic_size = midsplited_point[0].isotropic_size;
		lst_mp3d->push_back(mp3d);

		if (ideal_edge_number>NP_THRE && times<2) 
		{ 
			// only take the middle point (had better not be more, or the error of edge length would be larger)
			// it's necessary, some problems could happen
			mp3d.point = (midsplited_point[0].point + midsplited_point[midsplited_point.size()-1].point) / 2.0;
			uv[0] = mp3d.point[0];
			uv[1] = mp3d.point[1];
			lts->evaluate(uv, pt);
			mp3d.isotropic_size = mesh_size_spec_->get_final_spacing(pt);
			lst_mp3d->push_back(mp3d);
		}
		else
		{
			int j, k = 0, m = 0;
			double interval, c, f0, fc, temp_ideal_riemannian_length = 0.0, temp_mid_riemannian_length = 0.0;
			const double threshold_length = 1.e-6; // _tol!
			for (i=1; i<=ideal_edge_number-1; ++i)
			{
				for (j=k; j<(int)midsplited_point.size()-1; ++j)
				{
					if (midsplited_point[j].riemannian_length < i * ideal_riemannian_length 
						&& i * ideal_riemannian_length <= midsplited_point[j+1].riemannian_length)
					{
						k = j;
						break;
					}
				}
				assert(midsplited_point.size()-1 != j); // find the index

				// get the final discretized point by recursion
				// attension and todo_clg: the efficiency should be quite low because of the 
				// recursion, it need not to be such precise, it should be treated more roughly
				if (fabs(i*ideal_riemannian_length-midsplited_point[j].riemannian_length) < threshold_length)
				{
					mp3d.point = midsplited_point[j].point;
					mp3d.isotropic_size = midsplited_point[j].isotropic_size;
					lst_mp3d->push_back(mp3d);
					continue;
				}
				else if (fabs(midsplited_point[j+1].riemannian_length-i*ideal_riemannian_length) < threshold_length)
				{
					mp3d.point = midsplited_point[j].point;
					mp3d.isotropic_size = midsplited_point[j].isotropic_size;
					lst_mp3d->push_back(mp3d);
					continue;
				}

				temp_ideal_riemannian_length = i * ideal_riemannian_length - midsplited_point[j].riemannian_length;
				c = interval = 0.5;
				f0 = midsplited_point[j+1].riemannian_length - midsplited_point[j].riemannian_length 
					- 0.5*(midsplited_point[j].isotropic_size - midsplited_point[j+1].isotropic_size);
				fc = f0*c+0.5*c*c*(midsplited_point[j].isotropic_size - midsplited_point[j+1].isotropic_size);
				m = 1;
				while(fabs(temp_ideal_riemannian_length-fc) > threshold_length && m <= 20)
				{
					interval *= 0.5;
					if (temp_ideal_riemannian_length < fc)
						c -= interval; 
					else
						c += interval;
					fc = f0*c+0.5*c*c*(midsplited_point[j].isotropic_size - midsplited_point[j+1].isotropic_size);
					m++;
				}

				if (fabs(temp_ideal_riemannian_length-fc) > threshold_length && m > 20)
				{
					cout << "discretize_2d_mesh_inner_edges(): iteration exceeded!" << endl;
					assert(false);
					return; // exit(0)? return false?
				}

				mp3d.point = midsplited_point[j].point + c*(midsplited_point[j+1].point - midsplited_point[j].point);
				uv[0] = mp3d.point[0];
				uv[1] = mp3d.point[1];
				lts->evaluate(uv, pt);
				mp3d.isotropic_size = mesh_size_spec_->get_final_spacing(pt);
				lst_mp3d->push_back(mp3d);

			}
		}

		mp3d.point = midsplited_point[midsplited_point.size()-1].point;
        mp3d.isotropic_size = midsplited_point[midsplited_point.size()-1].isotropic_size;
	    lst_mp3d->push_back(mp3d);

 	    midsplited_point.clear();
	}

	++times;

	return;
}

void SurfaceMesher::filter_param_mesh_inner_edge_discretized_point(SurfaceRegion *_sr, 
		MyMesh &_mesh, double _tol)
{
    int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);

	VertexHandle   vh; // for debug
	Vec3d pt0, pt1, pt;
	int sz0, sz1;
	double rl;
	const double min_riemannian_dist = 0.70710678;

    MyMesh::ConstVertexIter  v_it(_mesh.vertices_begin()),
			  		         v_end(_mesh.vertices_end());
	std::vector<MyMesh::EdgeHandle> vec_unit_edge; // no inserting points, include boundary edges
	std::vector<MyMesh::EdgeHandle> vec_non_unit_edge;
	std::vector<MyMesh::EdgeHandle>::iterator vec_unit_edge_iter, vec_non_unit_edge_iter, vec_non_unit_edge_iter1;
	MyMesh::EdgeHandle eh;
    MyMesh::ConstVertexOHalfedgeIter cvoh_it, temp_cvoh_it;
	std::list<MeshPoint3D> *lst_mp3d0 = NULL, *lst_mp3d1 = NULL;
	std::list<MeshPoint3D>::iterator lst_iter0, lst_iter1, lst_iter00, lst_iter11;

	for (; v_it!=v_end; ++v_it)
	{
		vh = v_it.handle();
		pt = _mesh.point(v_it.handle());

        for (cvoh_it=_mesh.cvoh_iter(v_it.handle()); cvoh_it; ++cvoh_it)
		{
			eh = _mesh.edge_handle(cvoh_it.handle());
            lst_mp3d0 = _mesh.property(discretized_point, eh);

			if (is_equal(pt, lst_mp3d0->back().point))
				lst_mp3d0->reverse();
			else if (!is_equal(pt, lst_mp3d0->front().point)) // for debug
				assert(false);

			if (2==lst_mp3d0->size())
				vec_unit_edge.push_back(eh);
			else
				vec_non_unit_edge.push_back(eh);
		}

        // first compare with inserted end-points of edges in vec_unit_edge, be aware of efficiency
        for (vec_non_unit_edge_iter=vec_non_unit_edge.begin(); vec_non_unit_edge_iter!=vec_non_unit_edge.end(); 
			++vec_non_unit_edge_iter)
		{
			for (vec_unit_edge_iter=vec_unit_edge.begin(); vec_unit_edge_iter!=vec_unit_edge.end();
				++vec_unit_edge_iter)
			{
				lst_mp3d0 = _mesh.property(discretized_point, *vec_unit_edge_iter);
				pt0 = lst_mp3d0->back().point; // here only compare with the other end point

                lst_mp3d1 = _mesh.property(discretized_point, *vec_non_unit_edge_iter);
                lst_iter0 = lst_mp3d1->begin();
			    lst_iter1 = lst_mp3d1->end();
			    ++lst_iter0;
			    --lst_iter1;
			    for (; lst_iter0!=lst_iter1; )
			    {
					pt1 = lst_iter0->point;
                    riemannian_length_2d_line(_sr, pt0, pt1, rl);
					if (rl<min_riemannian_dist)
				  		lst_iter0 = lst_mp3d1->erase(lst_iter0);
				    else
						break;
			    }
			}
		}

		// then compare among edges with more than 2 discretized points
        for (vec_non_unit_edge_iter=vec_non_unit_edge.begin(); vec_non_unit_edge_iter!=vec_non_unit_edge.end(); 
			++vec_non_unit_edge_iter)
		{
			lst_mp3d0 = _mesh.property(discretized_point, *vec_non_unit_edge_iter);
			lst_iter1 = lst_mp3d0->end();
			--lst_iter1;

			vec_non_unit_edge_iter1 = vec_non_unit_edge_iter + 1;
			for (; vec_non_unit_edge_iter1!=vec_non_unit_edge.end(); ++vec_non_unit_edge_iter1)
			{
				lst_iter0 = lst_mp3d0->begin();
				++lst_iter0;

				lst_mp3d1 = _mesh.property(discretized_point, *vec_non_unit_edge_iter1);
                
				for (; lst_iter0!=lst_iter1; ++lst_iter0)
				{
					pt0 = (*lst_iter0).point;

					lst_iter00 = lst_mp3d1->begin();
				    ++lst_iter00;
				    lst_iter11 = lst_mp3d1->end();
				    --lst_iter11;
					for (; lst_iter00!=lst_iter11; )
					{
                        pt1 = (*lst_iter00).point;
                        riemannian_length_2d_line(_sr, pt0, pt1, rl);
                        if (rl<min_riemannian_dist)
				  		    lst_iter00 = lst_mp3d1->erase(lst_iter00);
						else
							break;
					}
				}
			}
		}

		vec_unit_edge.clear();
		vec_non_unit_edge.clear();
	}


	MyMesh::ConstEdgeIter e_it(_mesh.edges_begin()),
		                  e_end(_mesh.edges_end());
	for (; e_it!=e_end; ++e_it)
	{
		if (_mesh.is_boundary(e_it))
			continue;

        lst_mp3d0 = _mesh.property(discretized_point, e_it);
		cout << "No. of discretized points after filtering: " << lst_mp3d0->size() << endl;
	}

	return;
}

// should return the area value "double type", for only return one value!
// this function is not convenient to use!
void SurfaceMesher::area2(Vec3d _pt0, Vec3d _pt1, Vec3d _pt2, double &_area)
{
    _area = (_pt1[0]-_pt0[0])*(_pt2[1]-_pt0[1]) - (_pt2[0]-_pt0[0])*(_pt1[1]-_pt0[1]);

	return;
}

// the result angle is in [-0.5*PI, 0.5*PI], not in [0.0, PI].
// this interface is not reasonable and versatile, be careful
double SurfaceMesher::angle(Vec3d _pt0, Vec3d _pt1, Vec3d _pt2)
{
	Vec3d v10 = _pt0 - _pt1, v12 = _pt2 - _pt1;
	double denom = v10.norm() * v12.norm();
	if (is_zero(denom))
		return 0;

	double ang, cos_a = (v10 | v12) / (v10.norm() * v12.norm());
	if (cos_a<-1.0)
		cos_a = -1.0;
	if (cos_a>1.0)
		cos_a = 1.0;
	ang = acos(cos_a);

	return ang;
}

bool SurfaceMesher::is_flip_legal(MyMesh &_mesh, MyMesh::EdgeHandle &_eh, MyMesh::VertexHandle &_vh)
{
    // assert(!_mesh.is_boundary(_eh));

	if (_mesh.is_boundary(_eh)) // could be bounary edge
	{
		cout << "Flipping is not legal because it's a boundary edge" << endl;
		return false;
	}
	if (_mesh.status(_eh).deleted())
	{
		cout << "Flipping is not legal because it's been deleted" << endl;
		return false;
	}
	if (!_mesh.is_flip_ok(_eh))
	{
        cout << "Flipping is not legal because is_flip_ok() returns false" << endl;
		return false;
	}

    Vec3d pt, pt0, pt1, pt2;
	MyMesh::HalfedgeHandle heh0, heh1;
	MyMesh::VertexHandle vh, vh0, vh1;
	double ar0, ar1, ang;
	const double ang_thre = 0.005;

    heh0 = _mesh.halfedge_handle(_eh, 0);
	heh1 = _mesh.halfedge_handle(_eh, 1);

	vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh0));
    vh1 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh1));

	if (_vh!=vh0)
	{
        assert(_vh==vh1);
        vh = vh0;
	}
	else
		vh = vh1;

	pt = _mesh.point(_vh);
	pt0 = _mesh.point(_mesh.to_vertex_handle(heh0));
	pt1 = _mesh.point(_mesh.to_vertex_handle(heh1));
	pt2 = _mesh.point(vh);

	// is visible?
    area2(pt, pt0, pt2, ar0);
	area2(pt, pt1, pt2, ar1);
	if ((ar0>-1.0e-7&&ar1>-1.0e-7) || (ar0<1.0e-7&&ar1<1.0e-7))
	{
		cout << "Flipping is not legal because it's not visible" << endl;
		return false;
	}

    // are new elements abnormal?
    ang = angle(pt0, pt, pt2);
	if (fabs(ang)<ang_thre)
	{
		cout << "Flipping is not legal because angle is too small" << endl;
		return false;
	}
	ang = angle(pt2, pt, pt1);
	if (fabs(ang)<ang_thre)
	{
		cout << "Flipping is not legal because angle is too small" << endl;
		return false;
	}
	ang = angle(pt, pt2, pt0);
	if (fabs(ang)<ang_thre)
	{
		cout << "Flipping is not legal because angle is too small" << endl;
		return false;
	}
    ang = angle(pt1, pt2, pt);
	if (fabs(ang)<ang_thre)
	{
		cout << "Flipping is not legal because angle is too small" << endl;
        return false;
	}

    // passed all the test
	cout << "Flipping is legal" << endl;
	return true;

    /**
	if ((ar0>1.0e-5 && ar1<-1.0e-5) || (ar0<-1.0e-5 && ar1>1.0e-5))
	{
		cout << "Flipping is legal" << endl;
		return true;
	}
	else
	{
		cout << "Flipping is not legal because it's not visible" << endl;
	    return false;
	}
	*/
}

void SurfaceMesher::calc_riemannian_circumcenter(SurfaceRegion *_sr, MyMesh &_mesh, MyMesh::FaceHandle &_fh, Vec3d &_cc)
{
	int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);
    Vec3d v0, v1;
	lts->calculate_first_derivative(Vec2d(0.0, 0.0), v0, v1); // can get 2d metric tensor from ru, rv
	Mat2by2d metric_tensor(v0 | v0, v0 | v1, v0 | v1, v1 | v1);

	Mat2by2d mt[3];
	Vec3d pt[3];
	double spa;
	int i;
	MyMesh::ConstFaceVertexIter cfv_it;
	for (cfv_it=_mesh.cfv_iter(_fh), i=0; cfv_it; ++cfv_it, ++i)
	{
		pt[i] = _mesh.point(cfv_it);
        spa = _mesh.property(iso3d_size, cfv_it.handle());
		mt[i] = metric_tensor * (1.0/(spa*spa));
	}

	// initialization
    _cc = (pt[0] + pt[1] + pt[2]) / 3.0; 
	Mat2by1d fx, du;
	Mat2by2d dfx;
	Mat1by1d val;
	bool succ = false;
	const double tol = 1.0e-6;
	const int max_iter_num = 100;
	for (i=0; i<max_iter_num; ++i)
	{
		// val = ~du * mt[0] * du;
	    // fx(0, 0) = val(0, 0); // should provide a type transformation from Mat1by1d to double

		// use newton method to solve the nonlinear equation set
        du = Mat2by1d(pt[0][0]-_cc[0], pt[0][1]-_cc[1]);
		val = ~du * mt[0] * du;
        fx(0, 0) = fx(1, 0) = val(0,0);
        du = Mat2by1d(pt[1][0]-_cc[0], pt[1][1]-_cc[1]);
		val = ~du * mt[1] * du;
        fx(0, 0) -= val(0,0);
        du = Mat2by1d(pt[2][0]-_cc[0], pt[2][1]-_cc[1]);
		val = ~du * mt[2] * du;
        fx(1, 0) -= val(0,0);

		dfx(0, 0) = 2*(mt[0](0,0)-mt[1](0,0))*_cc[0] + 2*(mt[0](0,1)-mt[1](0,1))*_cc[1]
		    - 2*(mt[0](0,0)*pt[0][0]-mt[1](0,0)*pt[1][0]) - 2*(mt[0](0,1)*pt[0][1]-mt[1](0,1)*pt[1][1]);
		dfx(0, 1) = 2*(mt[0](1,1)-mt[1](1,1))*_cc[1] + 2*(mt[0](0,1)-mt[1](0,1))*_cc[0]
		    - 2*(mt[0](1,1)*pt[0][1]-mt[1](1,1)*pt[1][1]) - 2*(mt[0](0,1)*pt[0][0]-mt[1](0,1)*pt[1][0]);
		dfx(1, 0) = 2*(mt[0](0,0)-mt[2](0,0))*_cc[0] + 2*(mt[0](0,1)-mt[2](0,1))*_cc[1]
		    - 2*(mt[0](0,0)*pt[0][0]-mt[2](0,0)*pt[2][0]) - 2*(mt[0](0,1)*pt[0][1]-mt[2](0,1)*pt[2][1]);
		dfx(1, 1) = 2*(mt[0](1,1)-mt[2](1,1))*_cc[1] + 2*(mt[0](0,1)-mt[2](0,1))*_cc[0]
		    - 2*(mt[0](1,1)*pt[0][1]-mt[2](1,1)*pt[2][1]) - 2*(mt[0](0,1)*pt[0][0]-mt[2](0,1)*pt[2][0]);

		du = dfx.inversion() * fx;
        _cc[0] -= du(0, 0);
		_cc[1] -= du(1, 0);

		if (sqrt(du(0,0)*du(0,0)+du(1,0)*du(1,0))<tol)
		{
			succ = true;
			break;
		}
	}
	if (succ) // check if in parametric space
	{
		// can be not in the parametric space, but it seems ok
        // assert(_cc[0]>=0.0 && _cc[0]<=1.0 && _cc[1]>=0.0 && _cc[1]<=1.0);
        // assert(_cc[0]+_cc[1]-1.0<=0.0);
	}
	else
	{
	    // assert(false); // todo_clg: should be valid
	    cout << "calc_riemannian_circumcenter: maximum number of iterations exceeded" << endl;
	}

	return;
}

bool SurfaceMesher::is_delaunay_broken(SurfaceRegion *_sr, MyMesh &_mesh, MyMesh::VertexHandle &_vh, MyMesh::FaceHandle &_fh)
{
    Vec3d cc;
    calc_riemannian_circumcenter(_sr, _mesh, _fh, cc);

    int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);
    Mat2by2d metric_tensor = lts->get_metric_tensor();
    Vec3d pt = _mesh.point(_vh), pt0 = _mesh.point(_mesh.cfv_iter(_fh).handle());
	double sz;
	Mat2by1d op = Mat2by1d(pt[0]-cc[0], pt[1]-cc[1]), op0 = Mat2by1d(pt0[0]-cc[0], pt0[1]-cc[1]);
	Mat1by1d val;
	double rl0, rl1, dm = 0.0; // delaunay measure
	
	// first deal with the inserting point
	sz = _mesh.property(iso3d_size, _vh);
	val = ~op * (metric_tensor/(sz*sz)) * op;
    rl0 = sqrt(val(0,0));
	val = ~op0 * (metric_tensor/(sz*sz)) * op0;
	rl1 = sqrt(val(0,0));
	dm += rl0/rl1;

	// then three face end-points
	MyMesh::ConstFaceVertexIter cfv_it;
	for (cfv_it=_mesh.cfv_iter(_fh); cfv_it; ++cfv_it)
	{
        sz = _mesh.property(iso3d_size, cfv_it.handle());
		val = ~op * (metric_tensor/(sz*sz)) * op;
		rl0 = sqrt(val(0, 0));
		val = ~op0 * (metric_tensor/(sz*sz)) * op0;
	    rl1 = sqrt(val(0, 0));
	    dm += rl0/rl1;
	}

	if (dm<4)
		return true;
	else
	    return false;
}

void SurfaceMesher::insert_one_point(SurfaceRegion *_sr, MeshPoint3D &_mp3d, MyMesh &_mesh)
{
	MyMesh::FaceIter f_it(_mesh.faces_begin()),
		             f_end(_mesh.faces_end());
	std::vector<FaceHandle> vec_base_fh;
	FaceHandle fh;
	Vec3d pt = _mp3d.point, pt0, pt1;
	MyMesh::FaceHalfedgeIter fh_it;
	double area;
	int cnt0;
	MyMesh::FaceFaceIter ff_it;
	bool bval;
	EdgeHandle eh, eh0, eh1;
	std::deque<EdgeHandle> deq_eh;
	HalfedgeHandle heh;
	VertexHandle vh, vh0;
	int fid, eid;

    // first find one of base elements
    for (; f_it!=f_end; ++f_it)
	{
		cnt0 = 0;
		fid = f_it.handle().idx(); // for debug

        for (fh_it=_mesh.fh_iter(f_it.handle()); fh_it; ++fh_it)
		{
			eid = _mesh.edge_handle(fh_it).idx(); // for debug

            pt0 = _mesh.point(_mesh.from_vertex_handle(fh_it.handle()));
			pt1 = _mesh.point(_mesh.to_vertex_handle(fh_it.handle()));
            area2(pt, pt0, pt1, area);

			if (area<-1.0e-8) // todo_clg: originally -1.0e-12, but some examples failed!
				break;

			cnt0++;
		}

		if (cnt0==3)
			vec_base_fh.push_back(f_it.handle());
	}
	cout << "No. of found enclosing base elements: " << vec_base_fh.size() << endl;
	// assert(vec_base_fh.size()==1 || vec_base_fh.size()==2); // should be valid
	if (1!=vec_base_fh.size() && 2!=vec_base_fh.size())
	 	return;

	bval = false;
	if (2==vec_base_fh.size()) // two base faces must be adjacent
	{
        for (ff_it=_mesh.ff_iter(vec_base_fh[0]); ff_it; ++ff_it)
			if (ff_it.handle()==vec_base_fh[1])
			{
                bval = true;
				break;
			}

		// assert(bval); // todo_clg: why here bval could be false?!
		if (!bval)
			return;
	}
	
	// start from the base element
    if (2==vec_base_fh.size())
	{
        for (fh_it=_mesh.fh_iter(vec_base_fh[0]); fh_it; ++fh_it)
		{
            fh = _mesh.face_handle(_mesh.opposite_halfedge_handle(fh_it.handle()));
			if (fh==vec_base_fh[1])
				break;
		}
		assert(fh_it);
		eh = _mesh.edge_handle(fh_it.handle()); // shared edge

        heh = _mesh.halfedge_handle(eh, 0);
		deq_eh.push_back(_mesh.edge_handle(_mesh.next_halfedge_handle(heh)));
		deq_eh.push_back(_mesh.edge_handle(_mesh.prev_halfedge_handle(heh)));
		heh = _mesh.halfedge_handle(eh, 1);
        deq_eh.push_back(_mesh.edge_handle(_mesh.next_halfedge_handle(heh)));
		deq_eh.push_back(_mesh.edge_handle(_mesh.prev_halfedge_handle(heh)));

		// are the previous edge handles above still valid after these kind of
		// topological transformations (and add new vertex) without applying 
		// garbage_collection() ?
		vh = _mesh.add_vertex(pt);
		_mesh.split(eh, vh);
		_mesh.property(iso3d_size, vh) = _mp3d.isotropic_size;
        
        while (!deq_eh.empty())
		{
			eh = deq_eh.front(); // and vh
			deq_eh.pop_front();
			heh = _mesh.halfedge_handle(eh, 0);
            vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));
			if (vh0!=vh)
			{
				fh = _mesh.face_handle(heh);
				eh0 = _mesh.edge_handle(_mesh.prev_halfedge_handle(heh));
				eh1 = _mesh.edge_handle(_mesh.next_halfedge_handle(heh));

				// heh = _mesh.halfedge_handle(eh, 1);
                // vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));
				// assert(vh0==vh);
			}
			else
			{
                heh = _mesh.halfedge_handle(eh, 1);
				eh0 = _mesh.edge_handle(_mesh.prev_halfedge_handle(heh));
				eh1 = _mesh.edge_handle(_mesh.next_halfedge_handle(heh));
                // vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));

				fh = _mesh.face_handle(heh);
			}

			// two things: is flip legal? is delaunay criteria broken?
            if (is_flip_legal(_mesh, eh, vh) && is_delaunay_broken(_sr, _mesh, vh, fh))
			{
                _mesh.flip(eh);

				// add new edges, be careful: heh has been changed!
				deq_eh.push_back(eh0);
				deq_eh.push_back(eh1);
			}
		}
	}
	else // 1==vec_base_fh.size()
	{
        for (fh_it=_mesh.fh_iter(vec_base_fh[0]); fh_it; ++fh_it)
            deq_eh.push_back(_mesh.edge_handle(fh_it.handle()));

        vh = _mesh.add_vertex(pt);
		_mesh.split(vec_base_fh[0], vh);
		_mesh.property(iso3d_size, vh) = _mp3d.isotropic_size;

		while (!deq_eh.empty())
		{
			eh = deq_eh.front(); // and vh
			deq_eh.pop_front();
			heh = _mesh.halfedge_handle(eh, 0);
            vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));
			if (vh0!=vh)
			{
				fh = _mesh.face_handle(heh);
				eh0 = _mesh.edge_handle(_mesh.prev_halfedge_handle(heh)); // eh0 and eh1 wouldn't be changed after edge flip
				eh1 = _mesh.edge_handle(_mesh.next_halfedge_handle(heh));

				// heh = _mesh.halfedge_handle(eh, 1);
                // vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));
				// assert(vh0==vh);
			}
			else
			{
                heh = _mesh.halfedge_handle(eh, 1);
				eh0 = _mesh.edge_handle(_mesh.prev_halfedge_handle(heh));
				eh1 = _mesh.edge_handle(_mesh.next_halfedge_handle(heh));
                // vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));

				fh = _mesh.face_handle(heh);
			}

			// two things: is flip legal? is delaunay criteria broken?
            if (is_flip_legal(_mesh, eh, vh) && is_delaunay_broken(_sr, _mesh, vh, fh))
			{
                _mesh.flip(eh);

				// add new edges
				deq_eh.push_back(eh0);
				deq_eh.push_back(eh1);
			}
		}
	}

	_mesh.garbage_collection();

	// assert(OpenMesh::IO::write_mesh(_mesh, "result_param_mesh.off"));
	// assert(OpenMesh::IO::write_mesh(mesh1, "result_physical_mesh.off"));

	return;
}

void SurfaceMesher::insert_filtered_discretized_point(SurfaceRegion *_sr, MyMesh &_mesh, double _tol)
{
	/** commented on 20080108-0910, it's valid, but the order of inserting points is not 
	* the cause of extremely abnormal mesh elements, so i'd rather use initial order of 
	* inserting points.

	int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);

	std::list<MeshPoint3D> inserting_point;
	std::list<MeshPoint3D> *lst_mp3d = NULL;
	std::list<MeshPoint3D>::iterator lst_iter;
	MeshPoint3D mp3d;
	MyMesh::ConstEdgeIter e_it, e_end(_mesh.edges_end());
	bool go_on = true;

	while (go_on)
	{
		go_on = false;

		e_it = _mesh.edges_begin();
		for (; e_it!=e_end; ++e_it)
		{
			if (_mesh.is_boundary(e_it))
				continue;

			lst_mp3d = _mesh.property(discretized_point, e_it);
			if (lst_mp3d->size()==2)
				continue;
	        
			lst_iter = lst_mp3d->begin();
			lst_iter++;
			inserting_point.push_back(*lst_iter);
			lst_mp3d->erase(lst_iter);
			if (!go_on)
			    go_on = true;
		}
	}

	for (e_it=_mesh.edges_begin(); e_it!=_mesh.edges_end(); ++e_it)
	{
        lst_mp3d = _mesh.property(discretized_point, e_it);
		lst_mp3d->clear();
		delete lst_mp3d;
		lst_mp3d = NULL;
	}

	cout << "No. of inserting points: " << inserting_point.size() << endl;

	// insert new generated/discretized points
    for (lst_iter=inserting_point.begin(); lst_iter!=inserting_point.end(); ++lst_iter)
        insert_one_point(_sr, *lst_iter, _mesh);

	return;
	*/

    int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);

	std::list<MeshPoint3D> inserting_point;
	std::list<MeshPoint3D> *lst_mp3d = NULL;
	std::list<MeshPoint3D>::iterator lst_iter, lst_iter0;
	MeshPoint3D mp3d;
	MyMesh::ConstEdgeIter e_it(_mesh.edges_begin()),
		                  e_end(_mesh.edges_end());

	for (; e_it!=e_end; ++e_it)
	{
		if (_mesh.is_boundary(e_it))
			continue;

        lst_mp3d = _mesh.property(discretized_point, e_it);
		if (lst_mp3d->size()==2)
			continue;
        
		lst_iter0 = lst_mp3d->end();
		lst_iter0--;
		lst_iter = lst_mp3d->begin();
		lst_iter++;
		for (; lst_iter!=lst_iter0; ++lst_iter)
			inserting_point.push_back(*lst_iter);
	}
	cout << "No. of inserting points: " << inserting_point.size() << endl;

	// insert new generated/discretized points
    for (lst_iter=inserting_point.begin(); lst_iter!=inserting_point.end(); ++lst_iter)
        insert_one_point(_sr, *lst_iter, _mesh);

	return;
}

void SurfaceMesher::calc_element_quality(SurfaceRegion *_sr, MyMesh &_mesh, MyMesh::VertexHandle _vh[3], double &_qua)
{
    int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);

	int i, j;
	Mat2by1d edge[3];
	Vec3d vec[3];
	vec[0] = _mesh.point(_vh[1]) - _mesh.point(_vh[0]), 
	vec[1] = _mesh.point(_vh[2]) - _mesh.point(_vh[0]),
	vec[2] = _mesh.point(_vh[2]) - _mesh.point(_vh[1]);
	for (i=0; i<3; ++i)
		for (j=0; j<2; ++j)
			edge[i](j, 0) = vec[i][j];

	double det0 = Mat2by2d(vec[0][0], vec[1][0], vec[0][1], vec[1][1]).determinant();
	Mat2by2d mtp, mtg = lts->get_metric_tensor();
	_qua = 1.0e+30;
	double sz, numerator, denominator;
	const double sqrt3 = 1.73205080757;
	Mat1by1d rl;
	for (i=0; i<3; ++i)
	{
        sz = _mesh.property(iso3d_size, _vh[i]);
        mtp = mtg / (sz*sz);
        numerator = 2.0*sqrt3*fabs(sqrt(mtp.determinant())*det0); // 20080101: todo_clg: originally there's no sqrt!

		denominator = 0.0;
        for (j=0; j<3; ++j)
		{
			rl = ~edge[j] * mtp * edge[j];
			denominator += rl(0, 0);
		}

		if (numerator/denominator<_qua)
			_qua = numerator / denominator;
	}

	return;
}

void SurfaceMesher::optimize_with_diagonal_swapping_given_threshold(SurfaceRegion *_sr, MyMesh &_mesh, double _qua_ratio_threshold)
{
    MyMesh::EdgeIter e_it;
    double qua_ratio, qua0, qua1, qua;
	MyMesh::VertexHandle vh[4], vh_param[3];
	MyMesh::HalfedgeHandle heh0, heh1;
	MyMesh::EdgeHandle eh;
	
	for (e_it=_mesh.edges_begin(); e_it!=_mesh.edges_end(); ++e_it)
	{
		if (_mesh.is_boundary(e_it.handle()))
			continue;
        
		eh     = e_it.handle();
        heh0   = _mesh.halfedge_handle(eh, 0);
		heh1   = _mesh.halfedge_handle(eh, 1);
		vh[0]  = _mesh.to_vertex_handle(heh0);
        vh[1]  = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh0));
		vh[2]  = _mesh.to_vertex_handle(heh1);
		vh[3]  = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh1));
        
        vh_param[0] = vh[0];
		vh_param[1] = vh[1];
        vh_param[2] = vh[2];
		calc_element_quality(_sr, _mesh, vh_param, qua0);   
        vh_param[0] = vh[0];
		vh_param[1] = vh[2];
        vh_param[2] = vh[3];
		calc_element_quality(_sr, _mesh, vh_param, qua);
		if (qua<qua0)
			qua0 = qua;

        vh_param[0] = vh[0];
		vh_param[1] = vh[1];
        vh_param[2] = vh[3];
		calc_element_quality(_sr, _mesh, vh_param, qua1);
		vh_param[0] = vh[1];
		vh_param[1] = vh[2];
        vh_param[2] = vh[3];
		calc_element_quality(_sr, _mesh, vh_param, qua);
		if (qua<qua1)
			qua1 = qua;

		if (qua1/qua0>_qua_ratio_threshold && is_flip_legal(_mesh, eh, vh[1]))
		{
			cout << "optimize with edge swapping ..." << endl;
		    _mesh.flip(eh);
		}
	}
    _mesh.garbage_collection();

	return;
}

void SurfaceMesher::optimize_with_diagonal_swapping(SurfaceRegion *_sr, MyMesh &_mesh)
{
	cout << "optimize mesh with diagonal swapping given threshold 2.0" << endl;
	optimize_with_diagonal_swapping_given_threshold(_sr, _mesh, 2.0);
	cout << "optimize mesh with diagonal swapping given threshold 1.5" << endl;
    optimize_with_diagonal_swapping_given_threshold(_sr, _mesh, 1.5);
	cout << "optimize mesh with diagonal swapping given threshold 1.0" << endl;
	optimize_with_diagonal_swapping_given_threshold(_sr, _mesh, 1.0);

	return;
}

void SurfaceMesher::map_mesh(SurfaceRegion *_sr, MyMesh &_mesh_src, MyMesh &_mesh_tar)
{
    int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	_mesh_tar.assign_connectivity(_mesh_src); // doesn't copy vertices coordinates

	Vec2d uv;
	Vec3d pt;
	MyMesh::VertexIter v_it_src(_mesh_src.vertices_begin()), v_it_tar(_mesh_tar.vertices_begin());
	
	for (; v_it_src!=_mesh_src.vertices_end() && v_it_tar!=_mesh_tar.vertices_end(); ++v_it_src, ++v_it_tar)
	{
        pt = _mesh_src.point(v_it_src.handle());
		uv[0] = pt[0];
		uv[1] = pt[1];
		lts->evaluate(uv, pt);
		_mesh_tar.point(v_it_tar.handle()) = pt;
	}

	return;
}

bool SurfaceMesher::seg_seg_int(Vec3d _p0, Vec3d _p1, Vec3d _p2, Vec3d _p3, double& _coef0, 
		double& _coef1, Vec3d& _pt_int)
{
	double num, denom;

	denom = _p0[0]*(_p3[1]-_p2[1]) + _p1[0]*(_p2[1]-_p3[1])
		  + _p2[0]*(_p0[1]-_p1[1]) + _p3[0]*(_p1[1]-_p0[1]);

	if (is_zero(denom, 1.0e-10)) // parallel
		return false;

	num = _p0[0]*(_p3[1]-_p2[1]) + _p2[0]*(_p0[1]-_p3[1]) + _p3[0]*(_p2[1]-_p0[1]);
	_coef0 = num / denom;
	num = -(_p0[0]*(_p2[1]-_p1[1]) + _p1[0]*(_p0[1]-_p2[1]) + _p2[0]*(_p1[1]-_p0[1]));
	_coef1 = num / denom;

	_pt_int = _p0 + _coef0*(_p1-_p0);

	return true;
}

void SurfaceMesher::smooth_with_unit_length(SurfaceRegion *_sr, MyMesh &_mesh)
{

	return;
}

void SurfaceMesher::smooth_with_optimal_shape(SurfaceRegion *_sr, MyMesh &_mesh)
{
	int surface_id = _sr->surface_id;
	LinearTriangleSurface *lts = (LinearTriangleSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	Mat2by2d mat, mtp, mtg = lts->get_metric_tensor();

	MyMesh::VertexIter v_it;
	MyMesh::VertexOHalfedgeIter voh_it;
	MyMesh::HalfedgeHandle heh;
	Mat2by1d ip[3], p0, p01;
	std::vector<Vec3d> vec_ip;
	Vec3d pt, pt_cen, pt_int, p[2];
	MyMesh::VertexHandle vh[3];
	double sz, d, ar;
	double coef, coef0, coef1;
	int i;

	for (v_it=_mesh.vertices_begin(); v_it!=_mesh.vertices_end(); ++v_it)
	{
		if (_mesh.is_boundary(v_it.handle()))
			continue;

		vh[2] = v_it.handle();
		pt_cen = _mesh.point(v_it.handle());

        for (voh_it=_mesh.voh_iter(v_it.handle()); voh_it; ++voh_it)
		{
            heh = _mesh.next_halfedge_handle(voh_it.handle());
			vh[0] = _mesh.from_vertex_handle(heh);
			vh[1] = _mesh.to_vertex_handle(heh);
            p[0] = _mesh.point(vh[0]);
			p[1] = _mesh.point(vh[1]);
			p0(0, 0) = p[0][0];
			p0(1, 0) = p[0][1];
			p01(0, 0) = p[1][0] - p[0][0];
			p01(1, 0) = p[1][1] - p[0][1];

            for (i=0; i<3; ++i)
			{
                sz = _mesh.property(iso3d_size, vh[i]);
                mtp = mtg / (sz*sz);
                d = sqrt((mtp(0, 0)*mtp(1, 1) - mtp(0, 1)*mtp(0, 1)) / 3.0);
				mat(0, 0) = d - mtp(0, 1);
				mat(0, 1) = -mtp(1, 1);
				mat(1, 0) = mtp(0, 0);
				mat(1, 1) = d + mtp(0, 1);
				mat /= (2.0*d);

				ip[i] = p0 + mat*p01;
			}

            pt[0] = (ip[0](0, 0) + ip[1](0, 0) + ip[2](0, 0)) / 3.0;
			pt[1] = (ip[0](1, 0) + ip[1](1, 0) + ip[2](1, 0)) / 3.0;
			pt[2] = 0.0;

			vec_ip.push_back(pt);
		}

		pt = Vec3d(0.0, 0.0, 0.0);
		for (i=0; i<vec_ip.size(); ++i)
            pt += vec_ip[i];
		pt /= vec_ip.size();
		vec_ip.clear();

		// avoid generating invalid elements
		// compute the intersection points with the coefficient of each bounding edge and 
		// the edge (pt_cen, pt)
		coef = 1.0;
        for (voh_it=_mesh.voh_iter(v_it.handle()); voh_it; ++voh_it)
		{
            heh = _mesh.next_halfedge_handle(voh_it.handle());
			vh[0] = _mesh.from_vertex_handle(heh);
			vh[1] = _mesh.to_vertex_handle(heh);
            p[0] = _mesh.point(vh[0]);
			p[1] = _mesh.point(vh[1]);
            // assert(seg_seg_int(pt_cen, pt, p[0], p[1], coef0, coef1, pt_int)); // todo_clg: commented on 20080106-1506
			if (!seg_seg_int(pt_cen, pt, p[0], p[1], coef0, coef1, pt_int)) // parallel
				continue;
			else if (coef0>=0.0 && coef0<coef)
				coef = coef0;
		}
		if (coef<1.0)
			coef *= 0.5;

		_mesh.point(v_it.handle()) = pt_cen + coef*(pt-pt_cen);
	}

	return;
}

void SurfaceMesher::generate_final_parametric_mesh_for_linear_triangle_surface_region(SurfaceRegion *_sr, 
		SurfaceMesh *_sm, double _tol)
{
    MyMesh mesh, *mesh1 = new MyMesh;

	// OpenMesh::EPropHandleT<double> riemmanian_length; // if put here, it couldn't referenced by other functions, 
	mesh.add_property(riemannian_length);
	// mesh1.add_property(riemmanian_length);

	// OpenMesh::EPropHandleT<std::deque<MeshPoint3D> *> discretized_point; // Vec2d is enough
    mesh.add_property(discretized_point);
	// mesh1.add_property(discretized_point);

	// OpenMesh::VPropHandleT<double> iso3d_size; // set as global variable
	mesh.add_property(iso3d_size);
	// mesh.add_property(iso3d_size);

	mesh1->add_property(global_id);

    assert(OpenMesh::IO::read_mesh(mesh, "boundary_param_mesh.off"));

	// calculate the spacing value on each initial boundary vertex
	Vec3d pt;
	double spa;
	MyMesh::ConstVertexIter v_it(mesh.vertices_begin()),
		                    v_end(mesh.vertices_end());
	for (; v_it!=v_end; ++v_it)
	{
        pt = mesh.point(v_it.handle());
        spa = mesh_size_spec_->get_final_spacing(pt);
		mesh.property(iso3d_size, v_it.handle()) = spa;
	}

	int npt_new = 0, npt_old = mesh.n_vertices();
	while (true)
	{
		// collapse_2d_short_edges(_sr, mesh, _tol);

	    // discretize inner edges of parametric mesh just like discretize line curve
        discretize_2d_mesh_inner_edges(_sr, mesh, _tol);

	    // filter the discretized points on inner mesh edges
        filter_param_mesh_inner_edge_discretized_point(_sr, mesh, _tol);

	    // insert filtered discretized points
        insert_filtered_discretized_point(_sr, mesh, _tol);

		optimize_with_diagonal_swapping(_sr, mesh);

		for (int i=0; i<3; ++i)
		    smooth_with_optimal_shape(_sr, mesh);

        // map_mesh(_sr, mesh, mesh1);

		// assert(OpenMesh::IO::write_mesh(mesh, "result_param_mesh.off"));
		// assert(OpenMesh::IO::write_mesh(mesh1, "result_physical_mesh.off"));

		npt_new = mesh.n_vertices();

		if (npt_new==npt_old)
			break;
		else
		    npt_old = npt_new;
	}

	map_mesh(_sr, mesh, *mesh1);
	all_surface_region_mesh_.push_back(mesh1);
	assert(OpenMesh::IO::write_mesh(mesh, "result_param_mesh.off"));
	assert(OpenMesh::IO::write_mesh(*mesh1, "result_physical_mesh.off"));

	/**
	// store the generated independent one-surface-region mesh
	_sm->surface_region_id = _sr->surface_id; // todo_clg (20080114-1024): two different items 
	for (v_it=mesh1.vertices_begin(); v_it!=mesh1.vertices_end(); ++v_it)
		_sm->mesh_point_3d.push_back(mesh1.point(v_it.handle()));

	// faces (indices starting at 0)
	MyMesh::ConstFaceIter cf_it(mesh1.faces_begin());
	MyMesh::ConstFaceVertexIter cfv_it;
	Vec3i elem;
    for (; cf_it!=mesh1.faces_end(); ++cf_it)
	{
        cfv_it = mesh1.cfv_iter(cf_it.handle());
		elem[0] = cfv_it.handle().idx();
		++cfv_it;
		elem[1] = cfv_it.handle().idx();
		++cfv_it;
		elem[2] = cfv_it.handle().idx();

		_sm->mesh_element.push_back(elem);
	}
	*/

	return;
}

// generalized delauney triangulation in parametric space
void SurfaceMesher::dt_for_linear_triangle_surface_region(SurfaceRegion *_sr, double _tol)
{	
	SurfaceMesh *sm = new SurfaceMesh;
	// sm->surface_region_id = 0; // will be assigned outside
	// assumption: only one loop here
	sm->boundary_point_number = (int)one_initial_front_param_point_[0]->size();

	// because of the round-off error, the discrete boundary points
	// don't exactly on the parametric triangle, which need a minor modification
	// for (i=0; i<(int)one_initial_front_param_point_[0]->size(); ++i)
		// sm->mesh_point_2d.push_back(one_initial_front_param_point_[0]->at(i));

    // start generalized delauney triangulation 
    // generate boundary mesh
    generate_boundary_mesh_for_linear_triangle_surface_region(_sr, sm, _tol);

	// generate unit mesh of the domain
    generate_final_parametric_mesh_for_linear_triangle_surface_region(_sr, 
		sm, _tol);

	// all_surface_region_mesh_.push_back(sm);

	return;
}


// functions below are for ferguson surface
void SurfaceMesher::generate_boundary_mesh_for_ferguson_surface_region()
{
	int i;
	DTIso2D generator;
	assert(generator.crtEnv());
	generator.initBouInfo(one_initial_front_param_point_);
	generator.bndPntInst();
	generator.recoverBnds();
	generator.clrOuterEles();
	generator.recDistPnts();

	// for (i=0; i<2; ++i)
		// generator.innerPntInst();
	// generator.smooth();

	generator.rmvEmpNods();
	generator.rmvEmpEles();	
	generator.updateBndParent();

	generator.output();
	generator.writePl2("boundary_param_mesh.pl2");
	generator.writeOff("boundary_param_mesh.off");

	// clear one_initial_front_param_point_
	for (i=0; i<one_initial_front_param_point_.size(); ++i)
	{
		one_initial_front_param_point_[i]->clear();
        delete one_initial_front_param_point_[i];
	}
	one_initial_front_param_point_.clear();

	// write into _sm;

	return;
}

void SurfaceMesher::riemannian_length_2d_line_for_ferguson_surface_region(SurfaceRegion *_sr, Vec3d& _pt1, Vec3d& _pt2, double &_rl)
{
	int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);
    
	Vec2d uv;
	Vec3d ru, rv, pt;
	double h;
	Mat2by2d mt;
	Mat1by1d cons;
	Mat2by1d vec(_pt2[0]-_pt1[0], _pt2[1]-_pt1[1]); // todo_clg: originally vec(_pt2[1]-_pt1[1], _pt2[0]-_pt1[0]), seriously wrong!!!

	uv[0] = _pt1[0];
	uv[1] = _pt1[1];
	fs->evaluate(uv, pt);
	h = mesh_size_spec_->get_final_spacing(pt);
	fs->calculate_first_derivative(uv, ru, rv);
	mt = Mat2by2d(ru|ru, ru|rv, ru|rv, rv|rv);
	cons = ~vec * mt * vec;

	_rl = sqrt(cons(0, 0)) / h;

    uv[0] = _pt2[0];
	uv[1] = _pt2[1];
	fs->evaluate(uv, pt);
	h = mesh_size_spec_->get_final_spacing(pt);
	fs->calculate_first_derivative(uv, ru, rv);
	mt = Mat2by2d(ru|ru, ru|rv, ru|rv, rv|rv);
	cons = ~vec * mt * vec;

	_rl += sqrt(cons(0, 0)) / h;

	// added on 20080110-1603, use adaptive simpson formula to calculate the riemannian length of the mesh edge
	uv[0] = (_pt1[0] + _pt2[0]) / 2.0;
	uv[1] = (_pt1[1] + _pt2[1]) / 2.0;
	fs->evaluate(uv, pt);
	h = mesh_size_spec_->get_final_spacing(pt);
	fs->calculate_first_derivative(uv, ru, rv);
	mt = Mat2by2d(ru|ru, ru|rv, ru|rv, rv|rv);
	cons = ~vec * mt * vec;

	_rl += 4.0 * sqrt(cons(0, 0)) / h;

	_rl /= 6.0;

	return;
}

void SurfaceMesher::create_binary_tree_for_ferguson_surface_region(SurfaceRegion *_sr, MidsplitedEdge *_mse)
{
    const double threshold_length = 0.1; // _tol!
	if (_mse->riemannian_length <= threshold_length)  // recursive termination condition
		return;

	_mse->left = new MidsplitedEdge;
	_mse->right = new MidsplitedEdge;

	_mse->left->point[0] = _mse->point[0];
	_mse->left->point[1] = (_mse->point[0]+_mse->point[1])/2.0;
	_mse->left->left = _mse->left->right = NULL;
	_mse->left->parent = _mse;
	riemannian_length_2d_line_for_ferguson_surface_region(_sr, _mse->left->point[0], _mse->left->point[1], _mse->left->riemannian_length);

    _mse->right->point[0] = (_mse->point[0]+_mse->point[1])/2.0;
	_mse->right->point[1] = _mse->point[1];
	_mse->right->left = _mse->right->right = NULL;
	_mse->right->parent = _mse;
	riemannian_length_2d_line_for_ferguson_surface_region(_sr, _mse->right->point[0], _mse->right->point[1], _mse->right->riemannian_length);

	create_binary_tree_for_ferguson_surface_region(_sr, _mse->left);
	create_binary_tree_for_ferguson_surface_region(_sr, _mse->right);
}

void SurfaceMesher::discretize_2d_mesh_inner_edges_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh, double _tol)
{
    int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);

	HalfedgeHandle heh0, heh1;
	VertexHandle   vh0, vh1;
	EdgeHandle     eh; // for debug
	Vec3d pt(0.0, 0.0, 0.0);
	Vec2d uv(0.0, 0.0);
	double total_riemannian_length = 0.0, ideal_riemannian_length = 0.0;
	int i, ideal_edge_number;
	MeshPoint3D mp3d;

	const int NP_THRE = 10;
	static int times = 0;
	static SurfaceRegion *sr_old = NULL;
	if (0==times)
		sr_old = _sr;
	if (sr_old!=_sr)
	{
		sr_old = _sr;
		times = 0;
	}

	std::vector<MidsplitedPoint> midsplited_point;
	MidsplitedPoint msp;
	std::list<MeshPoint3D> *lst_mp3d = NULL;

	MyMesh::EdgeIter e_it(_mesh.edges_begin()),
		             e_end(_mesh.edges_end());

	// compute the entire 2d riemmanian length and the discretized points on inner mesh edges
    for (; e_it!=e_end; ++e_it)
	{
		eh = e_it.handle();
        heh0   = _mesh.halfedge_handle(e_it.handle(), 0);
	    heh1   = _mesh.halfedge_handle(e_it.handle(), 1);
        vh0    = _mesh.to_vertex_handle(heh0);
	    vh1    = _mesh.to_vertex_handle(heh1);

		// maybe we should compute the riemannian length of the boundary edge for debugging purpose
		if (_mesh.is_boundary(e_it))
		{
            _mesh.property(discretized_point, e_it) = new std::list<MeshPoint3D>;
		    lst_mp3d = _mesh.property(discretized_point, e_it);
            mp3d.point = _mesh.point(vh0);
			mp3d.isotropic_size = _mesh.property(iso3d_size, vh0);
			lst_mp3d->push_back(mp3d);
			mp3d.point = _mesh.point(vh1);
			mp3d.isotropic_size = _mesh.property(iso3d_size, vh1);
			lst_mp3d->push_back(mp3d);

			continue; // todo_clg: commented temporarily for debugging purpose, check the riemannian length of the boundary edges
		}

		// codes below are for inner edges
        // create binary tree to generate intermediate auxiliary points
        MidsplitedEdge *root = new MidsplitedEdge; // where to free memory?!
	    root->point[0] = _mesh.point(vh0);
	    root->point[1] = _mesh.point(vh1);
	    root->left = root->right = NULL;
	    root->parent = NULL;
	    riemannian_length_2d_line_for_ferguson_surface_region(_sr, root->point[0], root->point[1], root->riemannian_length);
        create_binary_tree_for_ferguson_surface_region(_sr, root); 

	    // store the intermediate auxiliary points as a vector
	    // std::vector<MidsplitedPoint> midsplited_point;
	    // MidsplitedPoint msp;
	    msp.point = root->point[0];
	    msp.riemannian_length = 0.0;
	    midsplited_point.push_back(msp);
	    // may need to output the results to test the correctness
		// this function is shared with STL model
	    traverse_to_get_midsplited_point_2(root, midsplited_point);


		// determine the total riemannian length, ideal edge number, ideal riemannian length
        ideal_riemannian_length = 0.0;
		total_riemannian_length = 0.0;
	    ideal_edge_number = 1;
	    for (i=0; i<(int)midsplited_point.size(); ++i)
	    {
		    total_riemannian_length += midsplited_point[i].riemannian_length;
			uv[0] = midsplited_point[i].point[0];
			uv[1] = midsplited_point[i].point[1];
			fs->evaluate(uv, pt);
		    midsplited_point[i].isotropic_size = mesh_size_spec_->get_final_spacing(pt);
		    if (i > 0)
			    midsplited_point[i].riemannian_length += midsplited_point[i-1].riemannian_length;
	    }
        ideal_edge_number = (int)total_riemannian_length;
	    if (total_riemannian_length - ideal_edge_number >= 0.5)
	     	ideal_edge_number++;
	    if (0==ideal_edge_number)
		    ideal_edge_number = 1;
	    ideal_riemannian_length = total_riemannian_length / (1.0 * ideal_edge_number);

	    cout << "    Number of final discretized points: " << ideal_edge_number+1 << endl;

	    // generate final discretized points
		// actually it's MeshPoint2D! because of OpenMesh, we have to use 3D instead of 2D!
	    // std::list<MeshPoint3D> *lst_mp3d = NULL; // = new std::list<MeshPoint3D>; do not have "(ideal_edge_number+1)"!
		_mesh.property(discretized_point, e_it) = new std::list<MeshPoint3D>;
		lst_mp3d = _mesh.property(discretized_point, e_it);
		mp3d.point = midsplited_point[0].point;
		mp3d.isotropic_size = midsplited_point[0].isotropic_size;
		lst_mp3d->push_back(mp3d);

		if (ideal_edge_number>NP_THRE && times<2) 
		{ 
			// only take the middle point (three inner points)
			// it's necessary, some problems could happen
			mp3d.point = (midsplited_point[0].point + midsplited_point[midsplited_point.size()-1].point) / 2.0;
			uv[0] = mp3d.point[0];
			uv[1] = mp3d.point[1];
			fs->evaluate(uv, pt);
			mp3d.isotropic_size = mesh_size_spec_->get_final_spacing(pt);
			lst_mp3d->push_back(mp3d);
		}
		else
		{
			int j, k = 0, m = 0;
			double interval, c, f0, fc, temp_ideal_riemannian_length = 0.0, temp_mid_riemannian_length = 0.0;
			const double threshold_length = 1.e-6; // _tol!
			for (i=1; i<=ideal_edge_number-1; ++i)
			{
				for (j=k; j<(int)midsplited_point.size()-1; ++j)
				{
					if (midsplited_point[j].riemannian_length < i * ideal_riemannian_length 
						&& i * ideal_riemannian_length <= midsplited_point[j+1].riemannian_length)
					{
						k = j;
						break;
					}
				}
				assert(midsplited_point.size()-1 != j); // find the index

				// get the final discretized point by recursion
				// attension and todo_clg: the efficiency should be quite low because of the 
				// recursion, it need not to be such precise, it should be treated more roughly
				if (fabs(i*ideal_riemannian_length-midsplited_point[j].riemannian_length) < threshold_length)
				{
					mp3d.point = midsplited_point[j].point;
					mp3d.isotropic_size = midsplited_point[j].isotropic_size;
					lst_mp3d->push_back(mp3d);
					// lst_mp3d->at(i).point = midsplited_point[j].point;
					// lst_mp3d->at(i).isotropic_size = midsplited_point[j].isotropic_size;
					continue;
				}
				else if (fabs(midsplited_point[j+1].riemannian_length-i*ideal_riemannian_length) < threshold_length)
				{
					mp3d.point = midsplited_point[j].point;
					mp3d.isotropic_size = midsplited_point[j].isotropic_size;
					lst_mp3d->push_back(mp3d);
					// lst_mp3d->at(i).point = midsplited_point[j+1].point;
					// lst_mp3d->at(i).isotropic_size = midsplited_point[j+1].isotropic_size;
					continue;
				}

				temp_ideal_riemannian_length = i * ideal_riemannian_length - midsplited_point[j].riemannian_length;
				c = interval = 0.5;
				f0 = midsplited_point[j+1].riemannian_length - midsplited_point[j].riemannian_length 
					- 0.5*(midsplited_point[j].isotropic_size - midsplited_point[j+1].isotropic_size);
				fc = f0*c+0.5*c*c*(midsplited_point[j].isotropic_size - midsplited_point[j+1].isotropic_size);
				m = 1;
				while(fabs(temp_ideal_riemannian_length-fc) > threshold_length && m <= 20) // originally ||
				{
					interval *= 0.5;
					if (temp_ideal_riemannian_length < fc)
						c -= interval; 
					else
						c += interval;
					fc = f0*c+0.5*c*c*(midsplited_point[j].isotropic_size - midsplited_point[j+1].isotropic_size);
					m++;
				}

				if (fabs(temp_ideal_riemannian_length-fc) > threshold_length && m > 20)
				{
					cout << "discretize_2d_mesh_inner_edges(): iteration exceeded!" << endl;
					assert(false);
					return; // exit(0)? return false?
				}

				mp3d.point = midsplited_point[j].point + c*(midsplited_point[j+1].point - midsplited_point[j].point);
				uv[0] = mp3d.point[0];
				uv[1] = mp3d.point[1];
				fs->evaluate(uv, pt);
				mp3d.isotropic_size = mesh_size_spec_->get_final_spacing(pt);
				lst_mp3d->push_back(mp3d);

			}
		}

		mp3d.point = midsplited_point[midsplited_point.size()-1].point;
        mp3d.isotropic_size = midsplited_point[midsplited_point.size()-1].isotropic_size;
	    lst_mp3d->push_back(mp3d);

 	    midsplited_point.clear();
	}

	++times;

	return;
}

void SurfaceMesher::filter_param_mesh_inner_edge_discretized_point_for_ferguson_surface_region(SurfaceRegion *_sr, 
		MyMesh &_mesh, double _tol)
{
	int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);

	VertexHandle   vh; // for debug
	Vec3d pt0, pt1, pt;
	int sz0, sz1;
	double rl;
	const double min_riemannian_dist = 0.70710678;

    MyMesh::ConstVertexIter  v_it(_mesh.vertices_begin()),
			  		         v_end(_mesh.vertices_end());
	std::vector<MyMesh::EdgeHandle> vec_unit_edge; // no inserting points, include boundary edges
	std::vector<MyMesh::EdgeHandle> vec_non_unit_edge;
	std::vector<MyMesh::EdgeHandle>::iterator vec_unit_edge_iter, vec_non_unit_edge_iter, vec_non_unit_edge_iter1;
	MyMesh::EdgeHandle eh;
    MyMesh::ConstVertexOHalfedgeIter cvoh_it, temp_cvoh_it;
	std::list<MeshPoint3D> *lst_mp3d0 = NULL, *lst_mp3d1 = NULL;
	std::list<MeshPoint3D>::iterator lst_iter0, lst_iter1, lst_iter00, lst_iter11;

	for (; v_it!=v_end; ++v_it)
	{
		vh = v_it.handle();
		pt = _mesh.point(v_it.handle());

        for (cvoh_it=_mesh.cvoh_iter(v_it.handle()); cvoh_it; ++cvoh_it)
		{
			eh = _mesh.edge_handle(cvoh_it.handle());
            lst_mp3d0 = _mesh.property(discretized_point, eh);

			if (is_equal(pt, lst_mp3d0->back().point))
				lst_mp3d0->reverse();
			else if (!is_equal(pt, lst_mp3d0->front().point)) // for debug
				assert(false);

			if (2==lst_mp3d0->size())
				vec_unit_edge.push_back(eh);
			else
				vec_non_unit_edge.push_back(eh);
		}

        // first compare with inserted end-points of edges in vec_unit_edge, be aware of efficiency
        for (vec_non_unit_edge_iter=vec_non_unit_edge.begin(); vec_non_unit_edge_iter!=vec_non_unit_edge.end(); 
			++vec_non_unit_edge_iter)
		{
			for (vec_unit_edge_iter=vec_unit_edge.begin(); vec_unit_edge_iter!=vec_unit_edge.end();
				++vec_unit_edge_iter)
			{
				lst_mp3d0 = _mesh.property(discretized_point, *vec_unit_edge_iter);
				pt0 = lst_mp3d0->back().point; // here only compare with the other end point

                lst_mp3d1 = _mesh.property(discretized_point, *vec_non_unit_edge_iter);
                lst_iter0 = lst_mp3d1->begin();
			    lst_iter1 = lst_mp3d1->end();
			    ++lst_iter0;
			    --lst_iter1;
			    for (; lst_iter0!=lst_iter1; )
			    {
					pt1 = lst_iter0->point;
                    riemannian_length_2d_line_for_ferguson_surface_region(_sr, pt0, pt1, rl);
					if (rl<min_riemannian_dist)
				  		lst_iter0 = lst_mp3d1->erase(lst_iter0);
				    else
						break;
			    }
			}
		}

		// then compare among edges with more than 2 discretized points
        for (vec_non_unit_edge_iter=vec_non_unit_edge.begin(); vec_non_unit_edge_iter!=vec_non_unit_edge.end(); 
			++vec_non_unit_edge_iter)
		{
			lst_mp3d0 = _mesh.property(discretized_point, *vec_non_unit_edge_iter);
			lst_iter1 = lst_mp3d0->end();
			--lst_iter1;

			vec_non_unit_edge_iter1 = vec_non_unit_edge_iter + 1;
			for (; vec_non_unit_edge_iter1!=vec_non_unit_edge.end(); ++vec_non_unit_edge_iter1)
			{
				lst_iter0 = lst_mp3d0->begin();
				++lst_iter0;

				lst_mp3d1 = _mesh.property(discretized_point, *vec_non_unit_edge_iter1);
                
				for (; lst_iter0!=lst_iter1; ++lst_iter0)
				{
					pt0 = (*lst_iter0).point;

					lst_iter00 = lst_mp3d1->begin();
				    ++lst_iter00;
				    lst_iter11 = lst_mp3d1->end();
				    --lst_iter11;
					for (; lst_iter00!=lst_iter11; )
					{
                        pt1 = (*lst_iter00).point;
                        riemannian_length_2d_line_for_ferguson_surface_region(_sr, pt0, pt1, rl);
                        if (rl<min_riemannian_dist)
				  		    lst_iter00 = lst_mp3d1->erase(lst_iter00);
						else
							break;
					}
				}
			}
		}

		vec_unit_edge.clear();
		vec_non_unit_edge.clear();
	}


	MyMesh::ConstEdgeIter e_it(_mesh.edges_begin()),
		                  e_end(_mesh.edges_end());
	for (; e_it!=e_end; ++e_it)
	{
		if (_mesh.is_boundary(e_it))
			continue;

        lst_mp3d0 = _mesh.property(discretized_point, e_it);
		cout << "No. of discretized points after filtering: " << lst_mp3d0->size() << endl;
	}

	return;
}

void SurfaceMesher::calc_riemannian_circumcenter_for_ferguson_surface_region(SurfaceRegion *_sr, 
		MyMesh &_mesh, MyMesh::FaceHandle &_fh, Vec3d &_cc)
{
    int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);

	Vec2d uv;
	Vec3d ru, rv;
	Mat2by2d mt[3];
	Vec3d pt[3];
	double spa;
	int i;
	MyMesh::ConstFaceVertexIter cfv_it;
	for (cfv_it=_mesh.cfv_iter(_fh), i=0; cfv_it; ++cfv_it, ++i)
	{
		pt[i] = _mesh.point(cfv_it);
		uv = Vec2d(pt[i][0], pt[i][1]);
		fs->calculate_first_derivative(uv, ru, rv);
        spa = _mesh.property(iso3d_size, cfv_it.handle());
		mt[i] = Mat2by2d(ru|ru, ru|rv, ru|rv, rv|rv);
		mt[i] = mt[i] / (spa*spa);
	}

	// initialization
    _cc = (pt[0] + pt[1] + pt[2]) / 3.0; 
	Mat2by1d fx, du;
	Mat2by2d dfx;
	Mat1by1d val;
	bool succ = false;
	const double tol = 1.0e-6;
	const int max_iter_num = 100;
	for (i=0; i<max_iter_num; ++i)
	{
		// val = ~du * mt[0] * du;
	    // fx(0, 0) = val(0, 0); // should provide a type transformation from Mat1by1d to double

		// use newton method to solve the nonlinear equation set
        du = Mat2by1d(pt[0][0]-_cc[0], pt[0][1]-_cc[1]);
		val = ~du * mt[0] * du;
        fx(0, 0) = fx(1, 0) = val(0,0);
        du = Mat2by1d(pt[1][0]-_cc[0], pt[1][1]-_cc[1]);
		val = ~du * mt[1] * du;
        fx(0, 0) -= val(0,0);
        du = Mat2by1d(pt[2][0]-_cc[0], pt[2][1]-_cc[1]);
		val = ~du * mt[2] * du;
        fx(1, 0) -= val(0,0);

		dfx(0, 0) = 2*(mt[0](0,0)-mt[1](0,0))*_cc[0] + 2*(mt[0](0,1)-mt[1](0,1))*_cc[1]
		    - 2*(mt[0](0,0)*pt[0][0]-mt[1](0,0)*pt[1][0]) - 2*(mt[0](0,1)*pt[0][1]-mt[1](0,1)*pt[1][1]);
		dfx(0, 1) = 2*(mt[0](1,1)-mt[1](1,1))*_cc[1] + 2*(mt[0](0,1)-mt[1](0,1))*_cc[0]
		    - 2*(mt[0](1,1)*pt[0][1]-mt[1](1,1)*pt[1][1]) - 2*(mt[0](0,1)*pt[0][0]-mt[1](0,1)*pt[1][0]);
		dfx(1, 0) = 2*(mt[0](0,0)-mt[2](0,0))*_cc[0] + 2*(mt[0](0,1)-mt[2](0,1))*_cc[1]
		    - 2*(mt[0](0,0)*pt[0][0]-mt[2](0,0)*pt[2][0]) - 2*(mt[0](0,1)*pt[0][1]-mt[2](0,1)*pt[2][1]);
		dfx(1, 1) = 2*(mt[0](1,1)-mt[2](1,1))*_cc[1] + 2*(mt[0](0,1)-mt[2](0,1))*_cc[0]
		    - 2*(mt[0](1,1)*pt[0][1]-mt[2](1,1)*pt[2][1]) - 2*(mt[0](0,1)*pt[0][0]-mt[2](0,1)*pt[2][0]);

		du = dfx.inversion() * fx;
        _cc[0] -= du(0, 0);
		_cc[1] -= du(1, 0);

		if (sqrt(du(0,0)*du(0,0)+du(1,0)*du(1,0))<tol)
		{
			succ = true;
			break;
		}
	}
	if (succ) // check if in parametric space
	{
		// can be not in the parametric space, but it seems ok
        // assert(_cc[0]>=0.0 && _cc[0]<=1.0 && _cc[1]>=0.0 && _cc[1]<=1.0);
        // assert(_cc[0]+_cc[1]-1.0<=0.0);
	}
	else
	{
	    // assert(false); // todo_clg: should be valid
	    cout << "calc_riemannian_circumcenter: maximum number of iterations exceeded" << endl;
	}

	return;
}

bool SurfaceMesher::is_delaunay_broken_for_ferguson_surface_region(SurfaceRegion *_sr, 
		MyMesh &_mesh, MyMesh::VertexHandle &_vh, MyMesh::FaceHandle &_fh)
{
    Vec3d cc;
    calc_riemannian_circumcenter_for_ferguson_surface_region(_sr, _mesh, _fh, cc);

    int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);
    // Mat2by2d metric_tensor = lts->get_metric_tensor();
	Mat2by2d mt;
	Vec2d uv;
	Vec3d ru, rv;
    Vec3d pt = _mesh.point(_vh), pt0 = _mesh.point(_mesh.cfv_iter(_fh).handle());
	double sz;
	Mat2by1d op = Mat2by1d(pt[0]-cc[0], pt[1]-cc[1]), op0 = Mat2by1d(pt0[0]-cc[0], pt0[1]-cc[1]);
	Mat1by1d val;
	double rl0, rl1, dm = 0.0; // delaunay measure
	
	// first deal with the inserting point
	uv = Vec2d(pt[0], pt[1]);
	fs->calculate_first_derivative(uv, ru, rv);
	mt = Mat2by2d(ru|ru, ru|rv, ru|rv, rv|rv);
	sz = _mesh.property(iso3d_size, _vh);
	mt = mt / (sz*sz);
	val = ~op * mt * op;
    rl0 = sqrt(val(0,0));
	val = ~op0 * mt * op0;
	rl1 = sqrt(val(0,0));
	dm += rl0/rl1;

	// then three face end-points
	MyMesh::ConstFaceVertexIter cfv_it;
	for (cfv_it=_mesh.cfv_iter(_fh); cfv_it; ++cfv_it)
	{
        pt0 = _mesh.point(cfv_it.handle());
		uv = Vec2d(pt0[0], pt0[1]);
		fs->calculate_first_derivative(uv, ru, rv);
		mt = Mat2by2d(ru|ru, ru|rv, ru|rv, rv|rv);
        sz = _mesh.property(iso3d_size, cfv_it.handle());
		mt = mt / (sz*sz);
		val = ~op * mt * op;
		rl0 = sqrt(val(0, 0));
		val = ~op0 * mt * op0;
	    rl1 = sqrt(val(0, 0));
	    dm += rl0/rl1;
	}

	if (dm<4)
		return true;
	else
	    return false;
}

void SurfaceMesher::insert_one_point_for_ferguson_surface_region(SurfaceRegion *_sr, MeshPoint3D &_mp3d, MyMesh &_mesh)
{
	MyMesh::FaceIter f_it(_mesh.faces_begin()),
		             f_end(_mesh.faces_end());
	std::vector<FaceHandle> vec_base_fh;
	FaceHandle fh;
	Vec3d pt = _mp3d.point, pt0, pt1;
	MyMesh::FaceHalfedgeIter fh_it;
	double area;
	int cnt0;
	MyMesh::FaceFaceIter ff_it;
	bool bval;
	EdgeHandle eh, eh0, eh1;
	std::deque<EdgeHandle> deq_eh;
	HalfedgeHandle heh;
	VertexHandle vh, vh0;

	double rl;
	const double RL_THRE = 0.1;

    // first find one of base elements
    for (; f_it!=f_end; ++f_it)
	{
		cnt0 = 0;

        for (fh_it=_mesh.fh_iter(f_it.handle()); fh_it; ++fh_it)
		{
            pt0 = _mesh.point(_mesh.from_vertex_handle(fh_it.handle()));
			pt1 = _mesh.point(_mesh.to_vertex_handle(fh_it.handle()));
            area2(pt, pt0, pt1, area);

			if (area<-1.0e-8) // todo_clg: originally -1.0e-12, but some examples failed!; originally 0.0, invalid elements generated
				break; // use -1.0e-8 to guarantee find the right base elements

			cnt0++;
		}

		if (cnt0==3)
			vec_base_fh.push_back(f_it.handle());
	}
	cout << "No. of found enclosing base elements: " << vec_base_fh.size() << endl;
	// assert(vec_base_fh.size()==1 || vec_base_fh.size()==2); // todo_clg: should be valid
	if (1!=vec_base_fh.size() && 2!=vec_base_fh.size())
		return;

	bval = false;
	if (2==vec_base_fh.size()) // two base faces must be adjacent
	{
        for (ff_it=_mesh.ff_iter(vec_base_fh[0]); ff_it; ++ff_it)
			if (ff_it.handle()==vec_base_fh[1])
			{
                bval = true;
				break;
			}

		assert(bval);
	}
	
	// start from the base element
    if (2==vec_base_fh.size())
	{
        for (fh_it=_mesh.fh_iter(vec_base_fh[0]); fh_it; ++fh_it)
		{
            fh = _mesh.face_handle(_mesh.opposite_halfedge_handle(fh_it.handle()));
			if (fh==vec_base_fh[1])
				break;
		}
		assert(fh_it);
		eh = _mesh.edge_handle(fh_it.handle()); // shared edge

		// added on 20080108-2158, check if the new edges have appropriate length
        heh = _mesh.halfedge_handle(eh, 0);
		pt0 = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh)));
		riemannian_length_2d_line_for_ferguson_surface_region(_sr, pt, pt0, rl);
		if (rl<RL_THRE)
			return;

		heh = _mesh.halfedge_handle(eh, 1);
		pt1 = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh)));
		riemannian_length_2d_line_for_ferguson_surface_region(_sr, pt, pt1, rl);
		if (rl<RL_THRE)
			return;

        heh = _mesh.halfedge_handle(eh, 0);
		deq_eh.push_back(_mesh.edge_handle(_mesh.next_halfedge_handle(heh)));
		deq_eh.push_back(_mesh.edge_handle(_mesh.prev_halfedge_handle(heh)));
		heh = _mesh.halfedge_handle(eh, 1);
        deq_eh.push_back(_mesh.edge_handle(_mesh.next_halfedge_handle(heh)));
		deq_eh.push_back(_mesh.edge_handle(_mesh.prev_halfedge_handle(heh)));

		// are the previous edge handles above still valid after these kind of
		// topological transformations (and add new vertex) without applying 
		// garbage_collection() ?
		vh = _mesh.add_vertex(pt);
		_mesh.split(eh, vh);
		_mesh.property(iso3d_size, vh) = _mp3d.isotropic_size;
        
        while (!deq_eh.empty())
		{
			eh = deq_eh.front(); // and vh
			deq_eh.pop_front();
			heh = _mesh.halfedge_handle(eh, 0);
            vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));
			if (vh0!=vh)
			{
				fh = _mesh.face_handle(heh);
				eh0 = _mesh.edge_handle(_mesh.prev_halfedge_handle(heh));
				eh1 = _mesh.edge_handle(_mesh.next_halfedge_handle(heh));

				// heh = _mesh.halfedge_handle(eh, 1);
                // vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));
				// assert(vh0==vh);
			}
			else
			{
                heh = _mesh.halfedge_handle(eh, 1);
				eh0 = _mesh.edge_handle(_mesh.prev_halfedge_handle(heh));
				eh1 = _mesh.edge_handle(_mesh.next_halfedge_handle(heh));
                // vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));

				fh = _mesh.face_handle(heh);
			}

			// two things: is flip legal? is delaunay criteria broken?
			// the function is_flip_legal() is shared with STL model
            if (is_flip_legal(_mesh, eh, vh) && is_delaunay_broken_for_ferguson_surface_region(_sr, _mesh, vh, fh))
			{
                _mesh.flip(eh);

				// add new edges, be careful: heh has been changed!
				deq_eh.push_back(eh0);
				deq_eh.push_back(eh1);
			}
		}
	}
	else // 1==vec_base_fh.size()
	{
        for (fh_it=_mesh.fh_iter(vec_base_fh[0]); fh_it; ++fh_it)
            deq_eh.push_back(_mesh.edge_handle(fh_it.handle()));

        vh = _mesh.add_vertex(pt);
		_mesh.split(vec_base_fh[0], vh);
		_mesh.property(iso3d_size, vh) = _mp3d.isotropic_size;

		while (!deq_eh.empty())
		{
			eh = deq_eh.front(); // and vh
			deq_eh.pop_front();
			heh = _mesh.halfedge_handle(eh, 0);
            vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));
			if (vh0!=vh)
			{
				fh = _mesh.face_handle(heh);
				eh0 = _mesh.edge_handle(_mesh.prev_halfedge_handle(heh));
				eh1 = _mesh.edge_handle(_mesh.next_halfedge_handle(heh));

				// heh = _mesh.halfedge_handle(eh, 1);
                // vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));
				// assert(vh0==vh);
			}
			else
			{
                heh = _mesh.halfedge_handle(eh, 1);
				eh0 = _mesh.edge_handle(_mesh.prev_halfedge_handle(heh));
				eh1 = _mesh.edge_handle(_mesh.next_halfedge_handle(heh));
                // vh0 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));

				fh = _mesh.face_handle(heh);
			}

			// two things: is flip legal? is delaunay criteria broken?
            if (is_flip_legal(_mesh, eh, vh) && is_delaunay_broken_for_ferguson_surface_region(_sr, _mesh, vh, fh))
			{
                _mesh.flip(eh);

				// add new edges, be careful: heh has been changed!
				deq_eh.push_back(eh0);
				deq_eh.push_back(eh1);
			}
		}
	}

	_mesh.garbage_collection();

	// assert(OpenMesh::IO::write_mesh(_mesh, "result_param_mesh.off"));

	// MyMesh mesh1;
	// map_mesh_for_ferguson_surface_region(_sr, _mesh, mesh1);
	// assert(OpenMesh::IO::write_mesh(mesh1, "result_physical_mesh.off"));

	return;
}

void SurfaceMesher::insert_filtered_discretized_point_for_ferguson_surface_region(SurfaceRegion *_sr, 
		MyMesh &_mesh, double _tol)
{
	int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);

	std::list<MeshPoint3D> inserting_point;
	std::list<MeshPoint3D> *lst_mp3d = NULL;
	std::list<MeshPoint3D>::iterator lst_iter, lst_iter0;
	MeshPoint3D mp3d;
	MyMesh::ConstEdgeIter e_it(_mesh.edges_begin()),
		                  e_end(_mesh.edges_end());

	for (; e_it!=e_end; ++e_it)
	{
		if (_mesh.is_boundary(e_it))
			continue;

        lst_mp3d = _mesh.property(discretized_point, e_it);
		if (lst_mp3d->size()==2)
			continue;
        
		lst_iter0 = lst_mp3d->end();
		lst_iter0--;
		lst_iter = lst_mp3d->begin();
		lst_iter++;
		for (; lst_iter!=lst_iter0; ++lst_iter)
			inserting_point.push_back(*lst_iter);
	}
	cout << "No. of inserting points: " << inserting_point.size() << endl;

	// insert new generated/discretized points
    for (lst_iter=inserting_point.begin(); lst_iter!=inserting_point.end(); ++lst_iter)
        insert_one_point_for_ferguson_surface_region(_sr, *lst_iter, _mesh);

	return;
}

void SurfaceMesher::calc_element_quality_for_ferguson_surface_region(SurfaceRegion *_sr, 
		MyMesh &_mesh, MyMesh::VertexHandle _vh[3], double &_qua)
{
    int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);

	int i, j;
	Mat2by1d edge[3];
	Vec3d vec[3];
	vec[0] = _mesh.point(_vh[1]) - _mesh.point(_vh[0]), 
	vec[1] = _mesh.point(_vh[2]) - _mesh.point(_vh[0]),
	vec[2] = _mesh.point(_vh[2]) - _mesh.point(_vh[1]);
	for (i=0; i<3; ++i)
		for (j=0; j<2; ++j)
			edge[i](j, 0) = vec[i][j];

	double det0 = Mat2by2d(vec[0][0], vec[1][0], vec[0][1], vec[1][1]).determinant();
	// Mat2by2d mtp, mtg = lts->get_metric_tensor();
	Mat2by2d mt;
	Vec2d uv;
	Vec3d pt, ru, rv;
	_qua = 1.0e+30;
	double sz, numerator, denominator;
	const double sqrt3 = 1.73205080757;
	Mat1by1d rl;
	for (i=0; i<3; ++i)
	{
        sz = _mesh.property(iso3d_size, _vh[i]);
        pt = _mesh.point(_vh[i]);
        uv = Vec2d(pt[0], pt[1]);
		fs->calculate_first_derivative(uv, ru, rv);
		mt = Mat2by2d(ru|ru, ru|rv, ru|rv, rv|rv);
        mt = mt / (sz*sz);
        numerator = 2.0*sqrt3*fabs(sqrt(mt.determinant())*det0);

		denominator = 0.0;
        for (j=0; j<3; ++j)
		{
			rl = ~edge[j] * mt * edge[j];
			denominator += rl(0, 0);
		}

		if (numerator/denominator<_qua)
			_qua = numerator / denominator;
	}

	return;
}

void SurfaceMesher::optimize_with_diagonal_swapping_given_threshold_for_ferguson_surface_region(SurfaceRegion *_sr, 
		MyMesh &_mesh, double _qua_ratio_threshold)
{
    MyMesh::EdgeIter e_it;
    double qua_ratio, qua0, qua1, qua;
	MyMesh::VertexHandle vh[4], vh_param[3];
	MyMesh::HalfedgeHandle heh0, heh1;
	MyMesh::EdgeHandle eh;
	
	for (e_it=_mesh.edges_begin(); e_it!=_mesh.edges_end(); ++e_it)
	{
		if (_mesh.is_boundary(e_it.handle()))
			continue;
        
		eh     = e_it.handle();
        heh0   = _mesh.halfedge_handle(eh, 0);
		heh1   = _mesh.halfedge_handle(eh, 1);
		vh[0]  = _mesh.to_vertex_handle(heh0);
        vh[1]  = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh0));
		vh[2]  = _mesh.to_vertex_handle(heh1);
		vh[3]  = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh1));
        
        vh_param[0] = vh[0];
		vh_param[1] = vh[1];
        vh_param[2] = vh[2];
		calc_element_quality_for_ferguson_surface_region(_sr, _mesh, vh_param, qua0);   
        vh_param[0] = vh[0];
		vh_param[1] = vh[2];
        vh_param[2] = vh[3];
		calc_element_quality_for_ferguson_surface_region(_sr, _mesh, vh_param, qua);
		if (qua<qua0)
			qua0 = qua;

        vh_param[0] = vh[0];
		vh_param[1] = vh[1];
        vh_param[2] = vh[3];
		calc_element_quality_for_ferguson_surface_region(_sr, _mesh, vh_param, qua1);
		vh_param[0] = vh[1];
		vh_param[1] = vh[2];
        vh_param[2] = vh[3];
		calc_element_quality_for_ferguson_surface_region(_sr, _mesh, vh_param, qua);
		if (qua<qua1)
			qua1 = qua;

		if (qua1/qua0>_qua_ratio_threshold && is_flip_legal(_mesh, eh, vh[1]))
		{
			cout << "optimize with edge swapping ..." << endl;
		    _mesh.flip(eh);
		}
	}
    _mesh.garbage_collection();

	return;
}

void SurfaceMesher::optimize_with_diagonal_swapping_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh)
{
	cout << "optimize mesh with diagonal swapping given threshold 2.0" << endl;
	optimize_with_diagonal_swapping_given_threshold_for_ferguson_surface_region(_sr, _mesh, 2.0);
	cout << "optimize mesh with diagonal swapping given threshold 1.5" << endl;
    optimize_with_diagonal_swapping_given_threshold_for_ferguson_surface_region(_sr, _mesh, 1.5);
	cout << "optimize mesh with diagonal swapping given threshold 1.0" << endl;
	optimize_with_diagonal_swapping_given_threshold_for_ferguson_surface_region(_sr, _mesh, 1.0);

	return;
}

void SurfaceMesher::smooth_with_optimal_shape_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh)
{
    int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	Mat2by2d mat, mt;

	MyMesh::VertexIter v_it;
	MyMesh::VertexOHalfedgeIter voh_it;
	MyMesh::HalfedgeHandle heh;
	Mat2by1d ip[3], p0, p01;
	std::vector<Vec3d> vec_ip;
	Vec2d uv;
	Vec3d ru, rv;
	Vec3d pt, pt_tmp, pt_cen, pt_int, p[2];
	MyMesh::VertexHandle vh[3];
	double sz, d, ar;
	double coef, coef0, coef1;
	int i, np = 3;

	for (v_it=_mesh.vertices_begin(); v_it!=_mesh.vertices_end(); ++v_it)
	{
		if (_mesh.is_boundary(v_it.handle()))
			continue;

		vh[2] = v_it.handle();
		pt_cen = _mesh.point(v_it.handle());

        for (voh_it=_mesh.voh_iter(v_it.handle()); voh_it; ++voh_it)
		{
            heh = _mesh.next_halfedge_handle(voh_it.handle());
			vh[0] = _mesh.from_vertex_handle(heh);
			vh[1] = _mesh.to_vertex_handle(heh);
            p[0] = _mesh.point(vh[0]);
			p[1] = _mesh.point(vh[1]);
			p0(0, 0) = p[0][0];
			p0(1, 0) = p[0][1];
			p01(0, 0) = p[1][0] - p[0][0];
			p01(1, 0) = p[1][1] - p[0][1];

			np = 3;
			pt = Vec3d(0.0, 0.0, 0.0);
            for (i=0; i<3; ++i)
			{
                sz = _mesh.property(iso3d_size, vh[i]);
				pt_tmp = _mesh.point(vh[i]);
				uv = Vec2d(pt_tmp[0], pt_tmp[1]);
                fs->calculate_first_derivative(uv, ru, rv);
				if (is_zero(ru.norm(), 1.0e-5) || is_zero(rv.norm(), 1.e-5))
				{
					--np;
					// ip[i](0, 0) = pt_cen[0];
					// ip[i](1, 0) = pt_cen[1];
					continue;
				}
				mt = Mat2by2d(ru|ru, ru|rv, ru|rv, rv|rv);
				mt = mt / (sz*sz);
                d = sqrt((mt(0, 0)*mt(1, 1) - mt(0, 1)*mt(0, 1)) / 3.0);
				mat(0, 0) = d - mt(0, 1);
				mat(0, 1) = -mt(1, 1);
				mat(1, 0) = mt(0, 0);
				mat(1, 1) = d + mt(0, 1);
				mat /= (2.0*d);

				ip[i] = p0 + mat*p01;
				pt += Vec3d(ip[i](0, 0), ip[i](1, 0), 0);
			}

            // pt[0] = (ip[0](0, 0) + ip[1](0, 0) + ip[2](0, 0)) / 3.0;
			// pt[1] = (ip[0](1, 0) + ip[1](1, 0) + ip[2](1, 0)) / 3.0;
			// pt[2] = 0.0;

			assert(np>0);	
		    pt /= (double)np;

			vec_ip.push_back(pt);
		}

		pt = Vec3d(0.0, 0.0, 0.0);
		for (i=0; i<vec_ip.size(); ++i)
            pt += vec_ip[i];
		pt /= vec_ip.size();
		vec_ip.clear();

		// avoid generating invalid elements
		// compute the intersection points with the coefficient of each bounding edge and 
		// the edge (pt_cen, pt)
		coef = 1.0;
        for (voh_it=_mesh.voh_iter(v_it.handle()); voh_it; ++voh_it)
		{
            heh = _mesh.next_halfedge_handle(voh_it.handle());
			vh[0] = _mesh.from_vertex_handle(heh);
			vh[1] = _mesh.to_vertex_handle(heh);
            p[0] = _mesh.point(vh[0]);
			p[1] = _mesh.point(vh[1]);
            
			// todo_clg, how to deal with this failed situation
			// assert(seg_seg_int(pt_cen, pt, p[0], p[1], coef0, coef1, pt_int)); // function seg_seg_int() is shared with STL model
			if (!seg_seg_int(pt_cen, pt, p[0], p[1], coef0, coef1, pt_int)) // parallel
				continue;  // todo_clg: all right?

			if (coef0>=0.0 && coef0<coef)
				coef = coef0;
		}
		if (coef<1.0)
			coef *= 0.5;

		_mesh.point(v_it.handle()) = pt_cen + coef*(pt-pt_cen);
	}

	return;
}

void SurfaceMesher::map_mesh_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh_src, MyMesh &_mesh_tar)
{
    int surface_id = _sr->surface_id;
	FergusonSurface *fs = (FergusonSurface*)geometry_->get_surface_geometry().at(surface_id-1);
	_mesh_tar.assign_connectivity(_mesh_src); // doesn't copy vertices coordinates

	Vec2d uv;
	Vec3d pt;
	MyMesh::VertexIter v_it_src(_mesh_src.vertices_begin()), v_it_tar(_mesh_tar.vertices_begin());
	
	for (; v_it_src!=_mesh_src.vertices_end() && v_it_tar!=_mesh_tar.vertices_end(); ++v_it_src, ++v_it_tar)
	{
        pt = _mesh_src.point(v_it_src.handle());
		uv = Vec2d(pt[0], pt[1]);
		// uv[0] = pt[0];
		// uv[1] = pt[1];
		fs->evaluate(uv, pt);
		_mesh_tar.point(v_it_tar.handle()) = pt;
	}

	return;
}

void SurfaceMesher::generate_final_parametric_mesh_for_ferguson_surface_region(SurfaceRegion *_sr, SurfaceMesh *_sm, double _tol)
{
	MyMesh mesh, *mesh1 = new MyMesh;

	// OpenMesh::EPropHandleT<double> riemmanian_length; // if put here, it couldn't referenced by other functions, 
	mesh.add_property(riemannian_length);
	// mesh1.add_property(riemmanian_length);

	// OpenMesh::EPropHandleT<std::deque<MeshPoint3D> *> discretized_point; // Vec2d is enough
    mesh.add_property(discretized_point);
	// mesh1.add_property(discretized_point);

	// OpenMesh::VPropHandleT<double> iso3d_size; // set as global variable
	mesh.add_property(iso3d_size);
	// mesh.add_property(iso3d_size);

	mesh1->add_property(global_id);

    assert(OpenMesh::IO::read_mesh(mesh, "boundary_param_mesh.off"));

	// calculate the spacing value on each initial boundary vertex
	Vec3d pt;
	double spa;
	MyMesh::ConstVertexIter v_it(mesh.vertices_begin()),
		                    v_end(mesh.vertices_end());
	for (; v_it!=v_end; ++v_it)
	{
        pt = mesh.point(v_it.handle());
        spa = mesh_size_spec_->get_final_spacing(pt);
		mesh.property(iso3d_size, v_it.handle()) = spa;
	}

	int npt_new = 0, npt_old = mesh.n_vertices();
	while (true)
	{
	    // discretize inner edges of parametric mesh just like discretize line curve
        discretize_2d_mesh_inner_edges_for_ferguson_surface_region(_sr, mesh, _tol);

	    // filter the discretized points on inner mesh edges
        filter_param_mesh_inner_edge_discretized_point_for_ferguson_surface_region(_sr, mesh, _tol);

	    // insert filtered discretized points
        insert_filtered_discretized_point_for_ferguson_surface_region(_sr, mesh, _tol);

		// map_mesh_for_ferguson_surface_region(_sr, mesh, *mesh1);
		// assert(OpenMesh::IO::write_mesh(mesh, "result_param_mesh.off"));
		// assert(OpenMesh::IO::write_mesh(*mesh1, "result_physical_mesh.off"));

		optimize_with_diagonal_swapping_for_ferguson_surface_region(_sr, mesh);

		for (int i=0; i<3; ++i)
		    smooth_with_optimal_shape_for_ferguson_surface_region(_sr, mesh);

        // map_mesh_for_ferguson_surface_region(_sr, mesh, *mesh1);
		// assert(OpenMesh::IO::write_mesh(mesh, "result_param_mesh.off"));
		// assert(OpenMesh::IO::write_mesh(*mesh1, "result_physical_mesh.off"));

		npt_new = mesh.n_vertices();

		if (npt_new==npt_old)
			break;
		else
		    npt_old = npt_new;
	}

	all_surface_region_mesh_.push_back(mesh1);

	map_mesh_for_ferguson_surface_region(_sr, mesh, *mesh1);
	assert(OpenMesh::IO::write_mesh(mesh, "result_param_mesh.off"));
	assert(OpenMesh::IO::write_mesh(*mesh1, "result_physical_mesh.off"));

	return;
}

void SurfaceMesher::dt_for_ferguson_surface_region(SurfaceRegion *_sr, double _tol)
{
    SurfaceMesh *sm = new SurfaceMesh;
	// sm->surface_region_id = 0; // will be assigned outside
	// assumption: only one loop here
	//sm->boundary_point_number = (int)one_initial_front_param_point_[0]->size();

	// because of the round-off error, the discrete boundary points
	// don't exactly on the parametric triangle, which need a minor modification
	// for (i=0; i<(int)one_initial_front_param_point_[0]->size(); ++i)
		// sm->mesh_point_2d.push_back(one_initial_front_param_point_[0]->at(i));

    // start generalized delauney triangulation 
    // generate boundary mesh
     generate_boundary_mesh_for_ferguson_surface_region();

	// generate unit mesh of the domain
    generate_final_parametric_mesh_for_ferguson_surface_region(_sr, sm, _tol);

	return;
}

// the actual dt in parametric space should be unified in one function in the future 
void SurfaceMesher::dt_in_parametric_space(SurfaceRegion *_sr, double _tol)
{
    assert(NULL!=_sr && _tol>0.0);

	SurfaceType surface_type = geometry_->get_surface_geometry().at(_sr->surface_id-1)->get_surface_type();
	switch (surface_type)
	{
	case SurfaceLinearTriangle:
		dt_for_linear_triangle_surface_region(_sr, _tol);
		break;
	case SurfaceQuadraticTriangle:
		// todo_clg
		break;
	case SurfaceFerguson:
		dt_for_ferguson_surface_region(_sr, _tol);
		break;
	case SurfaceNurbs:
		// todo_clg
		break;
	default:
		;
	}

	return;
}

void SurfaceMesher::form_one_surface_mesh(SurfaceRegion *_sr, double _tol)
{


	return;
}

// parameters? Don't have corresponding data structures!
void SurfaceMesher::enhance_one_surface_mesh() 
{


	return;
}

void SurfaceMesher::discretize_one_surface_region(SurfaceRegion *_sr, double _tol)
{
	project_discrete_boundary_point(_sr, _tol);
	dt_in_parametric_space(_sr, _tol);
	form_one_surface_mesh(_sr, _tol);
	// enhance_one_surface_mesh(); // haven't designed parameters yet

	return;
}

void SurfaceMesher::discretize_all_surface_region()
{
	int i;
	for (i=0; i<(int)geometry_->get_surface_region().size(); ++i)
	{
		cout << "SurfaceMesher -> discretizing the " << i+1 << "th surface region ..." << endl;
		discretize_one_surface_region(geometry_->get_surface_region().at(i), 1.e-6);
	}

	return;
}

void SurfaceMesher::form_global_mesh()
{
    // EPS1 is the tolerance for identifying the corner points
	const double EPS0 = 1.E-4, EPS1 = EPS0 * EPS0;
	int np, nco; // neg
	int i, j, k, sz, cur_idx;
	bool existed0, existed1;
	Vec3d v0, v1, v;

	MyMesh::ConstVertexIter cv_it;
	MyMesh::ConstFaceIter cf_it;
	MyMesh::ConstFaceVertexIter cfv_it;
	MyMesh *pmsh = NULL;
	SurfaceRegion *psr = NULL;
	Vec3i elem;

	global_surface_mesh_.surface_region_id = -1;
	global_surface_mesh_.boundary_point_number = 0;
	global_surface_mesh_.mesh_point_3d.clear();
	global_surface_mesh_.mesh_element.clear();

	// Now we deal with the vector containing the points generated on the
    // edges of the surface.

	// First identifies the corner points by looping over the edges, starts the
    // numbering and stores the 3D coordinates
	np = 0;
	for (i=0; i<all_discretized_curve_point_.size(); ++i)
	{
		sz = all_discretized_curve_point_[i]->size();
		v0 = all_discretized_curve_point_[i]->front().point;
		v1 = all_discretized_curve_point_[i]->back().point;
        
		// Checks if the end points have been included in the list
		existed0 = existed1 = false;
		if (global_surface_mesh_.mesh_point_3d.size()>0)
		{
			for (j=0; j<global_surface_mesh_.mesh_point_3d.size(); ++j)
			{
				v = global_surface_mesh_.mesh_point_3d[j];

				if (is_equal(v0, v, EPS0))
				{
                    all_discretized_curve_point_[i]->front().global_id = j;
					existed0 = true;
				}
                if (is_equal(v1, v, EPS0))
				{
                    all_discretized_curve_point_[i]->back().global_id = j;
					existed1 = true;
				}
				if (existed0 && existed1)
					break;
			}
		}
		if (!existed0)
		{
			global_surface_mesh_.mesh_point_3d.push_back(v0);
			++np;
            all_discretized_curve_point_[i]->front().global_id = np - 1;
		}
		if (!existed1)
		{
            global_surface_mesh_.mesh_point_3d.push_back(v1);
			++np;
            all_discretized_curve_point_[i]->back().global_id = np - 1; 
		}
	}
	// Number of corner points
	nco = np;

	// Numbers the rest of the points on the edges
	for (i=0; i<all_discretized_curve_point_.size(); ++i)
	{
        sz = all_discretized_curve_point_[i]->size();
		for (j=1; j<sz-1; ++j)
		{
			global_surface_mesh_.mesh_point_3d.push_back(all_discretized_curve_point_[i]->at(j).point);
			++np;
            all_discretized_curve_point_[i]->at(j).global_id = np - 1;
		} 
	}

	// Finally numbers the interior points on the faces and computes 
    // their 3D coordinates
	for (i=0; i<all_surface_region_mesh_.size(); ++i)
	{
		pmsh = all_surface_region_mesh_[i];
        for (cv_it=pmsh->vertices_begin(); cv_it!=pmsh->vertices_end(); ++cv_it)
		{
			if (!pmsh->is_boundary(cv_it))
			{
				global_surface_mesh_.mesh_point_3d.push_back(pmsh->point(cv_it.handle()));
				++np;
				pmsh->property(global_id, cv_it.handle()) = np - 1;
			}
		}
	}

	// Gets the connectivity array for the surface triangulation for both stl and fli models
	assert(all_surface_region_mesh_.size()==geometry_->get_surface_region().size());

	if ("stl"==geometry_->get_model_type())
	{
		for (i=0; i<all_surface_region_mesh_.size(); ++i)
		{
			psr = geometry_->get_surface_region()[i]; // todo_clg: ok? private. ok!
			pmsh = all_surface_region_mesh_[i];
			cv_it = pmsh->vertices_begin();
			for (j=0; j<psr->curve_number; ++j)
			{
				cur_idx = abs(psr->curve_loop[j]) - 1;
				if (psr->curve_loop[j]>0)
					for (k=0; k<all_discretized_curve_point_[cur_idx]->size()-1; ++k, ++cv_it)
						pmsh->property(global_id, cv_it.handle()) = all_discretized_curve_point_[cur_idx]->at(k).global_id;
				else
					for (k=all_discretized_curve_point_[cur_idx]->size()-1; k>0; --k, ++cv_it)
						pmsh->property(global_id, cv_it.handle()) = all_discretized_curve_point_[cur_idx]->at(k).global_id;
			}

			for (cf_it=pmsh->faces_begin(); cf_it!=pmsh->faces_end(); ++cf_it)
			{
				for (j=0, cfv_it=pmsh->cfv_iter(cf_it.handle()); j<3 && cfv_it; ++j, ++cfv_it)
					elem[j] = pmsh->property(global_id, cfv_it.handle());

				global_surface_mesh_.mesh_element.push_back(elem);
			}
		}
	}
	else if ("fli"==geometry_->get_model_type()) // todo_clg: the same as above!
	{
        for (i=0; i<all_surface_region_mesh_.size(); ++i)
		{
			psr = geometry_->get_surface_region()[i]; // todo_clg: ok? private.
			pmsh = all_surface_region_mesh_[i];
			cv_it = pmsh->vertices_begin();
			for (j=0; j<psr->curve_number; ++j)
			{
				cur_idx = abs(psr->curve_loop[j]) - 1;
				if (psr->curve_loop[j]>0)
					for (k=0; k<all_discretized_curve_point_[cur_idx]->size()-1; ++k, ++cv_it)
						pmsh->property(global_id, cv_it.handle()) = all_discretized_curve_point_[cur_idx]->at(k).global_id;
				else
					for (k=all_discretized_curve_point_[cur_idx]->size()-1; k>0; --k, ++cv_it)
						pmsh->property(global_id, cv_it.handle()) = all_discretized_curve_point_[cur_idx]->at(k).global_id;
			}

			for (cf_it=pmsh->faces_begin(); cf_it!=pmsh->faces_end(); ++cf_it)
			{
				for (j=0, cfv_it=pmsh->cfv_iter(cf_it.handle()); j<3 && cfv_it; ++j, ++cfv_it)
					elem[j] = pmsh->property(global_id, cfv_it.handle());

				global_surface_mesh_.mesh_element.push_back(elem);
			}
		}
	}
	else
		assert(false);

	return;
}

void SurfaceMesher::write_global_mesh()
{
	int np = global_surface_mesh_.mesh_point_3d.size(), ne = global_surface_mesh_.mesh_element.size();
    int i;
	Vec3d v;
	Vec3i elem;

	// off file
    ofstream off_file("global_mesh.off");
	if (!off_file)
	{
		cerr << "error: unable to open file global_mesh.off!" << endl;
		assert(false);
	}
    
	off_file << "OFF" << endl;
	off_file << left << setw(10) << np << " "<< setw(10) << ne << " " << 0 << endl;
	off_file.precision(8);
	off_file << fixed;
	for (i=0; i<np; ++i)
	{
		v = global_surface_mesh_.mesh_point_3d[i];
		off_file << setw(15) << v[0] << setw(15) << v[1] << setw(15) << v[2] << endl;
	}
	off_file << right;
	for (i=0; i<ne; ++i)
	{
		elem = global_surface_mesh_.mesh_element[i];
		off_file << 3 << setw(10) << elem[0] << setw(10) << elem[1] << setw(10) << elem[2] << endl;
	}
	off_file.close();

	// pls file
	ofstream pls_file("global_mesh.pls");
	if (!pls_file)
	{
		cerr << "error: unable to open file global_mesh.pls!" << endl;
		assert(false);
	}
    
	pls_file << left << setw(10) << ne << " " << setw(10) << np << " " << 0 << endl;
	pls_file.precision(8);
	pls_file << fixed;
	for (i=0; i<np; ++i)
	{
		v = global_surface_mesh_.mesh_point_3d[i];
		pls_file << setw(8) << i+1 << setw(15) << v[0] << setw(15) << v[1] << setw(15) << v[2] << endl;
	}
	pls_file << right;
	for (i=0; i<ne; ++i)
	{
		elem = global_surface_mesh_.mesh_element[i];
		pls_file << setw(8) << i+1 << setw(10) << elem[0]+1 << setw(10) << elem[1]+1 << setw(10) 
			<< elem[2]+1 << setw(6) << 0 << setw(6) << 0 << endl;
	}
	pls_file.close();

	return;
}
