#ifndef SURFACEMESHER_H
#define SURFACEMESHER_H

//== INCLUDES =================================================================
#include "GeometryDefinition.h"
#include "MeshSizeSpecification.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "SMDefs.h"   // for MyMesh
#include "MatrixT.h"

//== SURFACE MESH CLASS DEFINITION ============================================

/** \class SurfaceMesher  SurfaceMesher.h

    class for surface mesher
*/ 

#define NULL_ELEM  -1
#define NULL_NEIG  -1

typedef struct MeshPoint3D
{
	Vec3d point;
	double isotropic_size;
	int global_id;
	// Mat3by3d metric_tensor; // not needed if isotropic, cost memory
} MeshPoint3D;

typedef struct MeshPoint2D
{
	Vec2d point;
	Mat2by2d metric_tensor; // not needed for linear triangle surface
	double isotropic_size_3d; // for linear triangle surface, only this size is enough
} MeshPoint2D;

typedef struct Element
{
    int element[3];
	int neighbour[3]; // neighbour[i] is the index of element opposite to node element[i], -1 represents null element!
	// int surface_region_id;
} Element;

typedef struct SurfaceMesh
{
    int surface_region_id; // -1 for global surface mesh
	int boundary_point_number; // 0 represents closed global surface mesh
	// std::vector<MeshPoint2D> mesh_point_2d;
	// std::vector<Element>     mesh_element;
	std::vector<Vec3d> mesh_point_3d; 
	std::vector<Vec3i> mesh_element;  // is it ok with the dynamic growth
} SurfaceMesh;

class SurfaceMesher
{
public:

	//== special data structure  ==============================================

	// auxiliary binary tree for discretizing line curve
    typedef struct MidsplitedEdge
    {
	    Vec3d point[2];
	    double riemannian_length;
	    MidsplitedEdge *left, *right, *parent;
    } MidsplitedEdge;

	typedef struct MidsplitedPoint
	{
		Vec3d point;
		double isotropic_size;
		double riemannian_length;
	}MidsplitedPoint;

	//== constructors & destructors & basic functions =========================

	SurfaceMesher() : geometry_(NULL), mesh_size_spec_(NULL)      {}

	SurfaceMesher(Geometry* _geom, MeshSizeSpecification* _mss) 
		: geometry_(_geom), mesh_size_spec_(_mss)                 {}

	// it's not upto surface mesher to release geometry_ and mesh_size_spec_
	~SurfaceMesher() 
	{
        int i;
		for (i=0; i<(int)all_discretized_curve_point_.size(); ++i)
			delete all_discretized_curve_point_[i]; // std::vector<MeshPoint3D>*, is it ok?
		for (i=0; i<(int)all_surface_region_mesh_.size(); ++i)
			delete all_surface_region_mesh_[i];
	}

	void set_geometry(Geometry* _geom) { geometry_ = _geom; }

	Geometry* get_geometry()           { return geometry_; }

	void set_mesh_size_specification(MeshSizeSpecification* _mss)
	{ mesh_size_spec_ = _mss; }

	MeshSizeSpecification* get_mesh_size_specification()
	{ return mesh_size_spec_; }

	//== member functions =====================================================

	// could have different discretization methods
    void discretize_one_curve(Curve* _cv, double _tol=1.e-6);  // private member?

	void discretize_all_curve();

	void discretize_one_surface_region(SurfaceRegion* _sr, double _tol=1.e-6);  // DT or AFT

	void discretize_all_surface_region();

	void enhance_one_surface_mesh(); // parameters? swap diagonal and smoothing

	void enhance_global_surface_mesh(); // parameters? swap diagonal and smoothing

	void form_global_mesh();

	void write_global_mesh();

private:

    // result in ordered initial parametric domain front points
    void project_discrete_boundary_point(SurfaceRegion *_sr, double _tol=1.e-6);

	void dt_in_parametric_space(SurfaceRegion *_sr, double _tol=1.e-6);

	void form_one_surface_mesh(SurfaceRegion *_sr, double _tol=1.e-6);


	void swap_diagonal_for_one_surface_mesh();

	void smooth_one_surface_mesh();

	void swap_diagonal_for_global_surface_mesh();

	void smooth_global_surface_mesh();

	//== specially private member functions ===================================

	void discretize_line_curve(LineCurve *_lc, double _tol=1.e-6);

	void discretize_ferguson_curve_with_psue(FergusonCurve *_fc, double _tol=1.e-6);

	void discretize_ferguson_curve_with_linearization(FergusonCurve *_fc, double _tol=1.e-6);

	void discretize_nurbs_curve(NurbsCurve *_nc, double _tol=1.e-6);


	//== specially private member functions for discretize_line_curve()========

	// use adaptive Simpson quadrature formula to calculate a 3D line's Riemmanian length
    void riemannian_length_3d_line(Vec3d& _pt1, Vec3d& _pt2, double& _rl) const;

	// midsplit line curve, called by function discretize_line_curve
	void create_binary_tree(MidsplitedEdge *_mse) const; // bool could as return type

	// traverse the binary tree to get midsplited/intermediate point, need to tested
	// called by discretize_line_curve
	void traverse_to_get_midsplited_point(MidsplitedEdge *_root, 
		std::vector<MidsplitedPoint> &_midsplited_point);

	void project_discrete_bounary_point_for_linear_triangle_surface_region(SurfaceRegion *_sr, double _tol=1.e-6);

	// assume the triangle surface and the surface region are the same for simplicity
	// called by dt_in_parametric_space()
	void dt_for_linear_triangle_surface_region(SurfaceRegion *_sr, double _tol=1.e-6);

	//== specially private member functions for discretize_ferguson_curve_with_psue()====

	// sample uniformly the spacing of a ferguson curve
	void sample_ferguson_curve_spacing(FergusonCurve *_fc, std::vector<double>& sam_spa_);

	// calculate accumulated finally discretized arc length for a ferguson curve given the sampled spacing
	void calculate_discretized_arc_length_for_ferguson_curve(const std::vector<double>& _sam_spa, 
		std::vector<double>& _arc_len, const double& _sam_len, double _tol=1.0e-6);

	// generate final discrete point
	void generate_final_discrete_point_for_ferguson_curve(FergusonCurve *_fc, 
		const std::vector<double>& _arc_len, std::vector<MeshPoint3D> *_vec_mp3d, double _tol=1.0e-6);

    void project_discrete_bounary_point_for_ferguson_surface_region(SurfaceRegion *_sr, double _tol=1.e-6);

	void dt_for_ferguson_surface_region(SurfaceRegion *_sr, double _tol=1.e-6);

	//== special private member functions for dt_for_linear_triangle_surface_region()==
	void generate_boundary_mesh_for_linear_triangle_surface_region(SurfaceRegion *_sr, 
		SurfaceMesh *_pm, double _tol);
	void generate_final_parametric_mesh_for_linear_triangle_surface_region(SurfaceRegion *_sr, 
		SurfaceMesh *_sm, double _tol);

	void collapse_2d_short_edges(SurfaceRegion *_sr, MyMesh &_mesh, double _tol);
	void discretize_2d_mesh_inner_edges(SurfaceRegion *_sr, MyMesh &_mesh, double _tol);
	void riemannian_length_2d_line(SurfaceRegion *_sr, Vec3d& _pt1, Vec3d& _pt2, double &_rl);
	void create_binary_tree_2(SurfaceRegion *_sr, MidsplitedEdge *_mse);
	void traverse_to_get_midsplited_point_2(MidsplitedEdge *_mse, 
	    std::vector<MidsplitedPoint> &_midsplited_point);
	void filter_param_mesh_inner_edge_discretized_point(SurfaceRegion *_sr, MyMesh &_mesh, double _tol);
	void insert_filtered_discretized_point(SurfaceRegion *_sr, MyMesh &_mesh, double _tol);
	void insert_one_point(SurfaceRegion *_sr, MeshPoint3D &_mp3d, MyMesh &_mesh);
	void area2(Vec3d _pt0, Vec3d _pt1, Vec3d _pt2, double &_area);
	double angle(Vec3d _pt0, Vec3d _pt1, Vec3d _pt2);
	bool is_flip_legal(MyMesh &_mesh, MyMesh::EdgeHandle &_eh, MyMesh::VertexHandle &_vh);
	bool is_delaunay_broken(SurfaceRegion *_sr, MyMesh &_mesh, MyMesh::VertexHandle &_vh, MyMesh::FaceHandle &_fh);
	void calc_riemannian_circumcenter(SurfaceRegion *_sr, MyMesh &_mesh, MyMesh::FaceHandle &_fh, Vec3d &_cc);
	void optimize_with_diagonal_swapping(SurfaceRegion *_sr, MyMesh &_mesh);
	void calc_element_quality(SurfaceRegion *_sr, MyMesh &_mesh, MyMesh::VertexHandle _vh[3], double &_qua);
	void optimize_with_diagonal_swapping_given_threshold(SurfaceRegion *_sr, MyMesh &_mesh, double _qua_ratio_threshold);
	void map_mesh(SurfaceRegion *_sr, MyMesh &_mesh_src, MyMesh &_mesh_tar); // map parametric mesh to physical mesh
	void smooth_with_optimal_shape(SurfaceRegion *_sr, MyMesh &_mesh);
	void smooth_with_unit_length(SurfaceRegion *_sr, MyMesh &_mesh);
	bool seg_seg_int(Vec3d _p0, Vec3d _p1, Vec3d _p2, Vec3d _p3, double& _coef0, double& _coef1, Vec3d& _pt_int);

	void optimize_final_parametric_mesh_for_linear_triangle_surface_region(SurfaceRegion *_sr, 
		SurfaceMesh *_pm, double _tol);

	void form_final_surface_mesh_for_linear_triangle_surface_region(SurfaceRegion *_sr, 
		SurfaceMesh *_sm, double _tol);

	// == special private member functions for project_discrete_bounary_point_for_ferguson_surface_region() ===
	void project_one_loop_for_ferguson_surface_region(SurfaceRegion *_sr, std::vector<int> &_loop, double _tol);
	void boundary_value(SurfaceRegion *_sr, Vec2d _box[2], Vec2d _up[4], Vec2d _us[4], Vec2i _kp[4], 
		Vec2d _xd[4], Vec2d _bl[4]);
	void localm(SurfaceRegion *_sr, Vec3d _cur_pt, Vec2d _uv0, Vec2d _dir, double a, double b, 
		double _tol, double &_xmin, double &_fmin);
	void fguess(SurfaceRegion *_sr, Vec3d _cur_pt, Vec2d &_uk);
	void locu2(FergusonSurface *_fs, Vec3d _cur_pt, Vec2d &_uk, Vec3d &_eva_pt, double &_dis, int &_flag);
	bool feasib(Vec2d _uk, Vec2d _dr, Vec2d _box[2], double _tol);
	bool newdr(Vec2d _uk, Vec2d _gk, Vec2d &_dr, Vec2d _box[2], double _tol);
	void sort(double _a[], int _sz);
	void orientate_initial_parametric_front(SurfaceRegion *_sr);

	//== special private member functions for dt_for_ferguson_surface_region()==
    void generate_boundary_mesh_for_ferguson_surface_region();
	void generate_final_parametric_mesh_for_ferguson_surface_region(SurfaceRegion *_sr, SurfaceMesh *_sm, double _tol);

	void discretize_2d_mesh_inner_edges_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh, double _tol);
	void riemannian_length_2d_line_for_ferguson_surface_region(SurfaceRegion *_sr, Vec3d& _pt1, Vec3d& _pt2, double &_rl);
	void create_binary_tree_for_ferguson_surface_region(SurfaceRegion *_sr, MidsplitedEdge *_mse); // modified version of create_binary_tree_2(...)
	void filter_param_mesh_inner_edge_discretized_point_for_ferguson_surface_region(SurfaceRegion *_sr, 
		MyMesh &_mesh, double _tol);
	void insert_filtered_discretized_point_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh, double _tol);
	void insert_one_point_for_ferguson_surface_region(SurfaceRegion *_sr, MeshPoint3D &_mp3d, MyMesh &_mesh);
	bool is_delaunay_broken_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh, MyMesh::VertexHandle &_vh, MyMesh::FaceHandle &_fh);
	void calc_riemannian_circumcenter_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh, MyMesh::FaceHandle &_fh, Vec3d &_cc);
	void optimize_with_diagonal_swapping_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh);
	void optimize_with_diagonal_swapping_given_threshold_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh, double _qua_ratio_threshold);
	void calc_element_quality_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh, MyMesh::VertexHandle _vh[3], double &_qua);
	void smooth_with_optimal_shape_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh);
	void map_mesh_for_ferguson_surface_region(SurfaceRegion *_sr, MyMesh &_mesh_src, MyMesh &_mesh_tar);

private:

    //== intermediate important result ========================================

	// size==number of all curve segments, free memory in the destructor
	// cann't use std::vector<MeshPoint3D*>, or it cann't get the number of 
	// discretized curve points, which will be used later.
	std::vector< std::vector<MeshPoint3D>* > all_discretized_curve_point_;

	// size==number of loops in the current surface region, free memory in the function
    // remember to free memory after one surface region mesh is generated
	std::vector< std::vector<MeshPoint2D>* > one_initial_front_param_point_;

	std::vector< std::vector<int>* > loops;

	// size==number of surface regions, all surface patch meshes
	// std::vector<SurfaceMesh*> all_surface_region_mesh_;
	std::vector<MyMesh*> all_surface_region_mesh_;

private:
	// input geometry and mesh size data
    Geometry* geometry_;
	MeshSizeSpecification* mesh_size_spec_;

	// output mesh data
	SurfaceMesh global_surface_mesh_;
};

#endif