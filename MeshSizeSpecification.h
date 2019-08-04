#ifndef MESHSIZESPECIFICATION_H
#define MESHSIZESPECIFICATION_H

//== INCLUDES =================================================================
#include <OpenMesh\Core\Geometry\VectorT.hh>

#include <vector>
#include <string>

//== MESH SIZE SPECIFICATION CLASS DEFINITION =================================

/** \class MeshSizeSpecification  MeshSizeSpecification.h

    Class for defining mesh element characteristics, i.e. element size and shape.
	Both user and geometric specifications are included.
*/

using namespace OpenMesh;

class MeshSizeSpecification
{
public:

	//== definition of mesh sources and background mesh =======================

	// point source
	typedef struct PointSource
	{
		Vec3d center;				        // center position
		double inner_radius, outer_radius;	// outer & inner radius
		double spacing;                     // space value 
		Vec3d bmin, bmax;		            // box of the impacted region
	} PointSource;

	// line source
	typedef struct LineSource
	{
		PointSource point_source[2];		// two point sources
		Vec3f       bmin, bmax;     	    // box of the impacted region
	} LineSource;

	// triangle source
	typedef struct TriangleSource
	{
		PointSource point_source[3];		// three point sources
		Vec3f     bmin, bmax;	            // box of the impacted region
	} TriangleSource;

	// Node for background, kind of trouble, temporarily not used, may need to be modified
	typedef struct Node
	{
		Vec3d point;			// coords
		double spacing;		// space control
	} Node;

    // element (triangle in 2D & tetrahedral in 3D), here it is only used for background mesh 
	typedef struct Elem
	{
		int node[4];
	} Elem;

	//== constructors & destructors ===========================================

	MeshSizeSpecification()
	{
		constant_spacing_ = 1.0;
		min_spacing_ = 1.0;

		point_source_.clear();
		line_source_.clear();
		triangle_source_.clear();

		bkg_node_.clear();
		bkg_elem_.clear();
	}

	~MeshSizeSpecification()                   {}

	//== member functions =====================================================

	bool read_from_ba3(std::string& _filename);

	// modify the backgroun mesh and mesh sources, mesh gradation control
	// todo_clg

	double get_final_spacing(Vec3d &_pt); // cann't be const

	double get_min_spacing() { return min_spacing_; }

private:
	void smooth_mesh_gradation();
	void calculate_min_spacing();

	double get_spacing_from_mesh_source(Vec3d& _pt); // cann't be const
    double get_spacing_from_point_source(PointSource& _pnt_src, Vec3d& _pt) const;
	double get_spacing_from_line_source(LineSource& _ln_src, Vec3d& _pt) const;
	double get_spacing_from_triangle_source(TriangleSource& _tri_src, Vec3d& _pt) const;

private:
	// constant size, here determined by background mesh
    double constant_spacing_;

	double min_spacing_; // and max_spacing_ to trim the spacing

    // mesh sources, only isotropic
	std::vector<PointSource>    point_source_;
	std::vector<LineSource>     line_source_;
	std::vector<TriangleSource> triangle_source_;

	// background mesh, can be anisotropic
	std::vector<Node> bkg_node_;
	std::vector<Elem> bkg_elem_;
};

#endif  // MESHSIZESPECIFICATION_H defined