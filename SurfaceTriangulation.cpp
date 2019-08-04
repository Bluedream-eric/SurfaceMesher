//== INCLUDES =================================================================
#include "GeometryDefinition.h"
#include "MeshSizeSpecification.h"
#include "SurfaceMesher.h"

using namespace std;

int main(int argc, char** argv)
{
	std::string geometry_name(argv[1]);
	std::string::size_type idx = geometry_name.rfind('.');
	assert(idx!=std::string::npos);
	std::string problem_name(geometry_name.substr(0,idx).c_str()); 
	std::string ext_name(geometry_name.substr(idx+1));
	assert("stl"==ext_name || "STL"==ext_name || "fli"==ext_name 
		|| "FLI"==ext_name); // may need to be extended
	std::string ba3_name = problem_name + ".ba3";

	bool bres = false; // false

	Geometry geometry;
	MeshSizeSpecification mesh_size_spec;
	cout << "SurfaceMesher -> Reading geomery ..." << endl;
	bres = geometry.read(geometry_name);
	assert(bres);
	cout << "SurfaceMesher -> Reading mesh size specification ..." << endl;
	bres = mesh_size_spec.read_from_ba3(ba3_name);
	assert(bres);

    SurfaceMesher surface_mesher(&geometry, &mesh_size_spec);
	cout << "Discretizing all curves ..." << endl;
	surface_mesher.discretize_all_curve();
	cout << "Discretizing all surface regions ..." << endl;
	surface_mesher.discretize_all_surface_region();
	cout << "Forming global mesh ..." << endl;
	surface_mesher.form_global_mesh();
	cout << "Writing global mesh ..." << endl;
	surface_mesher.write_global_mesh();

	return 0;
}