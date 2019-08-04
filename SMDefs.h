#ifndef SMDEFS_H
#define SMDEFS_H

//== INCLUDES =================================================================
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>

// using namespace std;
// using namespace OpenMesh;

struct MyTraits : public OpenMesh::DefaultTraits
{
    typedef OpenMesh::Vec3d Point;

    VertexAttributes  ( OpenMesh::Attributes::Status       );
    EdgeAttributes    ( OpenMesh::Attributes::Status       );
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    FaceAttributes    ( OpenMesh::Attributes::Status       );
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

#endif 