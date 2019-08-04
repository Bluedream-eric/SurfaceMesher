#ifndef BBOX_H
#define BBOX_H

//== INCLUDES =================================================================
#include <OpenMesh\Core\Geometry\VectorT.hh>

//== BOUNDARY BOX CLASS DEFINITION ============================================

/** \class BBox  BBox.h

    Class for boundary box. 
*/

using namespace OpenMesh;

class BBox
{
public:

	//== constructors & basic functions =======================================

	BBox() : min_(Vec3d(0.0,0.0,0.0)), max_(Vec3d(0.0,0.0,0.0))              {}
	BBox(Vec3d _min, Vec3d _max) : min_(_min), max_(_max)                    {}
	~BBox()                                                                  {}

	Vec3d& min() { return min_; } 
	Vec3d& max() { return max_; }
	void set_min(Vec3d _min) { min_ = _min; }
	void set_max(Vec3d _max) { max_ = _max; }
	void set(Vec3d _min, Vec3d _max) { min_ = _min; max_ = _max; }
	void reset() {min_ = (Vec3d)(0.0,0.0,0.0); max_ = (Vec3d)(0.0,0.0,0.0);}

	Vec3d center() { return (min_ + max_) * 0.5; }
	double diagonal() { return (min_ - max_).norm(); }

private:
	Vec3d min_;
	Vec3d max_;
};

#endif