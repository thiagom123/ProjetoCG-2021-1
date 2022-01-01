#pragma once
#include "Vector3D.h"


class Ray{

	public:
		Point position;
		Vector3D direction;

	Point hitpoint(Ray ray,double t);	
        
};
