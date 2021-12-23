#pragma once
#include "Vector3D.h"
#include "Point.h"

class Ray{

	public:
		Point position;
		Vector3D direction;

		Point PointAtT(double t) {
			return SomarComVetor(position, t*direction);
		}
		
        
};
