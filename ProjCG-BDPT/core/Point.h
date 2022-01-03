#pragma once
//#include "Vector3D.h"

class Point{
public:
	float x, y, z;
	/*Point SomarComVetor(Point P, Vector3D v){
			Point P1;
			P1.x = P.x +v.x;
			P1.y = P.y + v.y;
			P1.z = P.z +v.z;
			return P1;
		}*/
	Point(float x, float y, float z);
	Point();	
};	
/*class Vertex : public Point{

};*/
using Vertice = Point;
using Vertex = Point;