#pragma once
#include "Vector3D.h"
class Color
{
public:
	float r, g, b;

	Color();
	Color(float r, float g, float b);
	Color(Vector3D v);

};