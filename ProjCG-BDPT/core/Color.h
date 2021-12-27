#pragma once
#include "Vector3D.h"
class Color
{
public:
	float r, g, b;

	Color();
	Color(float r, float g, float b);
	Color(Vector3D v);

	Vector3D toVetor();
	~Color();
};

Color csum(Color c1, Color c2);
Color KProdC(float k, Color c);
using Cor = Color;
