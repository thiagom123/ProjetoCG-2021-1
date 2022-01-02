#pragma once
#include <vector>
#include <cstring>
#include "Color.h"
#include "Object.h"
#include "Vector3D.h"


class Eye{
public:
	float x, y, z;
	Vector3D u,v,w;
	float view_dist;
};

class Window{
public:
	float x0, y0, x1, y1, nPixelX, nPixelY;
};

class Light {
public:
	Vector3D point;
	Color color;
	float lp;
	Object* object;
};

class Scene{
public:
	Color background;
	float ambient;
	int npaths;
	float tonemapping;
	int seed;

	Light light;
	vector<Object> objects;
	Eye eye;
	Window window;
};



bool LoadScene(const char* path, Scene &scene);
Vector3D get_direction(Eye eye, Point Lower_Left_Corner, double x, double y);
Eye compute_uvw(Eye eye);
