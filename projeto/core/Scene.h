#include <vector>
#include <cstring>
#include "Color.h"
#include "Object.h"
#include "Vector3D.h"


class Eye{
public:
	float x, y, z;
};

class Window{
public:
	float x0, y0, x1, y1, sizeX, sizeY;
};

class Light {
public:
	Vector3D point;
	Color color;
	float Ip;
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