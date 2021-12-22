#include <vector>
#include "Color.h"
#include "Object.h"
#include "Point.h"
#include "Quadric.h"

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
	Point point;
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
	vector<Quadric> quadrics;
	Eye eye;
	Window window;
};



bool LoadScene(const char* path, Scene &scene);