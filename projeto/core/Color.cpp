#include "Color.h"

Color::Color(){
	this->r = 0;
	this->g = 0;
	this->b = 0;
}
Color::Color(float r, float g, float b){
	this->r = r;
	this->g = g;
	this->b = b;
}
Color::Color(Vector3D v){
	this->r = v.x;
	this->g = v.y;
	this->b = v.z;
}


Color::~Color()
{
	// chamada quando delete é utilizado
}