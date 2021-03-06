#include "Color.h"
#include <iostream>

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

Vector3D Color::toVetor(){
	Vector3D out;
	out.x = this->r;
	out.y = this->g;
	out.z = this->b;

	return out;
}


Color csum(Color c1, Color c2){
	Color output;
	output.r = c1.r + c2.r;
	output.g = c1.g + c2.g;
	output.b = c1.b + c2.b;
	return output;
}

Color::~Color()
{
	// chamada quando delete é utilizado
}
