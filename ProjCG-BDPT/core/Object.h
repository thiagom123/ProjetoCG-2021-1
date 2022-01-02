#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <vector>

#include "Color.h"
#include "Face.h"
#include "Vector3D.h"
using namespace std;

class Object{
public:
	char path[100];
	Color color;
	float ka, kd, ks, kt, coeficienteEspecular;
	vector<Vertex> vertexs;
	vector<Face> faces;
	bool isLight = false;
	float coeficienteRefracao;
	float lp;

};

bool lerObjeto(const char* path, Object &Object);

using Objeto = Object;