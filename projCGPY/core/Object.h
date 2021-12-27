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
	//BoundingBox boundingBox;
	bool isLight = false;
	//bool isTexture = false;
	float coeficienteRefracao;
	float lp;
	//Vector3D normal;

	void normalVertice(); // Pre-processamento. Calcula e armazena a normal de todos os vertices nas faces

	Vector3D coordBaricentricas(Point p, Face f);

	Vector3D normalPoint(Point p, Face f);

};

bool lerObjeto(const char* path, Object &Object);

using Objeto = Object;