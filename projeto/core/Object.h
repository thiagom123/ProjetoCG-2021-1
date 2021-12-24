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
	bool isTexture = false;
	float coeficienteRefracao;


	void normalVertice(); // Pre-processamento. Calcula e armazena a normal de todos os vertices nas faces

	//Vetor coordBaricentricas(Ponto p, Face f);

	//Vetor normalPonto(Ponto p, Face f);

};

bool lerObjeto(const char* path, Object &Object);

using Objeto = Object;