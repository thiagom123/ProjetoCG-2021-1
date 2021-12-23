#pragma once
#include "Vector3D.h"

class Face{
public:
	Vector3D *v1, *v2, *v3; // Ponteiros para os vertices desse triangulo
	Vector3D n1, n2, n3; // Normais de cada vertice desse triangulo
	int v1Index, v2Index, v3Index;

	void calcularNormaisVertices();
};