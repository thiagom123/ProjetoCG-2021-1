#pragma once

#include <cmath>
#include <iostream>
#include <random>

using std::sqrt;
using std::fabs;

class Vector3D{
public:
	float x, y, z;

	float Norm();
};
// Type aliases for vec3
using Point = Vector3D;
using vec3 = Vector3D;   // 3D point
using Vertex = Vector3D;
using Vertice = Vector3D;
using Vetor = Vector3D;

Vector3D DefVector(Point ponto1, Point ponto2); // Definir um vetor a partir de dois pontos. Resultado aponta de ponto1 para ponto2
Vector3D KProd(float k, Vector3D c); // Produto de escalar por vetor
Vector3D Sumv(Vector3D v1, Vector3D v2); // Soma de vetores
Vector3D Subv(Vector3D v1, Vector3D v2); // Subtracao de vetores
Vector3D VecProd(Vector3D vetor1, Vector3D vetor2); // Produto vetorial v1 x v2
float ScalarProd(Vector3D vetor1, Vector3D vetor2); // Produto interno  <v1.v2>
Vector3D Normalize(Vector3D vetor); // Normalizar vetor. |vetor| = 1



