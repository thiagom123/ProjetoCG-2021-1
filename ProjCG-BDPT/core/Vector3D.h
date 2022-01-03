#pragma once
#include <cmath>
#include <iostream>
#include <random>
#include "Point.h"
using std::sqrt;
using std::fabs;

class Vector3D{
public:
	float x, y, z;

	float Norm();
	Vector3D();
	Vector3D(float x, float y, float z);
};

using vec3 = Vector3D;
using Vetor = Vector3D;

Vector3D DefVector(Point ponto1, Point ponto2); // Definir um vetor a partir de dois pontos. Resultado aponta de ponto1 para ponto2
Vector3D divisao(Vector3D a, float b);
Vector3D pointToVector(Point A);
Point vectorToPoint(Vector3D A);
float Length(Vector3D res);
Vector3D random_direction(double u1, double u2, Vector3D normal);
Vector3D flip_direction(Vector3D res);
Vector3D Normal(Vertex A ,Vertex B, Vertex C);
Vector3D KProd(float k, Vector3D c); // Produto de escalar por vetor
Vector3D Sumv(Vector3D v1, Vector3D v2); // Soma de vetores
Vector3D Subv(Vector3D v1, Vector3D v2); // Subtracao de vetores
Vector3D ProdVetorial(Vector3D vetor1, Vector3D vetor2); // Produto vetorial v1 x v2
float ProdEscalar(Vector3D vetor1, Vector3D vetor2); // Produto interno  <v1.v2>
Vector3D Normalize(Vector3D vetor); // Normalizar vetor. |vetor| = 1
Vector3D random_Hemisphere_direction(double u1, double u2, Vector3D N);
Vector3D sample_direction_hemisphere(double r1, double r2);

