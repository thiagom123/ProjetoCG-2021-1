#include "Face.h"

//Calcular normal vertice de um triangulo (feito em relaçao a p1) *nao normalizado*
Vector3D calcularNormal(Vertice* p1, Vertice* p2, Vertice* p3){
	Vector3D vetor1;
	Vector3D vetor2;
	vetor1.e[0] = p3->x() - p1->x();
	vetor1.e[1] = p3->y() - p1->y();
	vetor1.e[2] = p3->z() - p1->z();

	vetor2.e[0] = p2->x() - p1->x();
	vetor2.e[1] = p2->y() - p1->y();
	vetor2.e[2] = p2->z() - p1->z();

	return cross(vetor1, vetor2);
}

void Face::calcularNormaisVertices(){
	this->n1 = calcularNormal(this->v1, this->v3, this->v2);
	this->n2 = calcularNormal(this->v2, this->v1, this->v3);
	this->n3 = calcularNormal(this->v3, this->v2, this->v1);
}