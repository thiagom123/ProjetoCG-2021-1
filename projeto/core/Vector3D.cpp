#include <cmath>
#include <iostream>

#include "Vector3D.h"


float Vector3D::Norm(){
	float a = this->x*this->x;
	float b = this->y*this->y;
	float c = this->z*this->z;

	return sqrt(a + b + c);
	
}


// Soma de vetores
Vector3D Sumv(Vector3D v1, Vector3D v2){
	Vector3D out;
	out.x = v1.x + v2.x;
	out.y = v1.y + v2.y;
	out.z = v1.z + v2.z;
	return out;
}

//Funcoes matematicas importantes

// Produto escalar x vetor
Vector3D KProd(float k, Vector3D c){
	Vector3D out;
	out.x = c.x*k;
	out.y = c.y*k;
	out.z = c.z*k;
	return out;
}

Vector3D Subv(Vector3D v1, Vector3D v2){
	Vector3D out;

	out.x = v1.x - v2.x;
	out.y = v1.y - v2.y;
	out.z = v1.z - v2.z;
	return out;
}

//Definindo um vetor a partir de dois pontos
Vector3D DefVector(Point ponto1, Point ponto2){
	Vector3D retorno;

	retorno.x = ponto2.x - ponto1.x;
	retorno.y = ponto2.y - ponto1.y;
	retorno.z = ponto2.z - ponto1.z;

	return retorno;
};


//Produto Escalar
float ProdEscalar(Vector3D vetor1, Vector3D vetor2){
	float  escalar = vetor1.x*vetor2.x + vetor1.y*vetor2.y + vetor1.z*vetor2.z;
	return escalar;
}

//Produto Vetorial
Vector3D ProdVetorial(Vector3D vetor1, Vector3D vetor2){

	Vector3D resposta;
	resposta.x = (vetor1.y*vetor2.z) + (-1 * vetor1.z*vetor2.y);
	resposta.y = (vetor1.z*vetor2.x) + (-1 * vetor1.x*vetor2.z);
	resposta.z = (vetor1.x*vetor2.y) + (-1 * vetor1.y*vetor2.x);

	return resposta;
}

//Normalizar vetor
Vector3D Normalize(Vector3D vetor){

	float modulo = sqrt(vetor.x*vetor.x + vetor.y*vetor.y + vetor.z*vetor.z);

	Vector3D vetorNormalizado;
	vetorNormalizado.x = vetor.x / modulo;
	vetorNormalizado.y = vetor.y / modulo;
	vetorNormalizado.z = vetor.z / modulo;

	return vetorNormalizado;
}

Vector3D divisao(Vector3D a, float b){
	Vector3D retorno;

	retorno.x = a.x/b;
	retorno.y = a.y/b;
	retorno.z = a.z/b;
	//std::cout << retorno.x <<" "<< retorno.y <<" "<< retorno.z << std::endl;
	return retorno;
}

Vector3D pointToVector(Point A){
	Vector3D retorno;

	retorno.x = A.x;
	retorno.y = A.y;
	retorno.z = A.z;

	return retorno;
}

Point vectorToPoint(Vector3D A){
	Point retorno;

	retorno.x = A.x;
	retorno.y = A.y;
	retorno.z = A.z;

	return retorno;
}