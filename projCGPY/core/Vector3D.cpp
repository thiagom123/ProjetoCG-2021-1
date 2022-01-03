#include <cmath>
#include <iostream>
#include "Ray.h"
#include "Vector3D.h"
#include "Scene.h"
#include "Object.h"
const int range_from  = 0;
const int range_to    = 1;
std::random_device                  rand_dev;
std::mt19937                        generator(rand_dev());
std::uniform_real_distribution<double>  distr(range_from, range_to);

float Vector3D::Norm(){
	float a = this->x*this->x;
	float b = this->y*this->y;
	float c = this->z*this->z;

	return sqrt(a + b + c);
	
}

Vector3D::Vector3D(){
	this->x = 0;
	this->y = 0;
	this->z = 0;
}


Vector3D::Vector3D(float x, float y, float z){
	this->x = x;
	this->y = y;
	this->z = z;
}

Vector3D sample_direction_hemisphere(double r1, double r2){
	float sinTheta = sqrtf(1 - r1 * r1); 
    float phi = 2 * M_PI * r2; 
    float x = sinTheta * cosf(phi); 
    float z = sinTheta * sinf(phi); 
    return Vector3D(x, r1, z); 
}

Vector3D flip_direction(Vector3D res){
	return KProd(-1,res);
}
Vector3D random_Hemisphere_direction(double u1, double u2, Vector3D N){
	Vector3D Nb, Nt;
	if (std::fabs(N.x) > std::fabs(N.y)) {
        Nt = Vector3D(N.z, 0, -N.x); 
		float aux = (N.x * N.x + N.z * N.z);
		float aux_1 = 1/aux;
		Nt = KProd(aux_1, Nt);
	}
    else {
        Nt = Vector3D(0, -N.z, N.y);
		float aux = (N.y * N.y + N.z * N.z);
		float aux_1 = 1/aux;
		Nt = KProd(aux_1, Nt);
	}
	
    Nb = ProdVetorial(N, Nt); 
	Vector3D sample = sample_direction_hemisphere(u1, u2);
	Vector3D v;
	v.x = sample.x * Nb.x + sample.y * N.x + sample.z * Nt.x;
	v.y = sample.x * Nb.y + sample.y * N.y + sample.z * Nt.y;
	v.z = sample.x * Nb.z + sample.y * N.z + sample.z * Nt.z;
	v = Normalize(v);
	return v;
}
Vector3D sample_direction(double r1, double r2){
	float sinTheta = sqrtf(1 - r1 * r1); 
    float phi = 2 * M_PI * r2; 
    float x = sinTheta * cosf(phi); 
    float z = sinTheta * sinf(phi); 
    return Vector3D(x, r1, z); 
}

Vector3D random_direction(double u1, double u2, Vector3D normal){
	Vector3D p = sample_direction(u1,u2);

	Vector3D w = normal;
	Vector3D v = ProdVetorial(Vector3D(0.00319, 1.0, 0.0078),w);
	v = Normalize(v);
	Vector3D u = ProdVetorial(v,w);
	Vector3D hemi_dir = Sumv(Sumv((KProd(p.x,u)),(KProd(p.y,v))),(KProd(p.z,w)));
	return Normalize(hemi_dir);
}

float Length(Vector3D res){
	return sqrt(res.x * res.x + res.y*res.y + res.z*res.z);
}

Vector3D Normal(Vertex A ,Vertex B, Vertex C){
	Vector3D v = DefVector(A,B);
	Vector3D s = DefVector(A,C);

	Vector3D normal = ProdVetorial(v,s);
	return Normalize(normal);
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