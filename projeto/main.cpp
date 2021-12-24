#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include "core\Scene.h"
#include "core\Object.h"
#include "core\Point.h"
#include "core\Vector3D.h"
#include "core\Ray.h"

vector<Objeto> objetos;

Point* intersection(Ray ray,Vertex A, Vertex B, Vertex C){
    Vector3D ab = DefVector(A,B);//(B-A)
    Vector3D ac = DefVector(A,C);//(C-A)
    Vector3D bc = DefVector(B,C);//(C-B)
    Vector3D ca = DefVector(C,A);//(A-C)
     //std::cout << ab.x <<" "<< ab.y <<" "<< ab.z;
     //std::cout << ac.x <<" "<< ac.y <<" "<< ac.z;
    //Calculating normal for support plane
    Vector3D vet = ProdVetorial(ab,ac);
    //std::cout << vet.x <<" "<< vet.y <<" "<< vet.z<< endl;
    Vector3D n = divisao(vet,vet.Norm());
    if(ProdEscalar(n,ray.direction) == 0){
        return NULL;
    }
    float d = ProdEscalar(n,pointToVector(A));
    std::cout << d;
    float t = (d - ProdEscalar(n,pointToVector(ray.position)))/(ProdEscalar(n,ray.direction));
    static Point Q = vectorToPoint(Sumv(pointToVector(ray.position),KProd(t,ray.direction)));
    //Vetor qto = pointToVector(Q);
    if(ProdEscalar(ProdVetorial(ab,DefVector(A,Q)),n) >= 0 && ProdEscalar(ProdVetorial(bc,DefVector(B,Q)),n) >= 0 && ProdEscalar(ProdVetorial(ca,DefVector(C,Q)),n) >= 0){
        return &Q;
    }
    return NULL;
}

Vector3D* barycentricCoord(Ray ray,Vertex A, Vertex B, Vertex C){
    Vector3D ab = DefVector(A,B);//(B-A)
    Vector3D ac = DefVector(A,C);//(C-A)
    Vector3D bc = DefVector(B,C);//(C-B)
    Vector3D ca = DefVector(C,A);//(A-C)
    Vector3D cb = DefVector(C,B);//(B-C)
    Vector3D ba = DefVector(B,A);//(A-B)
     //std::cout << ab.x <<" "<< ab.y <<" "<< ab.z;
     //std::cout << ac.x <<" "<< ac.y <<" "<< ac.z;
    //Calculating normal for support plane
    Vector3D vet = ProdVetorial(ab,ac);
    std::cout << vet.x <<" "<< vet.y <<" "<< vet.z<< endl;
    Vector3D n = divisao(vet,vet.Norm());
    std::cout << n.x <<" "<< n.y <<" "<< n.z << endl;

    if(ProdEscalar(n,ray.direction) == 0){
        return NULL;
    }
    float d = ProdEscalar(n,pointToVector(A));
    //std::cout << d;
    float t = (d - ProdEscalar(n,pointToVector(ray.position)))/(ProdEscalar(n,ray.direction));
    static Point Q = vectorToPoint(Sumv(pointToVector(ray.position),KProd(t,ray.direction)));    
    //Vetor qto = pointToVector(Q);
    if(ProdEscalar(ProdVetorial(ab,DefVector(A,Q)),n) >= 0 && ProdEscalar(ProdVetorial(bc,DefVector(B,Q)),n) >= 0 && ProdEscalar(ProdVetorial(ca,DefVector(C,Q)),n) >= 0){
        Vector3D bq = DefVector(B,Q);//(Q-B)
        Vector3D cq = DefVector(C,Q);//(Q-C)
        Vector3D aq = DefVector(A,Q);//(Q-A)
        //Coordenada baricentricas
        float alpha = (ProdEscalar(ProdVetorial(bc,bq),n))/(ProdEscalar(ProdVetorial(ab,ac),n));
        float beta = (ProdEscalar(ProdVetorial(ca,cq),n))/(ProdEscalar(ProdVetorial(ab,ac),n));
        float gama = (ProdEscalar(ProdVetorial(ab,aq),n))/(ProdEscalar(ProdVetorial(ab,ac),n));
        std::cout << "alpha: " << alpha << " beta: " << beta << " gama: " << gama << std::endl;
        //Calcular normal do triângulo
        Vector3D na = (ProdVetorial(ab,ac));
        std::cout << na.x <<" "<< na.y <<" "<< na.z<< endl;
        Vector3D nb = (ProdVetorial(bc,ba));
        std::cout << nb.x <<" "<< nb.y <<" "<< nb.z<< endl;
        Vector3D nc = (ProdVetorial(ca,cb));
        std::cout << nc.x <<" "<< nc.y <<" "<< nc.z<< endl;
        Vector3D auxNQ = Sumv(Sumv(KProd(alpha,na),KProd(beta,nb)),KProd(gama,nc));
        static Vector3D nq = divisao(auxNQ,auxNQ.Norm());
        return &nq;
    }
    return NULL;
}



vector<Object> objects;
int main() {
    //Carregar os arquivos de entrada
    /*Cena scene;
    bool temp = lerCena("cornell_box\\cornellroom.sdl",scene);

    if(temp){
        std::cout << scene.objetos.size() << std::endl;
    }else{
        std::cout << "Não Passou" << std::endl;
    }
    objetos = scene.objetos;
    std::cout << objetos.at(1).path << std::endl;
    for (int i = 0; i < scene.objetos.size(); i++)
	{
        std::string objPath = "cornell_box\\";
		//char realPath [100]= "cornel_box\\";
		//strcat(objPath, objetos.at(i).path);
        objPath += objetos.at(i).path;
        std::cout << objPath << std::endl;
		lerObjeto(objPath.c_str(), objetos.at(i));
		objetos.at(i).normalVertice();
	}
    std::cout << objetos.at(1).vertices.at(1).z << std::endl;*/
    Ray r;
    Point p;
    p.x = 2;
    p.y= 2;
    p.z=-3;
    Vector3D n;
    n.x = 1;
    n.y=1;
    n.z=4;
    r.position = p;
    r.direction = n;
    Vertex A,B,C;
    A.x = 1;
    A.y=1;
    A.z=1;
    B.x = 3;
    B.y=4;
    B.z=1;
    C.x = 6;
    C.y=1;
    C.z=1;
    //Point* temp = intersection(r,A,B,C);
    Vector3D* res = barycentricCoord(r,A,B,C);
    //std::cout << "X: " << temp->x << " Y: " << temp->y << " Z: " << temp->z << std::endl;
    std::cout << " X: " << res->x << " Y: " << res->y << " Z: " << res->z << std::endl;
    return 0;
}