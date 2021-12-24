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
    std::cout << vet.x <<" "<< vet.y <<" "<< vet.z<< endl;
    Vector3D n = divisao(vet,vet.Norm());
    std::cout << n.x <<" "<< n.y <<" "<< n.z << endl;
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
    p.x = -3;
    p.y= 2;
    p.z=-3;
    Vector3D n;
    n.x = -1;
    n.y=1.5;
    n.z=4;
    r.position = p;
    r.direction = n;
    Vertex A,B,C;
    A.x = -5;
    A.y=-1;
    A.z=1;
    B.x = -2.23;
    B.y=3.94;
    B.z=1;
    C.x = -6.21;
    C.y=4.97;
    C.z=1;
    Point* temp = intersection(r,A,B,C);
    std::cout << "X: " << temp->x << "Y: " << temp->y << "Z: " << temp->z << std::endl;
    return 0;
}