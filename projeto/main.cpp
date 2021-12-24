#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include "core\Cena.h"
#include "core\Objeto.h"
#include "core\Ponto.h"
#include "core\Vetor.h"
#include "core\Quadric.h"
#include "core\Raio.h"

vector<Objeto> objetos;

Ponto* intersection(Raio ray,Vertice A, Vertice B, Vertice C){
    Vetor ab = defVetor(A,B);//(B-A)
    Vetor ac = defVetor(A,C);//(C-A)
    Vetor bc = defVetor(B,C);//(C-B)
    Vetor ca = defVetor(C,A);//(A-C)
     //std::cout << ab.x <<" "<< ab.y <<" "<< ab.z;
     //std::cout << ac.x <<" "<< ac.y <<" "<< ac.z;
    //Calculating normal for support plane
    Vetor vet = vetorial(ab,ac);
    std::cout << vet.x <<" "<< vet.y <<" "<< vet.z<< endl;
    Vetor n = divisao(vet,vet.norma());
    std::cout << n.x <<" "<< n.y <<" "<< n.z << endl;
    if(escalar(n,ray.direcao) == 0){
        return NULL;
    }
    float d = escalar(n,pointToVector(A));
    std::cout << d;
    float t = (d - escalar(n,ray.posicao))/(escalar(n,ray.direcao));
    static Ponto Q = vectorToPoint(vsum(ray.posicao,kprod(t,ray.direcao)));
    //Vetor qto = pointToVector(Q);
    if(escalar(vetorial(ab,defVetor(A,Q)),n) >= 0 && escalar(vetorial(bc,defVetor(B,Q)),n) >= 0 && escalar(vetorial(ca,defVetor(C,Q)),n) >= 0){
        return &Q;
    }
    return NULL;
}

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
    Raio r;
    Vetor p;
    p.x = -3;
    p.y= 2;
    p.z=-3;
    Vetor n;
    n.x = -1;
    n.y=1.5;
    n.z=4;
    r.posicao = p;
    r.direcao = n;
    Vertice A,B,C;
    A.x = -5;
    A.y=-1;
    A.z=1;
    B.x = -2.23;
    B.y=3.94;
    B.z=1;
    C.x = -6.21;
    C.y=4.97;
    C.z=1;
    Ponto* temp = intersection(r,A,B,C);
    std::cout << "X: " << temp->x << "Y: " << temp->y << "Z: " << temp->z << std::endl;
    return 0;
}