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
const int MAX_PATH_DEPTH = 5;

double DistEuclidiana(Point A, Point B){
	double diffx = A.x - B.x;
	double diffy = A.y - B.y;
	double diffz = A.z - B.z;
	double diff = diffx*diffx + diffy*diffy + diffz*diffz;
	return diff;
}

Ray Pixel_CameraRay(int i, int j, Window window, Eye eye){
    Ray ray;
    //Coordenadas x, y e z do pixel e as dimensões x e y da janela
    double x,y,z,SizeWindowX, SizeWindowY;
    SizeWindowX = window.x1 - window.x0;
    SizeWindowX = window.y1 - window.y0;
    x = ((float)x / window.sizeX)*SizeWindowX + window.x0;
	y = ((float)y / window.sizeY)*SizeWindowY + window.y1;
	z = 0;
    //Os Raios partem da câmera
    ray.position.x = eye.x;
    ray.position.y = eye.y;
    ray.position.z = eye.z;
    //Definido a direção dos raios
    ray.direction.x = x - eye.x;
    ray.direction.y = y - eye.y;
    ray.direction.x = z - eye.z;
    //Normalizar o vetor da direção
    ray.direction = Normalize(ray.direction);
    return ray;
}

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
Object* ClosestObject(Ray ray){
    Point* ClosestIntersection;
    Object* ClosestObj;
    Object CurrentObj;
    vector<Face> CurrentFaces;
    double minDist = INT_MAX;
    for(int i = 0; i< objetos.size(); i++){
        CurrentObj = objetos.at(i);
        CurrentFaces = CurrentObj.faces;
        for (int j=0; j< CurrentFaces.size(); j++){
            Face f = CurrentFaces.at(j);
            Vertex* P1 = f.v1;
            Vertex* P2 = f.v2;
            Vertex* P3 = f.v3;
            Point* Intersection = intersection(ray, *P1, *P2, *P3);
            double d = DistEuclidiana(*Intersection, ray.position);
            if (d>0 && d<minDist){
                minDist = d;
                ClosestObj = &CurrentObj;
            }
        }
    };
    return ClosestObj;


}


void print_color(Color PixelColor){
    float r = PixelColor.r;
    float g = PixelColor.g;
    float b = PixelColor.b;
    int ir = static_cast<int>(255.999 * r);
    int ig = static_cast<int>(255.999 * g);
    int ib = static_cast<int>(255.999 * b);
    std::cout << ir << ' ' << ig << ' ' << ib << '\n';
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