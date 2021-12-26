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
const int image_width = 256;
const int image_height = 256;


struct Intersec{
    public:
	Object objeto;
	Point p;
	Vector3D normal;
	bool hit;
    private:
	Face f;
};

struct PointNorm{
    public:
    Point p;
    Vector3D normal;
};

Color difuso(Color ip, float kd, Vector3D lightDir, Vector3D normal, Color corObjeto){

	Color retorno;
	float prodEscalar = ProdEscalar(lightDir, normal);		
	float aux = kd*prodEscalar;
	retorno.r = ip.r*aux*corObjeto.r; 
	retorno.g = ip.g*aux*corObjeto.g;
	retorno.b = ip.b*aux*corObjeto.b;
	
	return retorno;
}

Color especular(Ray ray, Light luz, Vector3D lightDir,Intersec intersection){

	Vector3D rVetor = Subv(KProd(2 * ProdEscalar(lightDir,intersection.normal), intersection.normal), lightDir);
	rVetor = Normalize(rVetor);
	Vector3D vVetor = KProd(-1, ray.direction);
	vVetor = Normalize(vVetor);
	Color especular;
	float aux = pow(ProdEscalar(rVetor, vVetor), intersection.objeto.coeficienteEspecular);
	aux = luz.Ip*intersection.objeto.ks*aux;
	especular.r = luz.color.r*aux;
	especular.g = luz.color.g*aux;
	especular.b = luz.color.b*aux;
    return especular;
}



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
void print_color(Color PixelColor){
    float r = PixelColor.r;
    float g = PixelColor.g;
    float b = PixelColor.b;
    int ir = static_cast<int>(255.999 * r);
    int ig = static_cast<int>(255.999 * g);
    int ib = static_cast<int>(255.999 * b);
    std::cout << ir << ' ' << ig << ' ' << ib << '\n';
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




PointNorm* barycentricCoord(Ray ray,Vertex A, Vertex B, Vertex C){
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
    Point Q = vectorToPoint(Sumv(pointToVector(ray.position),KProd(t,ray.direction)));    
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
        Vector3D nq = divisao(auxNQ,auxNQ.Norm());
        static PointNorm res;
        res.p = Q;
        res.normal = nq;
        return &res;
    }
    return NULL;
}

Intersec ClosestObject(Ray ray){
    Intersec res;
    Point* ClosestIntersection;
    Object ClosestObj;
    Object CurrentObj;
    vector<Face> CurrentFaces;
    double minDist = INT_MAX;
    for(int i=0; i < objetos.size();i++){
        CurrentObj = objetos.at(i);
        CurrentFaces = CurrentObj.faces;
        for(int j=0;j<CurrentFaces.size();j++){
            Face f = CurrentFaces.at(j);
            Vertex* P1 = f.v1;
            Vertex* P2 = f.v2;
            Vertex* P3 = f.v3;
            PointNorm* Intersection = barycentricCoord(ray, *P1,*P2,*P3);
            double d = DistEuclidiana(Intersection->p, ray.position);
            if (d > 0 && d < minDist){
                minDist = d;
                ClosestObj = CurrentObj;
                res.objeto = CurrentObj;
                res.p.x = Intersection->p.x;
                res.p.y = Intersection->p.y;
                res.p.z = Intersection->p.z;
                res.normal = Intersection->normal;
            }
        }
    }
    if (minDist != INT_MAX){
		res.hit = true;
	}else {
		res.hit = false;
	}
    return res;
}

float clamp(float x){ return x<0 ? 0 : x>1 ? 1 : x; };

float tColor(float x){ return pow(clamp(x), 1 / 2.2); };

Color trace_path(int depth, Ray ray, Scene scene, Light luz, int i, int j, int nSample){
    
    float bias = 1e-4;

    Intersec intersection = ClosestObject(ray);
    if (intersection.hit == false){
        return scene.background;
    }else if (intersection.objeto.isLight){
        return luz.color;
    }


    int triangulo = rand() % 2;
	Face triLuz = objetos.at(0).faces.at(triangulo);
	double alpha = rand() %100;
	double beta = rand() % 100;
	double gama = rand() % 100;
	double sum = alpha + beta + gama;
	alpha = alpha / sum;
	beta = beta / sum;
	gama = gama / sum;

	Vertex* v1 = triLuz.v1;
	Vertex* v2 = triLuz.v2;
	Vertex* v3 = triLuz.v3;

	Point lightRand;

	lightRand.x = alpha*v1->x + beta*v2->x+gama*v3->x;
	lightRand.y = v1->y;
	lightRand.z = alpha*v1->z + beta*v2->z + gama*v3->z;

	Vector3D toLight = Normalize(DefVector(intersection.p, vectorToPoint(scene.light.point)));
	float kd = intersection.objeto.kd, ks = intersection.objeto.ks, kt = intersection.objeto.kt;

    //Rambiente = Ia*kar
	float iA = scene.ambient;
	Vector3D ambiente = KProd(iA*intersection.objeto.ka, intersection.objeto.color.toVetor());

	//Rdifuso = Ip*kd(L.N)r
	Color compDifuso = difuso(luz.color, intersection.objeto.kd, toLight, intersection.normal, intersection.objeto.color);

	//Respecular = Ip*ks*(R.V)^n
    Color compEspecular = especular(ray, luz, toLight, intersection);
	/*Vector3D rVetor = Subv(KProd(2 * ProdEscalar(toLight,intersection.normal), intersection.normal), toLight);
	rVetor = Normalize(rVetor);
	Vector3D vVetor = KProd(-1, ray.direction);
	vVetor = Normalize(vVetor);
	Color especular;
	float aux = pow(ProdEscalar(rVetor, vVetor), intersection.objeto.coeficienteEspecular);
	aux = luz.Ip*intersection.objeto.ks*aux;
	especular.r = luz.color.r*aux;
	especular.g = luz.color.g*aux;
	especular.b = luz.color.b*aux;*/

}



/*
O verdadeiro main*/
int main(){
    Scene scene;
    bool temp = LoadScene("cornell_box\\cornellroom.sdl",scene);
    if(temp){
        std::cout <<"Objetos na cena:"<< scene.objects.size() << std::endl;
    }else{
        std::cout << "Não Funcionou" << std::endl;
    }

    objetos = scene.objects;

    for (int i = 0; i < scene.objects.size(); i++)
	{
        std::string objPath = "cornell_box\\";
		//char realPath [100]= "cornel_box\\";
		//strcat(objPath, objetos.at(i).path);
        objPath += objetos.at(i).path;
        std::cout << objPath << std::endl;
		lerObjeto(objPath.c_str(), objetos.at(i));
		objetos.at(i).normalVertice();
	}

    Color Pixel_Color;
    for (int j = image_height-1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            
            //Pixel_Color = Pegar a Cor do Pixel/Path Tracing
            print_color(Pixel_Color);
        }
    }

}



void testeLoadScene(Scene scene){

    bool temp = LoadScene("cornell_box\\cornellroom.sdl",scene);
    if(temp){
        std::cout <<"Objetos na cena:"<< scene.objects.size() << std::endl;
    }else{
        std::cout << "Não Funcionou" << std::endl;
    }
    std::cout << "Câmera:" << scene.eye.x << " " << scene.eye.y << " " << scene.eye.z;
}

void testeLoadObjects(vector<Object> objetos, Scene scene){

    std::cout << objetos.at(1).path << std::endl;
    for (int i = 0; i < scene.objects.size(); i++)
	{
        std::string objPath = "cornell_box\\";
		//char realPath [100]= "cornel_box\\";
		//strcat(objPath, objetos.at(i).path);
        objPath += objetos.at(i).path;
        std::cout << objPath << std::endl;
		lerObjeto(objPath.c_str(), objetos.at(i));
		objetos.at(i).normalVertice();
	}
    std::cout << objetos.at(1).vertexs.at(1).z << std::endl;
}

void testeIntersecao(){
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
    PointNorm* res = barycentricCoord(r,A,B,C);
    //std::cout << "X: " << temp->x << " Y: " << temp->y << " Z: " << temp->z << std::endl;
    //std::cout << " X: " << res->x << " Y: " << res->y << " Z: " << res->z << std::endl;
}