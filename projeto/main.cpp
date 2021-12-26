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
//Profundidade máxima de cada path
//Quantidade de raios secundários permitidos, por path
const int MAX_DEPTH = 5;
const int image_width = 200;
const int image_height = 200;
	
#define M_PI 3.14159265358979323846

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

/* a = b - c */
#define vector(a,b,c) \
	(a)[0] = (b)[0] - (c)[0];	\
	(a)[1] = (b)[1] - (c)[1];	\
	(a)[2] = (b)[2] - (c)[2];

float innerProduct(float* v1, float* v2){
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

//A x B = (a2b3 - a3b2, a3b1 - a1b3, a1b2 - a2b1); a vector quantity
void crossProduct(float* result, float* v1, float* v2){
	result[0] = v1[1] * v2[2] - v1[2] * v2[1];
	result[1] = v1[2] * v2[0] - v1[0] * v2[2];
	result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
float rayIntersectsTriangle(float *p, float *d,
	float *v0, float *v1, float *v2) {

	float e1[3], e2[3], h[3], s[3], q[3];
	float a, f, u, v, t;
	vector(e1, v1, v0);
	vector(e2, v2, v0);

	crossProduct(h, d, e2);
	a = innerProduct(e1, h);
    std::cout<< " A: " <<a << std::endl;
	if (a > -0.00001 && a < 0.00001)
		return(-1);

	f = 1 / a;
	vector(s, p, v0);
	u = f * (innerProduct(s, h));
    std::cout<< " u: " <<u << std::endl;
	if (u < 0.0 || u > 1.0)
		return(-1);

	crossProduct(q, s, e1);
	v = f * innerProduct(d, q);
    std::cout<< " v: " <<v << std::endl;
	if (v < 0.0 || u + v > 1.0)
		return(-1);

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	t = f * innerProduct(e2, q);
    std::cout<< " t: " <<t << std::endl;
	if (t > 0.00001) // ray intersection
		return(t);
	
	else // this means that there is a line intersection
		// but not a ray intersection
		return (-1);

}

float rand01(){	
    random_device rd;
    return ((float)rd()/rd.max());
	
}

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
    SizeWindowX = window.y0 - window.y1;
    x = ((float)i / window.sizeX)*SizeWindowX + window.x0;
	y = ((float)j / window.sizeY)*SizeWindowY + window.y1;
	z = 0;
    //Os Raios partem da câmera
    ray.position.x = eye.x;
    ray.position.y = eye.y;
    ray.position.z = eye.z;
    //Definido a direção dos raios
    ray.direction.x = x - eye.x;
    ray.direction.y = y - eye.y;
    ray.direction.z = z - eye.z;
    //Normalizar o vetor da direção
	std::cout << "DOrigem: " << ray.position.x <<" "<< ray.position.y <<" "<< ray.position.z << endl;
	std::cout << "DNormal: " << ray.direction.x <<" "<< ray.direction.y <<" "<< ray.direction.z << endl;
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
    //std::cout << ir << ' ' << ig << ' ' << ib << '\n';
} 



PointNorm* barycentricCoord(Ray ray,Vertex A, Vertex B, Vertex C){
    Vector3D ab = DefVector(A,B);//(B-A)
    Vector3D ac = DefVector(A,C);//(C-A)
    Vector3D bc = DefVector(B,C);//(C-B)
    Vector3D ca = DefVector(C,A);//(A-C)
    Vector3D cb = DefVector(C,B);//(B-C)
    Vector3D ba = DefVector(B,A);//(A-B)
	std::cout << "Origem: " << ray.position.x <<" "<< ray.position.y <<" "<< ray.position.z << endl;
	std::cout << "Normal: " << ray.direction.x <<" "<< ray.direction.y <<" "<< ray.direction.z << endl;
    std::cout << A.x <<" "<< A.y <<" "<< A.z << endl;
	std::cout << B.x <<" "<< B.y <<" "<< B.z << endl;
	std::cout << C.x <<" "<< C.y <<" "<< C.z << endl;
     //std::cout << ac.x <<" "<< ac.y <<" "<< ac.z;
    //Calculating normal for support plane
    Vector3D vet = ProdVetorial(ac,ab);
    //std::cout << vet.x <<" "<< vet.y <<" "<< vet.z<< endl;
    Vector3D n = divisao(vet,vet.Norm());
    std::cout << n.x <<" "<< n.y <<" "<< n.z << endl;

    if(ProdEscalar(n,ray.direction) == 0){
        return NULL;
    }
    float d = ProdEscalar(n,pointToVector(A));
    //std::cout << d;
    float t = (d - ProdEscalar(n,pointToVector(ray.position)))/(ProdEscalar(n,ray.direction));
    Point Q = vectorToPoint(Sumv(pointToVector(ray.position),KProd(t,ray.direction)));  
	std::cout << "Inter: "<< Q.x <<" "<< Q.y <<" " <<Q.z << std::endl;
	std::cout << " AQ: " << ProdEscalar(ProdVetorial(DefVector(A,Q),ab),n);
	std::cout << " BQ: " <<ProdEscalar(ProdVetorial(DefVector(B,Q),bc),n);
	std::cout << " CQ: " <<ProdEscalar(ProdVetorial(DefVector(C,Q),ca),n) << endl;
    //Vetor qto = pointToVector(Q);
    if(ProdEscalar(ProdVetorial(DefVector(A,Q),ab),n) >= 0 && ProdEscalar(ProdVetorial(DefVector(B,Q),bc),n) >= 0 && ProdEscalar(ProdVetorial(DefVector(C,Q),ca),n) >= 0){
        std::cout << "Passou" << std::endl;
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
			if(Intersection){
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

bool shadowRay(Ray ray, Scene scene){
	bool retorno = false;
	Intersec intersection = ClosestObject(ray);
	Object closest = intersection.objeto;
	if (closest.isLight) return false;	
	if (intersection.hit){
        retorno = true;
    }
	return retorno;
}

Vector3D calcularRefracao(float n1, float n2, Vector3D i, Vector3D n){
	float cosI = -ProdEscalar(i, n);
	float sen2t = pow(n1 / n2, 2)*(1 - pow(cosI, 2));
	Vector3D t = Sumv(KProd(n1 / n2, i),KProd(((n1 / n2)*cosI - sqrt(1 - sen2t)), n));
	return t;
}


Color trace_path(int depth, Ray ray, Scene scene, Light luz, int i, int j, int nSample){
    if (depth >= MAX_DEPTH) return Color(0,0,0);
    float bias = 1e-4;
    Color output;
    Intersec intersection = ClosestObject(ray);
    //std::cout << "teste" << std::endl;
    if (intersection.hit == false){
        return scene.background;
    }else if (intersection.objeto.isLight){
        return luz.color;
    }

    Vector3D normal = intersection.normal;
    normal = Normalize(normal);

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
	Color compDifuso = difuso(luz.color, intersection.objeto.kd, toLight, normal, intersection.objeto.color);

	//Respecular = Ip*ks*(R.V)^n
    Color compEspecular = especular(ray, luz, toLight, intersection);

    //Shadow Ray
	Ray ray2;
	ray2.direction =  DefVector(intersection.p, lightRand);

	//Walk a little bit in the normal direction in order to avoid self intersection
	Vector3D dist = KProd(bias, normal);
	
	ray2.position.x = intersection.p.x + dist.x;
	ray2.position.y = intersection.p.y + dist.y;
	ray2.position.z = intersection.p.z + dist.z;

	bool sombra = shadowRay(ray2,scene);

	////Definindo o valor da cor local
	Color corLocal;
	if (sombra){
		corLocal.r = 0;
		corLocal.g = 0;
		corLocal.b = 0;
	}else{
		corLocal = csum(csum(compDifuso, Color(ambiente)), compEspecular);
	}
	
	//Ray cast from the recursion
	Ray novoRaio;
	novoRaio.position.x = intersection.p.x;
	novoRaio.position.y = intersection.p.y;
	novoRaio.position.z = intersection.p.z;

	// -------------------------recursion for contribution from other objects---------------------------------
	float ktot = kd + ks + kt;
	float r = rand01()*ktot;
	Vector3D direcao, posicao;
	if (r < kd){
		// raio difuso
		float  r1 = 2 * M_PI * rand01();  // random angle around
		float r2 = rand01();           // random distance from center
		float r2s = sqrt(r2);          // square root of distance from center

		Vector3D w = normal;           // set first axis equal to normal
		Vector3D v1;
		v1.x = 0;
		v1.y = 1;
		v1.z = 0;
		Vector3D v2;
		v2.x = 1;
		v2.y = 0;
		v2.z = 0;
		Vector3D u = fabs(w.x) > 0.1 ? v1 : v2;
		u = (Normalize(ProdVetorial(u, w)));      // second axis
		Vector3D v = ProdVetorial(w, u);          // final axis

		// random direction 
		Vector3D psi;
		psi.x = u.x*cos(r1)*r2s + v.x*sin(r1)*r2s + w.x*sqrt(1 - r2);
		psi.y = u.y*cos(r1)*r2s + v.y*sin(r1)*r2s + w.y*sqrt(1 - r2);
		psi.z = u.z*cos(r1)*r2s + v.z*sin(r1)*r2s + w.z*sqrt(1 - r2);

		psi = Normalize(psi);

		direcao = psi;
	}
	else if (r < kd + ks){
		// raio especular
		// direcao: R=2N(NL) - L
		direcao = Subv(KProd(2 * ProdEscalar(normal, toLight), normal), toLight);

	}
	else {
		// raio transmitido
		// TODO
		// objeto opaco? nenhuma cor transmitida
		// caso contrario... verificar refracao


		float cos = ProdEscalar(ray.direction, normal);
		float n1, n2;
		if (cos > 0){
			//Inside the object
			n1 = intersection.objeto.coeficienteRefracao;
			n2 = 1;
			direcao = calcularRefracao(n1, n2, ray.direction, KProd(-1, normal));
			normal = KProd(-1, normal);
		}
		else {
			n1 = 1;
			n2 = intersection.objeto.coeficienteRefracao;
			direcao = calcularRefracao(n1, n2, ray.direction, normal);
		}

		Vector3D dist2 = KProd(bias, normal);
		novoRaio.position.x = intersection.p.x - dist2.x;
		novoRaio.position.y = intersection.p.y - dist2.y;
		novoRaio.position.z = intersection.p.z - dist2.y;
	}

	// Use the direction and position vector to make a ray

	novoRaio.direction = direcao;
	novoRaio.direction = Normalize(novoRaio.direction);
	float cos_theta = ProdEscalar(direcao, normal);

	Color emitance;
	float auxEmtR = ProdEscalar(normal, toLight)*scene.light.Ip*intersection.objeto.color.r;
	float auxEmtG = ProdEscalar(normal, toLight)*scene.light.Ip*intersection.objeto.color.g;
	float auxEmtB = ProdEscalar(normal, toLight)*scene.light.Ip*intersection.objeto.color.b;
	emitance = Color(auxEmtR, auxEmtG, auxEmtB);
	

	Color BRDF;
	BRDF.r = 2 * emitance.r * cos_theta;
	BRDF.g = 2 * emitance.g * cos_theta;
	BRDF.b = 2 * emitance.b * cos_theta;

	Color recursion = trace_path(depth + 1, novoRaio, scene, scene.light, i, j, nSample);
	
	// -----------------------------------output--------------------------------------
	//*****Testar diferentes pesos*****

	/*output.r = (recursion.r  + corLocal.r)*kd;
	output.g = (recursion.g + corLocal.g)*kd;
	output.b = (recursion.b + corLocal.b)*kd;*/
	/*
	output.r = recursion.r*.5 + corLocal.r*.5;
	output.g = recursion.g*.5 + corLocal.g*.5;
	output.b = recursion.b *.5 + corLocal.b*.5;
	*/

	float factor = max(max(kd, ks), kt);

	output.r = (recursion.r + corLocal.r)*factor;
	output.g = (recursion.g + corLocal.g)*factor;
	output.b = (recursion.b + corLocal.b)*factor; 
    std::cout << " RED: " << output.r << std::endl;
	return output;

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
	Ray ray = Pixel_CameraRay(0, 0, scene.window, scene.eye);
	Pixel_Color = trace_path(0, ray, scene, scene.light, 0, 0, scene.npaths);
    /*Color Pixel_Color;
    for (int j = 0; j < image_height; j++) {
        for (int i = 0; i < image_width; i++) {
            Ray ray = Pixel_CameraRay(i, j, scene.window, scene.eye);

            Pixel_Color = trace_path(0, ray, scene, scene.light, i, j, scene.npaths);
            //cout << "i = "<<i<<"\n";
            print_color(Pixel_Color);
        }
    }*/

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