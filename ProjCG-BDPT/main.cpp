#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <fstream>
#include <iomanip>
#include "core\Scene.h"
#include "core\Object.h"
#include "core\Point.h"
#include "core\Vector3D.h"
#include "core\Ray.h"
#include "core\Color.h"
#include "core\Face.h"
#define PI 3.14159265358979323846 /* pi */

vector<Objeto> objetos;
//Alguns parâmetros
const int mDepth = 5;
const int mBounces = 4;
const bool UsarShadowRay = false;
const bool ShadowRayEmTodos = false;
const int NShadow_Ray = 5;
const int nPaths = 20;
const int CornellBox = 2;
const bool ApplyTonemapping = false;


//Variável para determinar qual tipo de BPDT vai ser utilizado
//0 - Sem BiDirectional, Path tracing normal
//1 - Cada Ponto do LightPath é uma fonte de luz secundária
//2 - Lançamos o light Path e cada ponto envia um raio para a câmera
const int BiDirectionalPT = 2;

ofstream file;
Object camera;
Color colorLightPath[400][400];
Color colorEyePath[400][400];

struct Intersec
{
	//Struct que serve como Retorno da função de Intersecção
public:
	bool hit;
	float distance;
	Point hit_point;
	Vector3D normal;
	Intersec(bool hit, float distance, Point hit_point, Vector3D normal)
	{
		this->hit = hit;
		this->distance = distance;
		this->hit_point = hit_point;
		this->normal = normal;
	}
};
struct LightPathPoint
{
	//Pontos do Light Path (Path que sai da Luz no Bidirectional)
	Point p;
	Color c;
	int objectIndex;
	int faceIndex;
};

vector<LightPathPoint> LightPath;

float rand01(float lo, float hi)
{	
	//Essa função gera um número aleatório de 0 a 1
	std::random_device rd;	// Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(lo, hi);
	return dis(gen);
}

Vector3D calcularRefracao(float n1, float n2, Vetor i, Vetor n)
{
	//Calcula a direção do raio após a refração
	float cosI = -ProdEscalar(i, n);

	float sen2t = pow(n1 / n2, 2) * (1 - pow(cosI, 2));

	Vector3D t = Sumv(KProd(n1 / n2, i), KProd(((n1 / n2) * cosI - sqrt(1 - sen2t)), n));
	return t;
}

Intersec intersection(Ray ray, Vertex A, Vertex B, Vertex C)
{
	//Verificar se o raio intersecta o triângulo e calcula o ponto de intersecção
	Vector3D r = Normalize(ray.direction);
	Vector3D temp = Normal(A, B, C);
	//std::cout << "X: " << temp.x << " Y: " << temp.y << " Z: " << temp.z << endl;
	bool hit = false;
	float distance = 0.0;
	Point hit_point;
	if (ProdEscalar(r, temp) == 0.0)
	{
		return Intersec(hit, distance, hit_point, temp);
	}
	Vector3D ray_dir = Normalize(ray.direction);
	float t = ProdEscalar(temp, DefVector(ray.position, A)) / ProdEscalar(ray_dir, temp);

	//float d = ProdEscalar(temp, pointToVector(A));

	//float t = (d - ProdEscalar(temp,pointToVector(ray.position)))/(ProdEscalar(temp,ray.direction));
	hit_point = ray.hitpoint(ray, t);
	//std::cout << "X: " << hit_point.x << " Y: " << hit_point.y << " Z: " << hit_point.z << endl;
	//std::cout << "T: " << t << std::endl;
	if (t < 0.0001)
	{
		bool hit = false;
		float distance = 0.0;
		Point hit_point;
		return Intersec(hit, distance, hit_point, temp);
	}

	Vector3D vectorAB = DefVector(A, B); //(B-A)
	Vector3D ac = DefVector(A, C);		 //(C-A)
	Vector3D vectorBC = DefVector(B, C); //(C-B)
	Vector3D vectorCA = DefVector(C, A); //(A-C)
	Vector3D cb = DefVector(C, B);		 //(B-C)
	Vector3D ba = DefVector(B, A);		 //(A-B)

	Vector3D C0 = Subv(pointToVector(hit_point), pointToVector(A)); //(Q-A)
	Vector3D C1 = Subv(pointToVector(hit_point), pointToVector(B)); //(Q-B)
	Vector3D C2 = Subv(pointToVector(hit_point), pointToVector(C)); //(Q-C)
	//Vector3D VecNormalEmQ =
	//std::cout << ProdEscalar(temp,ProdVetorial(vectorAB,C0)) << endl;
	//std::cout << ProdEscalar(temp,ProdVetorial(vectorBC,C1)) << endl;
	//std::cout << ProdEscalar(temp,ProdVetorial(vectorCA,C2)) << endl;
	if (ProdEscalar(temp, ProdVetorial(vectorAB, C0)) > 0 && ProdEscalar(temp, ProdVetorial(vectorBC, C1)) > 0 && ProdEscalar(temp, ProdVetorial(vectorCA, C2)) > 0)
	{
		//std::cout << "Intersecção" << endl;
		hit = true;
		distance = t;
		hit_point = ray.hitpoint(ray, t);
		return Intersec(hit, distance, hit_point, temp);
		;
	}
	if (ProdEscalar(temp, ProdVetorial(vectorAB, C0)) < 0 && ProdEscalar(temp, ProdVetorial(vectorBC, C1)) < 0 && ProdEscalar(temp, ProdVetorial(vectorCA, C2)) < 0)
	{
		//std::cout << "Intersecção" << endl;
		hit = true;
		distance = t;
		hit_point = ray.hitpoint(ray, t);
		return Intersec(hit, distance, hit_point, temp);
		;
	}
	return Intersec(hit, distance, hit_point, temp);
}

Color trace_ray(Ray ray, Scene scene, int depth, float nRefractedInitial, int MaxDepth, Eye eye, vector<LightPathPoint> lightPaths)
{
	float bias = 1e-4;
	if (depth > MaxDepth)
		return Color(0, 0, 0);
	float lp = scene.light.lp;
	float dist = 100000000;
	float dist2 = 100000000;
	float dist3 = 100000000;
	Vector3D hit_point = Vector3D(0.0, 0.0, 0.0);
	Vector3D normal = Vector3D(0.0, 0.0, 0.0);
	Color ColorBiDirectional = Color(0, 0, 0);
	Color ColorShadow = Color(0, 0, 0);
	Color ColorLocal = Color(0, 0, 0);
	float temLuz = 1.0;
	bool hit = false;
	Object ClosestObj;
	Object CurrentObj;
	vector<Face> CurrentFaces;
	for (int i = 0; i < objetos.size(); i++)
	{
		CurrentObj = objetos.at(i);
		CurrentFaces = CurrentObj.faces;
		for (int j = 0; j < CurrentFaces.size(); j++)
		{
			Face f = CurrentFaces.at(j);
			Vertex *P1 = f.v1;
			Vertex *P2 = f.v2;
			Vertex *P3 = f.v3;
			Intersec inter = intersection(ray, *P1, *P2, *P3);
			if (inter.hit && inter.distance < dist)
			{
				dist = inter.distance;
				ClosestObj = CurrentObj;
				hit = inter.hit;
				hit_point = pointToVector(inter.hit_point);
				normal = inter.normal;
			}
		}
	}
	if (hit == false)
	{
		return scene.background;
	}

	if (ClosestObj.isLight)
	{
		return KProdC(ClosestObj.lp, ClosestObj.color);
	}

	//Shadow Ray

	bool shadow = true;
	float lpShadow = 0.0;

	//Lança Vários Shadow ray
	float xMin = objetos.at(0).vertexs.at(0).x;
	float xMax = objetos.at(0).vertexs.at(2).x;
	float zMin = objetos.at(0).vertexs.at(2).z;
	float zMax = objetos.at(0).vertexs.at(0).z;
	float ly = objetos.at(0).vertexs.at(0).y;
	bool hit2 = false;
	float lx, lz;

	if ((depth == 0 || ShadowRayEmTodos == true) && UsarShadowRay == true)
	{
		for (int k = 0; k < NShadow_Ray; k++)
		{
			//int kx = k%RaizNShadow_Ray;
			//int kz = floor(k/RaizNShadow_Ray);
			//std::cout<<"k: "<< k <<" kx "<< kx <<" kz "<< kz <<endl;
			//float lx = xMin + (kx+0.5)*(xMax-xMin)/RaizNShadow_Ray;
			//float lz = zMin + (kz+0.5)*(zMax-zMin)/RaizNShadow_Ray;
			lx = rand01(xMin, xMax);
			lz = rand01(zMin, zMax);
			//std::cout<< "x:" << " " << xMin << " " << xMax<<endl;
			//std::cout<< "z:" << " " << zMin << " " << zMax<<endl;
			//std::cout<< k <<" "<< kx <<" "<< kz <<endl;
			//std::cout << "lx:: "<<lx << " ly:"<< ly <<" lz:"<< lz <<endl;
			//lx = 0;
			//ly = 3.8360;
			//lz = -24.906;
			dist2 = 100000000;
			Vector3D luz = Normalize(Subv(Vector3D(lx, ly, lz), hit_point));
			Ray shadow_ray2;
			shadow_ray2.position = vectorToPoint(Sumv(KProd(bias, normal), Vector3D(hit_point.x, hit_point.y, hit_point.z)));
			shadow_ray2.direction = luz;
			//bool hit2 = false;
			Object ClosestObj2;
			ClosestObj2.isLight = false;
			//Obs: Se colocar para inicializar como false, ele vai dar falso para todos os objetos
			Object CurrentObj2;
			Vector3D Normal2;
			vector<Face> CurrentFaces2;
			for (int i = 0; i < objetos.size(); i++)
			{
				CurrentObj2 = objetos.at(i);
				CurrentFaces2 = CurrentObj2.faces;
				for (int j = 0; j < CurrentFaces2.size(); j++)
				{
					Face f2 = CurrentFaces2.at(j);
					Vertex *P12 = f2.v1;
					Vertex *P22 = f2.v2;
					Vertex *P32 = f2.v3;
					Intersec inter2 = intersection(shadow_ray2, *P12, *P22, *P32);
					if (inter2.hit && inter2.distance < dist2)
					{
						dist2 = inter2.distance;
						ClosestObj2 = CurrentObj2;
						hit2 = inter2.hit;
						Normal2 = inter2.normal;
						//Normal2 = KProd(bias,inter2.normal);
						//Armazena o objeto mais próximo
					}
				}
			}
			//std::cout << ClosestObj2.isLight <<endl;
			//Se o mais próximo for a luz
			if (ClosestObj2.isLight)
			{
				//Atenção: Pegamos o lp da luz, mas o kd do ponto em que estamos calculando
				float lp2 = scene.light.lp;
				float kd = ClosestObj.kd;
				float cossenoAng = ProdEscalar(normal, luz);
				if (cossenoAng < 0)
					cossenoAng = (-1) * cossenoAng;
				//ColorShadow vai ser a média dos lp2*kd*cossenoAng*scene.light.color.r/g/b
				ColorShadow.r += lp2 * kd * cossenoAng * scene.light.color.r * ClosestObj.color.r / ((float)NShadow_Ray);
				ColorShadow.g += lp2 * kd * cossenoAng * scene.light.color.g * ClosestObj.color.g / ((float)NShadow_Ray);
				ColorShadow.b += lp2 * kd * cossenoAng * scene.light.color.b * ClosestObj.color.b / ((float)NShadow_Ray);
			}
		}
	}

	//Shadow Ray para o LightPath

	if (BiDirectionalPT == 1)
	{	
		//std::cout<<"TESTE 1"<<endl;
		int NLightPoints = LightPath.size();
		//std::cout<<"NLIGHT "<<NLightPoints<<endl;
		for (int k = 0; k < NLightPoints; k++)
		{
			//std::cout<<"TESTE 2"<<endl;
			dist2 = 100000000;
			Vector3D LightPoint = pointToVector(LightPath.at(k).p);
			Vector3D luz = Normalize(Subv(LightPoint, hit_point));
			Ray shadow_ray2;
			shadow_ray2.position = vectorToPoint(Sumv(KProd(bias, normal), Vector3D(hit_point.x, hit_point.y, hit_point.z)));
			shadow_ray2.direction = luz;
			bool hit2 = false;
			Object ClosestObj2;
			//ClosestObj2.isLight = false;
			//Obs: Se colocar para inicializar como false, ele vai dar falso para todos os objetos
			Object CurrentObj2;
			Vector3D Normal2;
			vector<Face> CurrentFaces2;
			for (int i = 0; i < objetos.size(); i++)
			{
				CurrentObj2 = objetos.at(i);
				CurrentFaces2 = CurrentObj2.faces;
				for (int j = 0; j < CurrentFaces2.size(); j++)
				{
					if (!(LightPath.at(k).objectIndex == i && LightPath.at(k).faceIndex == j))
					{
						Face f2 = CurrentFaces2.at(j);
						Vertex *P12 = f2.v1;
						Vertex *P22 = f2.v2;
						Vertex *P32 = f2.v3;
						Intersec inter2 = intersection(shadow_ray2, *P12, *P22, *P32);
						if (inter2.hit)
						{
							hit2 = true;
						}
					}
				}
			}
			if (!hit2)
			{
				//Atenção: Pegamos o lp da luz, mas o kd do ponto em que estamos calculando
				float lp2 = scene.light.lp;
				float kd = ClosestObj.kd;
				float cossenoAng = ProdEscalar(normal, luz);
				if (cossenoAng < 0){
					cossenoAng = (-1) * cossenoAng;
				}
					

				ColorBiDirectional.r += lp2 * kd * cossenoAng * scene.light.color.r*ClosestObj.color.r / ((float)NLightPoints);
				ColorBiDirectional.g += lp2 * kd * cossenoAng * scene.light.color.g*ClosestObj.color.r / ((float)NLightPoints);
				ColorBiDirectional.b += lp2 * kd * cossenoAng * scene.light.color.b*ClosestObj.color.r / ((float)NLightPoints);
				//std:: cout<< "COLOR BD.r  "<<ColorBiDirectional.r<<endl;
			}
		}
	}

	//std::cout<<"Shadow Color:" << " "<<ColorShadow.r << " " <<ColorShadow.g << " "<<ColorShadow.b <<endl;
	Color ColorAmbiente = KProdC((scene.ambient * ClosestObj.ka), ClosestObj.color);
	if (CornellBox == 2)
		ColorAmbiente = Color(0, 0, 0);
	Color ColorDireta = csum(csum(ColorAmbiente, ColorShadow), ColorBiDirectional);
	//std::cout<<"Color:" << " "<<ColorDireta.r << " " <<ColorDireta.g << " "<<ColorDireta.b <<endl;
	//return ColorDireta;
	float ktot = ClosestObj.kd + ClosestObj.ks + ClosestObj.kt;
	float r = rand01(0, 1) * ktot;
	Color ColorIndireto = Color(0, 0, 0);
	Ray new_ray;
	if (depth < MaxDepth)
	{
		if (r < ClosestObj.kd)
		{
			float x = rand01(0, 1);
			float y = rand01(0, 1);
			Vector3D dir = random_Hemisphere_direction(x, y, normal);
			new_ray.position = vectorToPoint(Sumv(KProd(bias, normal), Vector3D(hit_point.x, hit_point.y, hit_point.z)));
			new_ray.direction = Normalize(dir);
			//Pode tirar do CSUM
			ColorIndireto = trace_ray(new_ray, scene, depth + 1, ClosestObj.coeficienteRefracao, MaxDepth, eye, lightPaths);
			//Tem que multiplicar pela cor do objeto
			float cossenoAng = ProdEscalar(normal, new_ray.direction);
			ColorIndireto.r = ColorIndireto.r * ClosestObj.color.r * ClosestObj.kd*cossenoAng;
			ColorIndireto.b = ColorIndireto.b * ClosestObj.color.b * ClosestObj.kd*cossenoAng;
			ColorIndireto.g = ColorIndireto.g * ClosestObj.color.g * ClosestObj.kd*cossenoAng;
		}
		else if (r < ClosestObj.kd + ClosestObj.ks)
		{
			Vector3D L = Normalize(flip_direction(ray.direction));
			Vector3D R = Subv( KProd(2*ProdEscalar(normal, L), normal), L);
			new_ray.position = vectorToPoint(Sumv(KProd(bias, normal), Vector3D(hit_point.x, hit_point.y, hit_point.z)));
			new_ray.direction = Normalize(R);
			ColorIndireto = trace_ray(new_ray, scene, depth + 1, ClosestObj.coeficienteRefracao, MaxDepth, eye, lightPaths);
			ColorIndireto.r = ColorIndireto.r * ClosestObj.ks;
			ColorIndireto.b = ColorIndireto.b * ClosestObj.ks;
			ColorIndireto.g = ColorIndireto.g * ClosestObj.ks;
		}
		else
		{
			float cos = ProdEscalar(ray.direction, normal);
            float n1, n2;
            Vector3D dir;
            if (cos > 0)
            {
                //Está saindo do objeto
                n1 = ClosestObj.coeficienteRefracao;
                n2 = 1;
                dir = calcularRefracao(n1, n2, ray.direction, KProd(-1, normal));
                normal = KProd(-1, normal);
            }
            else
            {
                //Entrando no Objeto
                n1 = 1;
                n2 = ClosestObj.coeficienteRefracao;
                dir = calcularRefracao(n1, n2, ray.direction, normal);
                
            }

            Vetor dist2 = KProd(bias, normal);
            new_ray.position = vectorToPoint(Subv(hit_point, dist2));
            new_ray.direction = Normalize(dir);
            ColorIndireto = trace_ray(new_ray, scene, depth + 1, ClosestObj.coeficienteRefracao, MaxDepth, eye, lightPaths);
		}
	}
	return csum(ColorDireta, ColorIndireto);
}

void print_color(Color PixelColor)
{
	float r = PixelColor.r;
	float g = PixelColor.g;
	float b = PixelColor.b;
	int ir = static_cast<int>(255.999 * r);
	int ig = static_cast<int>(255.999 * g);
	int ib = static_cast<int>(255.999 * b);
	file << ir << ' ' << ig << ' ' << ib << '\n';
}

Color Tonemapping(Color pixel, float tmapping)
{
	pixel.r = pixel.r / (pixel.r + tmapping);
	pixel.g = pixel.g / (pixel.g + tmapping);
	pixel.b = pixel.b / (pixel.b + tmapping);
	return pixel;
}

void CalcularLightPath(Scene scene, Ray lightRay, int bounces, int maxbounces, Color corAtualRaio)
{
	float bias = 1e-4;
	LightPathPoint LightPoint;
	if (bounces > maxbounces)
		return;
	float lp = scene.light.lp;
	float dist = 10000000;
	Vector3D hit_point = Vector3D(0.0, 0.0, 0.0);
	Vector3D normal = Vector3D(0.0, 0.0, 0.0);
	float temLuz = 1.0;
	bool hit = false;
	Object ClosestObj;
	Object CurrentObj;
	vector<Face> CurrentFaces;
	for (int i = 0; i < objetos.size(); i++)
	{
		CurrentObj = objetos.at(i);
		CurrentFaces = CurrentObj.faces;
		for (int j = 0; j < CurrentFaces.size(); j++)
		{
			Face f = CurrentFaces.at(j);
			Vertex *P1 = f.v1;
			Vertex *P2 = f.v2;
			Vertex *P3 = f.v3;
			Intersec inter = intersection(lightRay, *P1, *P2, *P3);
			if (inter.hit && inter.distance < dist)
			{
				dist = inter.distance;
				ClosestObj = CurrentObj;
				hit = inter.hit;
				hit_point = pointToVector(inter.hit_point);
				normal = inter.normal;
				LightPoint.objectIndex = i;
				LightPoint.faceIndex = j;
			}
		}
	}

	if (hit == false)
	{
		if (BiDirectionalPT == 2)
		{
			for (int j = 0; j < camera.faces.size(); j++)
			{
				Face f = camera.faces.at(j);
				Vertex *P1 = f.v1;
				Vertex *P2 = f.v2;
				Vertex *P3 = f.v3;
				Intersec inter = intersection(lightRay, *P1, *P2, *P3);
				//Calcular "i"
				if (inter.hit)
				{
					//std::cout << "Hit" << std::endl;
					//Calcular "i" e "j" da matriz de pixels e adicionar cores
					float imgHeight = scene.window.y1 - scene.window.y0;
					float imgWidth = scene.window.x1 - scene.window.x0;
					int i = floor((hit_point.x - scene.window.x0) * (scene.window.nPixelX - 1) / imgWidth);
					int j = floor((hit_point.y - scene.window.y0) * (scene.window.nPixelY - 1) / imgHeight);
					colorLightPath[i][j] = csum(colorLightPath[i][j], corAtualRaio);
					//Calculo das contribuições
					//Salver cor na matriz de pixels
				}
			}
		}
		return;
	}

	LightPoint.p = vectorToPoint(hit_point);
	//Possivelmente armazenar características do objeto em que houve intersecção
	//Fizemos apenas o difuso e a transmissão
	float ktot = ClosestObj.kd + ClosestObj.kt;
	float r = rand01(0, 1) * ktot;
	Ray new_ray;

	if (BiDirectionalPT == 2){
		Ray LightCameraRay;
		LightCameraRay.position = LightPoint.p;
		Vector3D olho = Vector3D(scene.eye.x, scene.eye.y, scene.eye.z);
		LightCameraRay.direction = Normalize(Subv(olho, pointToVector(LightCameraRay.position)));
		Object CurrentObj2;
		Vector3D Normal2;
		vector<Face> CurrentFaces2;
		bool hitLight = false;
		for (int i = 0; i < objetos.size(); i++){
			CurrentObj2 = objetos.at(i);
			CurrentFaces2 = CurrentObj2.faces;
			for (int j = 0; j < CurrentFaces2.size(); j++){
				Face f2 = CurrentFaces2.at(j);
				Vertex *P12 = f2.v1;
				Vertex *P22 = f2.v2;
				Vertex *P32 = f2.v3;
				Intersec inter2 = intersection(LightCameraRay, *P12, *P22, *P32);
				if (inter2.hit){
					hitLight = true;
					goto endOfLoop;
				}
			}
		}
		endOfLoop:
		if (!hitLight){
			for (int j = 0; j < camera.faces.size(); j++){
				Face f = camera.faces.at(j);
				Vertex *P1 = f.v1;
				Vertex *P2 = f.v2;
				Vertex *P3 = f.v3;
				Intersec inter = intersection(lightRay, *P1, *P2, *P3);
				//Calcular "i"
				if (inter.hit){
					//std::cout << "Hit" << std::endl;
					//Calcular "i" e "j" da matriz de pixels e adicionar cores
					float imgHeight = scene.window.y1 - scene.window.y0;
					float imgWidth = scene.window.x1 - scene.window.x0;
					int i = floor((hit_point.x - scene.window.x0) * (scene.window.nPixelX - 1) / imgWidth);
					int j = floor((hit_point.y - scene.window.y0) * (scene.window.nPixelY - 1) / imgHeight);
					colorLightPath[i][j] = csum(colorLightPath[i][j], corAtualRaio);
					//Calculo das contribuições
					//Salver cor na matriz de pixels
				}
			}
		}
	}
	if (r < ClosestObj.kd)
	{
		float x = rand01(0, 1);
		float y = rand01(0, 1);
		Vector3D dir = random_Hemisphere_direction(x, y, normal);
		new_ray.position = vectorToPoint(Sumv(KProd(bias, normal), Vector3D(hit_point.x, hit_point.y, hit_point.z)));
		new_ray.direction = Normalize(dir);
		float cossenoAng = ProdEscalar(normal, new_ray.direction);
		corAtualRaio.r = corAtualRaio.r * ClosestObj.color.r * ClosestObj.kd*cossenoAng;
		corAtualRaio.b = corAtualRaio.b * ClosestObj.color.b * ClosestObj.kd*cossenoAng;
		corAtualRaio.g = corAtualRaio.g * ClosestObj.color.g * ClosestObj.kd*cossenoAng;
		LightPoint.c = corAtualRaio;
		LightPath.push_back(LightPoint);
		CalcularLightPath(scene, new_ray, bounces + 1, maxbounces, corAtualRaio);
	}
	else
	{
		//Transmissão
		float cos = ProdEscalar(lightRay.direction, normal);
		float n1, n2;
		Vector3D dir;
		if (cos > 0)
		{
			//Inside the object
			n1 = ClosestObj.coeficienteRefracao;
			n2 = 1;
			dir = calcularRefracao(n1, n2, lightRay.direction, KProd(-1, normal));
			normal = KProd(-1, normal);
		}
		else
		{
			n1 = 1;
			n2 = ClosestObj.coeficienteRefracao;
			dir = calcularRefracao(n1, n2, lightRay.direction, normal);
		}

		Vetor dist2 = KProd(bias, normal);
		new_ray.position = vectorToPoint(Subv(hit_point, dist2));
		new_ray.direction = Normalize(dir);

		CalcularLightPath(scene, new_ray, bounces + 1, maxbounces, corAtualRaio);
	}
}

void render(Scene scene, int npaths, int maxDepth, int maxBounces)
{
	Eye eye = scene.eye;
	Window window = scene.window;
	float tonemapping = scene.tonemapping;
	file.open("out.ppm");
	Ray ray, lightRay;
	ray.position.x = eye.x;
	ray.position.y = eye.y;
	ray.position.z = eye.z;
	Vector3D origin = pointToVector(ray.position);
	float window_height = 2.0;
	float window_width = 2.0;

	double sampleX, sampleY;
	Color color = Color(0, 0, 0);
	Color colorAux;
	Vector3D Lower_Left;
	Lower_Left.x = scene.window.x0;
	Lower_Left.y = scene.window.y0;
	Lower_Left.z = 0;
	float xMin = objetos.at(0).vertexs.at(0).x;
	float xMax = objetos.at(0).vertexs.at(2).x;
	float zMin = objetos.at(0).vertexs.at(2).z;
	float zMax = objetos.at(0).vertexs.at(0).z;
	lightRay.position.y = objetos.at(0).vertexs.at(0).y;
	float u1, u2;
	//std::cout  <<"npath: " << npaths << std::endl;
	file << "P3\n"
		 << window.nPixelX << ' ' << window.nPixelY << "\n255\n";
	for (int j = window.nPixelY - 1; j >= 0; j--)
	{
		for (int i = 0; i < window.nPixelX; i++)
		{
			color.r = 0;
			color.g = 0;
			color.b = 0;
			for (int s = 0; s < npaths; s++)
			{
				sampleX = (i + rand01(0.0, 1.0)) * window_width / (window.nPixelX - 1) + Lower_Left.x;
				sampleY = (j + rand01(0.0, 1.0)) * window_height / (window.nPixelY - 1) + Lower_Left.y;
				//Tem que implementar essa função aqui

				Vector3D direction = Subv(Sumv(KProd(sampleX, Vector3D(1, 0, 0)), KProd(sampleY, Vector3D(0, 1, 0))), origin);
				ray.direction = Normalize(direction);

				//std::cout << ray.direction.x <<" "<< ray.direction.y << " "<<ray.direction.z << endl;
				//std::cout << "Teste" << std::endl;
				if (BiDirectionalPT > 0)
				{
					lightRay.position.x = rand01(xMin, xMax);
					lightRay.position.z = rand01(zMin, zMax);
					u1 = rand01(0, 1);
					u2 = rand01(0, 1);
					//Calcular normal da fonte de luz, deve ser a mesmo do ponto
					LightPath.clear();
					Vector3D dirLuz = Normal(*scene.light.object->faces.at(0).v1, *scene.light.object->faces.at(0).v2, *scene.light.object->faces.at(0).v3);
					lightRay.direction = random_Hemisphere_direction(u1, u2, dirLuz);
					CalcularLightPath(scene, lightRay, 0, maxBounces, scene.light.color);
				}
				colorAux = trace_ray(ray, scene, 0, 1.0, maxDepth, eye, LightPath);
				colorEyePath[i][j] = csum(colorEyePath[i][j], colorAux);
				//std::cout << "TesteDepois" << std::endl;
			}
			colorEyePath[i][j].r = colorEyePath[i][j].r / npaths;
			colorEyePath[i][j].g = colorEyePath[i][j].g / npaths;
			colorEyePath[i][j].b = colorEyePath[i][j].b / npaths;
			if (BiDirectionalPT == 2)
			{
				colorLightPath[i][j].r = colorLightPath[i][j].r / npaths;
				colorLightPath[i][j].g = colorLightPath[i][j].g / npaths;
				colorLightPath[i][j].b = colorLightPath[i][j].b / npaths;
			}
		}
		std::cout << "[";
		int pos = window.nPixelX * ((window.nPixelX - j) / window.nPixelX);
		for (int i = 0; i < window.nPixelX; ++i)
		{
			if (i < pos)
				std::cout << "=";
			else if (i == pos)
				std::cout << ">";
			else
				std::cout << " ";
		}
		std::cout << "] " << int(((window.nPixelX - j) / window.nPixelX) * 100.0) << " %\r";
		if (j > 0)
		{
			std::cout.flush();
		}
		else
		{
			std::cout << std::endl;
		}
	}
	Color colorPixel;
	for (int j = window.nPixelY - 1; j >= 0; j--)
	{
		for (int i = 0; i < window.nPixelX; i++)
		{
			colorPixel = csum(colorLightPath[i][j], colorEyePath[i][j]);
			if (ApplyTonemapping == true)
				print_color(Tonemapping(colorPixel, tonemapping));
			else
				print_color(colorPixel);
		}
	}
}

/*
O verdadeiro main*/
int main()
{
	Scene scene;
	if (CornellBox == 1)
		bool temp = LoadScene("cornell_box\\cornellroom.sdl", scene);
	else if (CornellBox == 2)
		bool temp = LoadScene("cornell_box_2\\cornellroom.sdl", scene);

	for (int i = 0; i < scene.window.nPixelX; i++)
	{
		for (int j = 0; j < scene.window.nPixelY; j++)
		{
			colorLightPath[i][j] = Color(0, 0, 0);
			colorEyePath[i][j] = Color(0, 0, 0);
		}
	}

	/*if(temp){
        std::cout <<"Objetos na cena:"<< scene.objects.size() << std::endl;
    }else{
        std::cout << "Não Funcionou" << std::endl;
    }*/
	//Acho que está sem usar
	//scene.eye = compute_uvw(scene.eye);
	objetos = scene.objects;
	std::cout << "Numero de Obj " << objetos.size();
	if (BiDirectionalPT == 2)
	{
		Face f1;
		Vertex f1v1 = Point(scene.window.x0, scene.window.y1, 0);
		Vertex f1v2 = Point(scene.window.x0, scene.window.y0, 0);
		Vertex f1v3 = Point(scene.window.x1, scene.window.y0, 0);
		Vertex f1v4 = Point(scene.window.x1, scene.window.y1, 0);
		f1.v1 = &f1v1;
		f1.v2 = &f1v2;
		f1.v3 = &f1v3;
		Face f2;
		f2.v1 = &f1v1;
		f2.v2 = &f1v3;
		f2.v3 = &f1v4;
		camera.faces.push_back(f1);
		camera.faces.push_back(f2);
	}

	std::cout << " ambient " << scene.ambient << " background " << scene.background.r << scene.background.g << scene.background.b << " npath " << scene.npaths << " tonemap " << scene.tonemapping << " seed " << scene.seed << endl;
	for (int i = 0; i < scene.objects.size(); i++)
	{
		std::string objPath;
		if (CornellBox == 1)
			objPath = "cornell_box\\";
		if (CornellBox == 2)
			objPath = "cornell_box_2\\";
		//char realPath [100]= "cornel_box\\";
		//strcat(objPath, objetos.at(i).path);
		objPath += objetos.at(i).path;
		//std::cout << objPath << std::endl;
		lerObjeto(objPath.c_str(), objetos.at(i));
		std::cout << "Objeto: " << i << endl;
		//objetos.at(i).normalVertice();
	}
	scene.light.object = &objetos.at(0);
	/*
	for(int i = 0; i < objetos.size();i++){
		std::cout << "Objeto: " << i << endl;
		std::cout << "Cor do objeto: " << objetos.at(i).color.r <<" "<< objetos.at(i).color.g <<" "<< objetos.at(i).color.b << endl;
		for (int j = 0; j < objetos.at(i).faces.size(); j++){
			std::cout << "Face: " << j << endl;
			std::cout << objetos.at(i).faces.at(j).v1->x << " "<<objetos.at(i).faces.at(j).v1->y <<" "<< objetos.at(i).faces.at(j).v1->z << endl;
			std::cout << objetos.at(i).faces.at(j).v2->x << " "<<objetos.at(i).faces.at(j).v2->y <<" "<< objetos.at(i).faces.at(j).v2->z << endl;
			std::cout << objetos.at(i).faces.at(j).v3->x << " "<<objetos.at(i).faces.at(j).v3->y <<" "<< objetos.at(i).faces.at(j).v3->z << endl;
		}
	}
	*/
	/*
	Vector3D normal = Vector3D(-1,0,0);
	for (int a=0; a< 100; a++){
		float x = rand01(0,1);
		float y = rand01(0,1);
		Vector3D dir = random_direction(x,y,normal);
		std::cout<< dir.x << " " << dir.y << " " << dir.z<<endl;
	}*/
	render(scene, nPaths, mDepth, mBounces);
	file.close();
	std::cout << "Finalizado" << std::endl;
}