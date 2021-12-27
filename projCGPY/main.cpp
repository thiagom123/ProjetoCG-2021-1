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
#include "core\Color.h"
#include "core\Face.h"
const float view_dist = 600.0;
vector<Objeto> objetos;
const int mDepth = 5;
struct Intersec
{
	public:
		bool hit;
		float distance;
		Point hit_point;
		Vector3D normal;
	Intersec(bool hit, float distance, Point hit_point, Vector3D normal){
		this->hit=hit;
		this->distance=distance;
		this->hit_point=hit_point;
		this->normal=normal;
	}
};

float rand01(float lo,float hi){	
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(lo,hi);
	return dis(gen);	
}


Intersec intersection(Ray ray, Vertex A, Vertex B, Vertex C){
	Vector3D r = Normalize(ray.direction);
	Vector3D temp = Normal(A,B,C);
	//std::cout << "X: " << temp.x << " Y: " << temp.y << " Z: " << temp.z << endl;
	bool hit = false;
	float distance = 0.0;
	Point hit_point;
	if(ProdEscalar(r,temp) == 0.0){
		return Intersec(hit,distance,hit_point,temp);
	}
	Vector3D ray_dir = Normalize(ray.direction);
	float t = ProdEscalar(temp,DefVector(ray.position,A))/ProdEscalar(ray_dir,temp);

	//float d = ProdEscalar(temp, pointToVector(A));

	//float t = (d - ProdEscalar(temp,pointToVector(ray.position)))/(ProdEscalar(temp,ray.direction));
	hit_point = ray.hitpoint(ray,t);
	//std::cout << "X: " << hit_point.x << " Y: " << hit_point.y << " Z: " << hit_point.z << endl;
	//std::cout << "T: " << t << std::endl;
	if(t < 0.0001){
		bool hit = false;
		float distance = 0.0;
		Point hit_point;
		return Intersec(hit,distance,hit_point,temp);
	}

	Vector3D vectorAB = DefVector(A,B);//(B-A)
    Vector3D ac = DefVector(A,C);//(C-A)
    Vector3D vectorBC = DefVector(B,C);//(C-B)
    Vector3D vectorCA = DefVector(C,A);//(A-C)
    Vector3D cb = DefVector(C,B);//(B-C)
    Vector3D ba = DefVector(B,A);//(A-B)

	Vector3D C0 = Subv(pointToVector(hit_point),pointToVector(A));//(Q-A)
	Vector3D C1 = Subv(pointToVector(hit_point),pointToVector(B));//(Q-B)
	Vector3D C2 = Subv(pointToVector(hit_point),pointToVector(C));//(Q-C)
	//std::cout << ProdEscalar(temp,ProdVetorial(vectorAB,C0)) << endl;
	//std::cout << ProdEscalar(temp,ProdVetorial(vectorBC,C1)) << endl;
	//std::cout << ProdEscalar(temp,ProdVetorial(vectorCA,C2)) << endl;
	if(ProdEscalar(temp,ProdVetorial(vectorAB,C0)) > 0 && ProdEscalar(temp,ProdVetorial(vectorBC,C1)) > 0 && ProdEscalar(temp,ProdVetorial(vectorCA,C2)) > 0){
		//std::cout << "Intersecção" << endl;
		hit = true;
        distance = t;
        hit_point = ray.hitpoint(ray,t);
        return Intersec(hit,distance,hit_point,temp);;
	}
	if(ProdEscalar(temp,ProdVetorial(vectorAB,C0)) < 0 && ProdEscalar(temp,ProdVetorial(vectorBC,C1)) < 0 && ProdEscalar(temp,ProdVetorial(vectorCA,C2)) < 0){
		//std::cout << "Intersecção" << endl;
		hit = true;
        distance = t;
        hit_point = ray.hitpoint(ray,t);
        return Intersec(hit,distance,hit_point,temp);;
	}
	return Intersec(hit,distance,hit_point,temp);
}

Color local_color(Object obj, Vector3D hit_normal, Ray ray, Eye eye, float lp){
	Color color = obj.color;

	Vector3D p1 = Normalize(flip_direction(ray.direction));
	Vector3D p2 = hit_normal;

	if(Length(p1) != 1){
		p1 = Normalize(p1);
	}
	if (Length(hit_normal) != 1){
    	p2 = Normalize(hit_normal);
	}

	float angulo = ProdEscalar(p1,p2);

	if(angulo < 0){
		angulo = angulo * -1;
	}

	Color lv = KProdC((angulo * float(obj.kd) * lp),obj.color);
	color = csum(color,lv);

	Vector3D L = flip_direction(ray.direction);
	p1 = Normalize(Vector3D(L.x * (-1), L.y, L.z));
    p2.x = eye.x;
	p2.y = eye.y;
	p2.z = eye.z;
	
	if (Length(p1) != 1){
        p1 = Normalize(p1);
	}
    if (Length(pointToVector(ray.position)) != 1){
        p2 = Normalize(pointToVector(ray.position));
	}

	angulo = ProdEscalar(p1,p2);

	lv = KProdC((obj.ks * pow(angulo, obj.coeficienteRefracao) * lp),Color(1.0,1.0,1.0));
	color = csum(color,lv);

    return color;
}

Color trace_ray(Ray ray, Scene scene, int depth, float nRefractedInitial, int MaxDepth, Eye eye){
	Color difuso = Cor(0,0,0);
	Color especular = Cor(0,0,0);
	Color transmitido = Cor(0,0,0);
	// result = Cor(0,0,0);
	float lp = 1.0;
	int NShadow_Ray = 9;
	int dist = 100;
	int dist2 = 100;
	int dist3 = 100;
	Vector3D hit_point = Vector3D(0.0, 0.0, 0.0);
    Vector3D normal = Vector3D(0.0, 0.0, 0.0);
	float temLuz = 1.0;
	bool hit = false;
	Object ClosestObj;
    Object CurrentObj;
	vector<Face> CurrentFaces;
	for(int i=0; i < objetos.size();i++){
        CurrentObj = objetos.at(i);
        CurrentFaces = CurrentObj.faces;
        for(int j=0;j<CurrentFaces.size();j++){
            Face f = CurrentFaces.at(j);
            Vertex* P1 = f.v1;
            Vertex* P2 = f.v2;
            Vertex* P3 = f.v3;
			Intersec inter = intersection(ray,*P1,*P2,*P3);
			if(inter.hit && inter.distance < dist){
				dist = inter.distance;
				ClosestObj = CurrentObj;
				hit = inter.hit;
				hit_point = pointToVector(inter.hit_point);
				normal = inter.normal;
			}
        }
    }
	if(hit == false){
		return scene.background;
	}
	
	if(ClosestObj.isLight){
			return KProdC(ClosestObj.lp, ClosestObj.color);
	}
	//Shadow Ray
	Color ColorLocal = Color(0,0,0);
	bool shadow = true;
	float lpShadow = 0.0;
	Color ColorShadow 	= Color(0,0,0);
	//Lança Vários Shadow ray
	for (int k = 0; k < NShadow_Ray; k++){
	
		float lx = rand01(scene.light.object->vertexs.at(0).x, scene.light.object->vertexs.at(2).x);
		float lz = rand01(scene.light.object->vertexs.at(0).z, scene.light.object->vertexs.at(2).z);
		Vector3D luz = Normalize(Subv(Vector3D(lx,scene.light.point.y,lz),hit_point));
		Ray shadow_ray2;
		shadow_ray2.position = vectorToPoint(Vector3D(hit_point.x, hit_point.y, hit_point.z));
		shadow_ray2.direction = luz;
		bool hit2 = false;
		Object ClosestObj2;
		Object CurrentObj2;
		Vector3D Normal2;
		vector<Face> CurrentFaces2;
		for(int i=0; i < objetos.size();i++){
			CurrentObj2 = objetos.at(i);
			CurrentFaces2 = CurrentObj.faces;
			for(int j=0;j<CurrentFaces2.size();j++){
				Face f2 = CurrentFaces2.at(j);
				Vertex* P12 = f2.v1;
				Vertex* P22 = f2.v2;
				Vertex* P32 = f2.v3;
				Intersec inter2 = intersection(shadow_ray2,*P12,*P22,*P32);
				if(inter2.hit && inter2.distance < dist2){
					dist2 = inter2.distance;
					ClosestObj2 = CurrentObj2;
					hit2 = inter2.hit;
					Normal2 = inter2.normal;
					//Armazena o objeto mais próximo
									
				}
			}

		}
		//Se o mais próximo for a luz
		if(ClosestObj2.isLight){
			//Atenção: Pegamos o lp da luz, mas o kd do ponto em que estamos calculando
			float lp2= ClosestObj2.lp;
			float kd = ClosestObj.kd;
			float cossenoAng = ProdEscalar(normal, luz);
			//ColorShadow vai ser a média dos lp2*kd*cossenoAng*scene.light.color.r/g/b
			ColorShadow.r += lp2*kd*cossenoAng*scene.light.color.r/((float)NShadow_Ray);
			ColorShadow.g += lp2*kd*cossenoAng*scene.light.color.g/((float)NShadow_Ray);
			ColorShadow.b += lp2*kd*cossenoAng*scene.light.color.b/((float)NShadow_Ray);
			
		
		}	
		
	}
	Color ColorAmbiente = KProdC( (scene.ambient*ClosestObj.ka), ClosestObj.color);			
	Color ColorDireta = csum(ColorAmbiente, ColorShadow);
	float ktot = ClosestObj.kd + ClosestObj.ks + ClosestObj.kt;
	float r = rand01(0,1)*ktot;
	if(depth < MaxDepth){
		if(r < ClosestObj.kd){
			float x = rand01(0,1);
			float y = rand01(0,1);
			Vector3D dir = random_direction(x,y,normal);
			Ray RaioSecundario;
			//float lx = rand01(-0.9100, 0.9100);
			//float lz = rand01(-23.3240, -26.4880);
			//Vector3D luz = Normalize(Subv(Vector3D(lx,3.8360,lz),dir));
			RaioSecundario.position = vectorToPoint(hit_point);
			RaioSecundario.direction = dir;
			difuso = trace_ray(RaioSecundario,scene, depth+1, ClosestObj.coeficienteRefracao, MaxDepth, eye);
			
		}/*else if(r < ClosestObj.kd + ClosestObj.ks){
			Vector3D L = Normalize(flip_direction(ray.direction));			
			Vector3D N = calcularNormal(ClosestObj.faces.at(0).v1,ClosestObj.faces.at(0).v2,ClosestObj.faces.at(0).v3);
			Vector3D R = KProd(2, (Subv(KProd(ProdEscalar(N,L),N),L))); 

			float lx = rand01(-0.9100, 0.9100);
			float lz = rand01(-23.3240, -26.4880);
			Vector3D luz = Normalize(Subv(Vector3D(lx,3.8360,lz),R));
			shadow_ray.position = vectorToPoint(Vector3D(R.x, R.y, R.z));
			shadow_ray.direction = luz;
			bool hit3 = false;
			Object ClosestObj3;
			Object CurrentObj3;
			vector<Face> CurrentFaces3;
			for (int i = 0; i < objetos.size(); i++)
			{
				CurrentObj3 = objetos.at(i);
				CurrentFaces3 = CurrentObj3.faces;
				for (int j = 0; j < CurrentFaces3.size(); j++)
				{
					Face f2 = CurrentFaces3.at(j);
					Vertex *P12 = f2.v1;
					Vertex *P22 = f2.v2;
					Vertex *P32 = f2.v3;
					Intersec inter3 = intersection(shadow_ray, *P12, *P22, *P32);
					if(ClosestObj.isLight){
								lp = ClosestObj.;
							}
					if (inter3.hit && inter3.distance < dist3)
					{
						dist3 = inter3.distance;
						ClosestObj3 = CurrentObj3;
						hit3 = inter3.hit;
						if (hit3)
						{
							if (ClosestObj3.isLight)
							{
								temLuz = 1;
							}
						}
					}
				}
			}
			if(temLuz == 0){
				especular = Color(0,0,0);
			}else{
				Ray new_ray;
				new_ray.position = vectorToPoint(hit_point);
				new_ray.direction = Normalize(R);
				especular =trace_ray(new_ray,depth +1, ClosestObj.coeficienteRefracao, minDepth, eye);
			}
			*/
		//}else{
			//Transmissão
			/*
			if (ClosestObj.kt > 0){
				Vector3D L = Normalize(ray.direction);
				Vector3D N = calcularNormal(ClosestObj.faces.at(0).v1,ClosestObj.faces.at(0).v2,ClosestObj.faces.at(0).v3);
				if(Length(N) != 1){
					N = Normalize(N);
				}
				float cos1 = ProdEscalar(N, flip_direction(L));
				float div = nRefractedInitial/ ClosestObj.coeficienteRefracao;
				float delta = 1-((pow(div,2))*(1-(pow(cos1,2))));
				if(delta >= 0){
					float cos2 = sqrt(delta);
					if(nRefractedInitial != ClosestObj.coeficienteRefracao){
						nRefractedInitial = ClosestObj.coeficienteRefracao;
					}else{
						nRefractedInitial = 1.0;
					}
					
					if(cos1 > 0){
						transmitido = Sumv(KProd(div,L),KProd((div*cos1)-cos2,N));
					}else{
						transmitido = Sumv(KProd(div,L),KProd((div*cos1)+cos2,N));
					}
					float lx = rand01(-0.9100, 0.9100);
					float lz = rand01(-23.3240, -26.4880);
					Vector3D luz = Normalize(Subv(Vector3D(lx,3.8360,lz),transmitido.toVetor()));
					shadow_ray.position = vectorToPoint(Vector3D(transmitido.toVetor().x, transmitido.toVetor().y, transmitido.toVetor().z));
					shadow_ray.direction = luz;
					bool hit3 = false;
					Object ClosestObj3;
					Object CurrentObj3;
					vector<Face> CurrentFaces3;
					for (int i = 0; i < objetos.size(); i++)
					{
						CurrentObj3 = objetos.at(i);
						CurrentFaces3 = CurrentObj3.faces;
						for (int j = 0; j < CurrentFaces3.size(); j++)
						{
							Face f2 = CurrentFaces3.at(j);
							Vertex *P12 = f2.v1;
							Vertex *P22 = f2.v2;
							Vertex *P32 = f2.v3;
							Intersec inter3 = intersection(shadow_ray, *P12, *P22, *P32);
							if(ClosestObj.isLight){
										lp = ClosestObj.;
									}
							if (inter3.hit && inter3.distance < dist3)
							{
								dist3 = inter3.distance;
								ClosestObj3 = CurrentObj3;
								hit3 = inter3.hit;
								if (hit3)
								{
									if (ClosestObj3.isLight)
									{
										temLuz = 1;
									}
								}
							}
						}
					}
					if(temLuz == 0){
						transmitido = Color(0,0,0);
					}else{
						Ray new_ray;
						new_ray.position = vectorToPoint(hit_point);
						new_ray.direction = Normalize(transmitido.toVetor());
						transmitido =trace_ray(new_ray,depth +1, ClosestObj.coeficienteRefracao, minDepth, eye);
					}
				}
			}
			*/
			
		//}
	}
	Color ColorIndireta = KProdC(ClosestObj.kd, difuso);
	return csum(ColorDireta, ColorIndireta);
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

Color Tonemapping(Color pixel, float tmapping){
    //if 0.9999 < pixel.r > 1.0001 :
        pixel.r = pixel.r / (pixel.r + tmapping);
    //if 0.9999 < pixel.g > 1.0001 :
        pixel.g = pixel.g / (pixel.g + tmapping);
    //if 0.9999 < pixel.b > 1.0001 :
        pixel.b = pixel.b / (pixel.b + tmapping);
		return pixel;
}

void render(Eye eye,Scene scene, Window window, int npaths, int depth, double tonemapping,int maxDepth){
        Ray ray;
        ray.position.x = eye.x;
        ray.position.y = eye.y;
        ray.position.z = eye.z;
        double sampleX, sampleY;
        Color color = Color(0,0,0); 
        Color colorAux;
		//std::cout  <<"npath: " << npaths << std::endl;
		std::cout << "P3\n" << window.sizeX<< ' ' << window.sizeY << "\n255\n";
        for (int j = window.sizeY - 1; j >=0 ; j--) {
            for (int i = 0; i < window.sizeX; i++) {
                color.r = 0;
                color.g = 0;
                color.b = 0;
                for (int s=0; s< npaths; s++){
                    sampleX = (i + rand01(0.0, 1.0)) - (window.sizeX / 2.0);
                    sampleY = (j + rand01(0.0, 1.0)) - (window.sizeY / 2.0);
                    //Tem que implementar essa função aqui
                    ray.direction = get_direction(eye, window, sampleX, sampleY);
					//std::cout << ray.direction.x <<" "<< ray.direction.y << " "<<ray.direction.z << endl;
					//std::cout << "Teste" << std::endl;
                    colorAux = trace_ray(ray, scene, depth, 1.0, maxDepth, eye);
					//std::cout << "TesteDepois" << std::endl;
                    color.r = color.r + colorAux.r;
                    color.g = color.g + colorAux.g;
                    color.b = color.b + colorAux.b;
                }
                color.r = color.r / npaths;
                color.g = color.g / npaths;
                color.b = color.b / npaths;
                color = Tonemapping(color, tonemapping);
            	print_color(color);
                //Função de save_pixel(pixel, x, y)
            }
        }
}


/*
O verdadeiro main*/
int main(){
    Scene scene;
    bool temp = LoadScene("cornell_box\\cornellroom.sdl",scene);
    /*if(temp){
        std::cout <<"Objetos na cena:"<< scene.objects.size() << std::endl;
    }else{
        std::cout << "Não Funcionou" << std::endl;
    }*/
	scene.eye = compute_uvw(scene.eye);
	scene.eye.view_dist = view_dist;
    objetos = scene.objects;
	
	//std::cout << " ambient " << scene.ambient << " background " << scene.background.r << scene.background.g << scene.background.b  << " npath "<< scene.npaths << " tonemap " << scene.tonemapping << "sda" <<scene.seed <<endl;
    for (int i = 0; i < scene.objects.size(); i++)
	{
        std::string objPath = "cornell_box\\";
		//char realPath [100]= "cornel_box\\";
		//strcat(objPath, objetos.at(i).path);
        objPath += objetos.at(i).path;
        //std::cout << objPath << std::endl;
		lerObjeto(objPath.c_str(), objetos.at(i));
		//objetos.at(i).normalVertice();
	}
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
	}*/
	/*Ray r;
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
	Intersec intae = intersection(r,A,B,C);
	if(intae.hit){
		std::cout << "Funciona" << std::endl;
	}else{
		std::cout << "Não Funciona" << std::endl;
	}*/
	int depth = mDepth;
	render(scene.eye,scene,scene.window,10,depth,scene.tonemapping,mDepth);

}