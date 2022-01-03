#include "Scene.h"
#include <stdio.h>
Vector3D get_direction(Eye eye, Point Lower_Left_Corner, double x, double y){
        Vector3D direction = Subv( Sumv(KProd(x,eye.u), KProd(y,eye.v)), KProd(eye.view_dist,eye.w));

		return Normalize(direction);
}

Eye compute_uvw(Eye eye){
        // w
        eye.w.x = eye.x; 
        eye.w.y = eye.y;
        eye.w.z = eye.z;

        eye.w = Normalize(eye.w);
        // u
        Vector3D upVec;
        upVec.x = 0;
        upVec.y= 1;
        upVec.z = 0;
        eye.u = ProdVetorial(upVec, eye.w);
        eye.u = Normalize(eye.u);
        // v
        eye.v = ProdVetorial(eye.w, eye.u);
        eye.v = Normalize(eye.v);
		return eye;
}

bool LoadScene(const char* path, Scene &Scene){

    //O método abaixo foi baseado no cógigo encontrado no tutorial de OpenGL:
	//http://www.opengl-tutorial.org/beginners-tutorials/tutorial-7-model-loading/


    FILE * file = fopen(path, "r");
	if (file == NULL){
        printf("Impossível ler o arquivo");
        return false;
    }
    while(true){
        char lineHeader[128];
        int res = fscanf(file, "%s", lineHeader);

        if( res == EOF){
            break;
        }
        if (strcmp(lineHeader, "eye") == 0){
			fscanf(file, "%f %f %f\n", &Scene.eye.x, &Scene.eye.y, &Scene.eye.z);
		}
        else if (strcmp(lineHeader, "size") == 0){
			fscanf(file, "%f %f\n", &Scene.window.nPixelX, &Scene.window.nPixelY);

		}
		else if (strcmp(lineHeader, "ortho") == 0){
			fscanf(file, "%f %f %f %f\n", &Scene.window.x0, &Scene.window.y0, &Scene.window.x1, &Scene.window.y1);
		}
		else if (strcmp(lineHeader, "background") == 0){
			fscanf(file, "%f %f %f\n", &Scene.background.r, &Scene.background.g, &Scene.background.b);
		}
		else if (strcmp(lineHeader, "ambient") == 0){
			fscanf(file, "%f\n", &Scene.ambient);
		}
		else if (strcmp(lineHeader, "light") == 0){
			Object o;
		
			fscanf(file, "%s %f %f %f %f\n", &o.path, &Scene.light.color.r, &Scene.light.color.g, &Scene.light.color.b, &Scene.light.lp);			
			Scene.light.object = &o;
			o.isLight = true;
			o.color.r = Scene.light.color.r;
			o.color.g = Scene.light.color.g;
			o.color.b = Scene.light.color.b;
			o.lp = Scene.light.lp;
			Scene.objects.push_back(o);
		}
		else if (strcmp(lineHeader, "npaths") == 0){
			fscanf(file, "%i", &Scene.npaths);
		}
		else if (strcmp(lineHeader, "tonemapping") == 0){
			fscanf(file, "%f", &Scene.tonemapping);
		}
		else if (strcmp(lineHeader, "seed") == 0){
			fscanf(file, "%i", &Scene.seed);
		}
		else if (strcmp(lineHeader, "object") == 0){
			Object o;
			fscanf(file, "%s %f %f %f %f %f %f %f %f %f\n", &o.path, &o.color.r, &o.color.g, &o.color.b, &o.ka, &o.kd, &o.ks, &o.kt, &o.coeficienteEspecular, &o.coeficienteRefracao);
			Scene.objects.push_back(o);
		}
		else{
			char Ignorar[1000];
			fgets(Ignorar, 1000, file);
		}

    }
	return true;
}