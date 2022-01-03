#include "Object.h"
#include "Vector3D.h"
#include "Ray.h"


bool lerObjeto(const char* path, Objeto &objeto){
	//O método abaixo foi baseado no cógigo encontrado no tutorial de OpenGL:
	//http://www.opengl-tutorial.org/beginners-tutorials/tutorial-7-model-loading/

	FILE * file = fopen(path, "r");
	if (file == NULL){
		printf("Não foi possível ler o arquivo\n");
		getchar();
		return false;
	}
	


	while (1){

		char lineHeader[128];
		// Leia a primeira palavra da linha
		int res = fscanf(file, "%s", lineHeader);

		if (res == EOF)
			break; 

		// Se não for EOF, continue

		if (strcmp(lineHeader, "v") == 0){
			Vertex vertex;
			fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
			objeto.vertexs.push_back(vertex);

		}
		else if (strcmp(lineHeader, "f") == 0){
			std::string vertex1, vertex2, vertex3;
			unsigned int vertexIndex[3];
			//ACHO QUE NÃO ESTÁ SENDO USADO esse matches
			int matches = fscanf(file, "%d %d %d\n", &vertexIndex[0], &vertexIndex[1], &vertexIndex[2]);

			Face t;
			t.v1Index = vertexIndex[0];
			t.v2Index = vertexIndex[1];
			t.v3Index = vertexIndex[2];
			objeto.faces.push_back(t);
		}
		else{
			// Ignorar o resto
			char ignorar[1000];
			fgets(ignorar, 1000, file);
		}

	}

	for (size_t i = 0; i < objeto.faces.size(); i++)
	{
		Face &currentFace = objeto.faces.at(i);
		currentFace.v1 = &objeto.vertexs.at(currentFace.v1Index - 1);
		currentFace.v2 = &objeto.vertexs.at(currentFace.v2Index - 1);
		currentFace.v3 = &objeto.vertexs.at(currentFace.v3Index - 1);
	}


	return true;
}