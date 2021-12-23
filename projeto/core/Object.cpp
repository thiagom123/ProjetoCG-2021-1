#include "Object.h"
#include "Vector3D.h"


void Objeto::normalVertice(){
	for (int i = 0; i < faces.size(); i++)
	{
		faces.at(i).calcularNormaisVertices();
	}
}

bool lerObjeto(const char* path, Objeto &objeto){
	//O método abaixo foi baseado no cógigo encontrado no tutorial de OpenGL:
	//http://www.opengl-tutorial.org/beginners-tutorials/tutorial-7-model-loading/

	FILE * file = fopen(path, "r");
	if (file == NULL){
		printf("Impossible to open the file ! Are you in the right path ? See Tutorial 1 for details\n");
		getchar();
		return false;
	}
	
	float gx = INT_MIN, lx = INT_MAX, gy = INT_MIN, ly = INT_MAX, gz = INT_MIN, lz = INT_MAX;


	while (1){

		char lineHeader[128];
		// read the first word of the line
		int res = fscanf(file, "%s", lineHeader);

		if (res == EOF)
			break; // EOF = End Of File. Quit the loop.

		// else : parse lineHeader

		if (strcmp(lineHeader, "v") == 0){
			Vertice vertex;
			fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
			objeto.vertexs.push_back(vertex);

			if (vertex.x > gx) gx = vertex.x;
			if (vertex.x < lx) lx = vertex.x;
			if (vertex.y > gy) gy = vertex.y;
			if (vertex.y < ly) ly = vertex.y;
			if (vertex.z > gz) gz = vertex.z;
			if (vertex.z < lz) lz = vertex.z;
		}
		else if (strcmp(lineHeader, "f") == 0){
			std::string vertex1, vertex2, vertex3;
			unsigned int vertexIndex[3];
			int matches = fscanf(file, "%d %d %d\n", &vertexIndex[0], &vertexIndex[1], &vertexIndex[2]);

			Face t;
			t.v1Index = vertexIndex[0];
			t.v2Index = vertexIndex[1];
			t.v3Index = vertexIndex[2];
			objeto.faces.push_back(t);
		}
		else{
			// Probably a comment, eat up the rest of the line
			char stupidBuffer[1000];
			fgets(stupidBuffer, 1000, file);
		}

	}

	for (size_t i = 0; i < objeto.faces.size(); i++)
	{
		Face &currentFace = objeto.faces.at(i);
		currentFace.v1 = &objeto.vertexs.at(currentFace.v1Index - 1);
		currentFace.v2 = &objeto.vertexs.at(currentFace.v2Index - 1);
		currentFace.v3 = &objeto.vertexs.at(currentFace.v3Index - 1);
	}

	//objeto.boundingBox.fill(gx, lx, gy, ly, gz, lz);

	return true;
}