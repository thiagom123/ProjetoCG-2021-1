#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "core\Scene.h"
#include "core\Object.h"





vector<Object> objects;
int main() {
    //Carregar os arquivos de entrada
    Scene scene;
    bool temp = LoadScene("cornell_box\\cornellroom.sdl",scene);

    if(temp){
        std::cout << scene.objects.size() << std::endl;
    }else{
        std::cout << "Não Passou" << std::endl;
    }
    objects = scene.objects;
    std::cout << objects.at(1).path << std::endl;
    for (int i = 0; i < scene.objects.size(); i++)
	{
        std::string objPath = "cornell_box\\";
		//char realPath [100]= "cornel_box\\";
		//strcat(objPath, objetos.at(i).path);
        objPath += objects.at(i).path;
        std::cout << objPath << std::endl;
		lerObjeto(objPath.c_str(), objects.at(i));
		objects.at(i).normalVertice();
	}
    std::cout << objects.at(1).vertexs.at(1).z << std::endl;
    return 0;
}