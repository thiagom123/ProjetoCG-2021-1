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

vector<Objeto> objetos;
int main() {
    //Carregar os arquivos de entrada
    Cena scene;
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
    std::cout << objetos.at(1).vertices.at(1).z << std::endl;
    return 0;
}