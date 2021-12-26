#include "Ray.h"
#include "Point.h"
#include "Vector3D.h"

Point Ray::hitpoint(Ray ray,double t){
    Vector3D sum = Sumv(pointToVector(ray.position),KProd(t,ray.direction));
    return vectorToPoint(sum);
}