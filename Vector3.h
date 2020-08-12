//
//  Vector3.H
//  
//

#ifndef Vector3_h
#define Vector3_h

class Vector3{
public:
    float x,y,z;
    
    Vector3(){}
	Vector3(const Vector3 &vector) {
		x = vector.x; y = vector.y; z = vector.z;
	}
    Vector3(float nx, float ny, float nz) {
		x = nx; y = ny; z = nz;
	}
    
Vector3 &operator =(const Vector3 &vector){
    x = vector.x; y = vector.y; z = vector.z;
    return *this;
}

bool operator ==(const Vector3 &vector){
    return x == vector.x && y == vector.y && z == vector.z;
}

bool operator !=(const Vector3 &vector){
    return x != vector.x || y != vector.y || z != vector.z;
}

void zero() { x = 0.0f;  y =0.0f; z = 0.0f; }

Vector3 operator -() const { return Vector3(-x,-y,-z); }

Vector3 operator +(const Vector3 &vector) const {
    return Vector3( x + vector.x, y + vector.y, z + vector.z);
}

Vector3 operator -(const Vector3 &vector) const{
    return Vector3( x - vector.x, y - vector.y, z - vector.z);
}

Vector3 operator *(float s){
    return Vector3( x*s, y*s, z*s);
}

Vector3 operator /(float s){
    float oneOverS = 1.0f / s;
    return Vector3( x*oneOverS, y*oneOverS, z*oneOverS);
}

Vector3 &operator +=(const Vector3 &vector){
    x += vector.x; y += vector.y; z += vector.z;
    return *this;
}

Vector3 &operator -=(const Vector3 &vector){
    x -= vector.x; y -= vector.y; z -= vector.z;
    return *this;
}

Vector3 &operator *=(float s){
    x *= s; y *= s; z *= s;
    return *this;
}

Vector3 &operator /=(float s){
    float oneOverS = 1.0f / s;
    x *= oneOverS; y *= oneOverS; z *= oneOverS;
    return *this;
    
}

//
//Normalize Vector
//Vnorm = V / ||v||
//
void normalize(){
    float magSqrd = x*x + y*y + z*z;
	if (magSqrd > 0.0f) {
		float oneOverMag = 1.0f / sqrt(magSqrd);
		x *= oneOverMag;
		y *= oneOverMag;
		z *= oneOverMag;
	}	
}

//
//Dot Product
// sum of products of vectors
//

float operator *(const Vector3 &vector){
    return (x * vector.x) + (y * vector.y) + (z * vector.z);
}

};

//
//Magnitude
//

inline float vectorMag(const Vector3  &vector){
    return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

//
//Cross Product
//

inline Vector3 crossProduct(const Vector3 &a, const Vector3 &b){
    return Vector3(
                   (a.y * b.z) - (a.z * b.y),
                   (a.z * b.x) - (a.x * b.z),
                   (a.x * b.y) - (a.y * b.x)
                 );
}

//
//multiply by scalar
inline Vector3 operator *(float s, const Vector3 &vector){
    return Vector3(s * vector.x, s * vector.y, s * vector.z);
}

//
//Distance between two points
//

inline float distance(const Vector3 &a, const Vector3 &b){
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    float dz = a.z - b.z;
    return sqrt(dx * dx + dy * dy + dz * dy);
}

extern const Vector3 kZeroVector;



#endif
