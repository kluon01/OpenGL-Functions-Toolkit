//---------------------------------------
// Program: project4.cpp
// Purpose: Implement a toolkit of graphics functions
// Author:  Karshin Luong, John Gauch
// Date:    3/27/2019
//---------------------------------------
#include <math.h>
#include <cstdlib>
#include <iostream>

using namespace std;

int * translate(int x, int y, int z, int Tx, int Ty, int Tz)
{
    int * trans_xyz = new int[3];
    int x_prime, y_prime, z_prime;

    x_prime = x + Tx;
    y_prime = y + Ty;
    z_prime = z + Tz;
    trans_xyz[0] = x_prime;
    trans_xyz[1] = y_prime;
    trans_xyz[2] = z_prime;

    return trans_xyz;
}

double * rotateX(int x, int y, int z, double theta)
{
    double * rotX_xyz = new double[3];
    double x_prime, y_prime, z_prime;

    x_prime = x;
    y_prime = y * cos(theta) - z * sin(theta);
    z_prime = y * sin(theta) + z * cos(theta);
    rotX_xyz[0] = x_prime;
    rotX_xyz[1] = y_prime;
    rotX_xyz[2] = z_prime;

    return rotX_xyz;
}

double * rotateY(int x, int y, int z, double theta)
{
    double * rotY_xyz = new double[3];
    double x_prime, y_prime, z_prime;

    x_prime = x * cos(theta) + z * sin(theta);
    y_prime = y;
    z_prime = -x * sin(theta) + z * cos(theta);
    rotY_xyz[0] = x_prime;
    rotY_xyz[1] = y_prime;
    rotY_xyz[2] = z_prime;

    return rotY_xyz;
}

double * rotateZ(int x, int y, int z, double theta)
{
    double * rotZ_xyz = new double[3];
    double x_prime, y_prime, z_prime;

    x_prime = x * cos(theta) - y * sin(theta);
    y_prime = x * sin(theta) + y * cos(theta);
    z_prime = z;
    rotZ_xyz[0] = x_prime;
    rotZ_xyz[1] = y_prime;
    rotZ_xyz[2] = z_prime;

    return rotZ_xyz;
}

int * scale(int x, int y, int z, int Sx, int Sy, int Sz)
{
    int * scale_xyz = new int[3];
    int x_prime, y_prime, z_prime;

    x_prime = x * Sx;
    y_prime = y * Sy;
    z_prime = z * Sz;
    scale_xyz[0] = x_prime;
    scale_xyz[1] = y_prime;
    scale_xyz[2] = z_prime;

    return scale_xyz;
}

double * project(int x, int y, int z)
{
    double * projectXpYp = new double[2];
    double Xp, Yp;
    double d = 0.5;

    Xp = x / (z/d);
    Yp = y / (z/d);

    projectXpYp[0] = Xp;
    projectXpYp[1] = Yp;

    return projectXpYp;
}

string outcode(double x, double y, double z)
{
    string code = "000000";
    int x_min = 0;
    int x_max = 1;
    int y_min = 0;
    int y_max = 1;
    int z_min = 0;
    int z_max = 1;

	if (x < x_min)           // to the left of clip window
		code[5] = '1';
	else if (x > x_max)      // to the right of clip window
		code[4] = '1';
	if (y < y_min)           // below the clip window
		code[3] = '1';
	else if (y > y_max)      // above the clip window
		code[2] = '1';
    if (z < z_min)
        code[1] = '1';
    else if (z > z_max)
        code[0] = '1';

	return code;
}

bool acceptOutcode(string outcode1, string outcode2)
{
    return (outcode1 == "000000" && outcode2 == "000000");
}

bool rejectOutcode(string outcode1, string outcode2)
{
    if(outcode1[0] == '1' && outcode2[0] == '1')
        return true;
    else if(outcode1[1] == '1' && outcode2[1] == '1')
        return true;
    else if(outcode1[2] == '1' && outcode2[2] == '1')
        return true;
    else if(outcode1[3] == '1' && outcode2[3] == '1')
        return true;
    else if(outcode1[4] == '1' && outcode2[4] == '1')
        return true;
    else if(outcode1[5] == '1' && outcode2[5] == '1')
        return true;
    else
        return false;
}

bool clipOutcode(bool accept, bool reject)
{
    if(accept == 0 && reject == 0)
        return true;
    else
        return false;
}

double * normalize(int Ax, int Ay, int Az)
{
    double * normalPts = new double[3];
    double length = 0;
    double Bx, By, Bz;
    length = sqrt((Ax * Ax) + (Ay * Ay) + (Az * Az));
    Bx = Ax / length;
    By = Ay / length;
    Bz = Az / length;

    normalPts[0] = Bx;
    normalPts[1] = By;
    normalPts[2] = Bz;

    return normalPts;
}

double dotProduct(double Ax, double Ay, double Az, double * normals)
{
    double dot = 0.0;
    dot = (Ax * normals[0]) + (Ay * normals[1]) + (Az * normals[2]);
    return dot;
}

double * crossProduct(double Ax, double Ay, double Az, double * normals)
{ 
    double * crossPts = new double[3];
    double Cx, Cy, Cz;
    Cx = (Ay * normals[2]) - (Az * normals[1]);
    Cy = (Az * normals[0]) - (Ax * normals[2]);
    Cz = (Ax * normals[1]) - (Ay * normals[0]);

    crossPts[0] = Cx;
    crossPts[1] = Cy;
    crossPts[2] = Cz;

    return crossPts;
}

double calcDiffuse(int Lx, int Ly, int Lz, double Nx, double Ny, double Nz)
{
    double dot;
    double diffuseMaterial = 1.0;
    double diffuseLightColor = 1.0;
    double * Npts = new double[3];
    Npts[0] = Nx;
    Npts[1] = Ny;
    Npts[2] = Nz;
    dot = dotProduct(Lx,Ly,Lz,Npts) * diffuseMaterial * diffuseLightColor;
    return dot;
}

double * calcIdealRef(int Lx, int Ly, int Lz, double Nx, double Ny, double Nz)
{
    double dot;
    double * Rpts = new double[3];
    double * Npts = new double[3];
    double Rx, Ry, Rz;
    Npts[0] = Nx;
    Npts[1] = Ny;
    Npts[2] = Nz;
    dot = dotProduct(Lx,Ly,Lz,Npts);
    Rx = (2*dot) * (Npts[0] - Lx); 
    Ry = (2*dot) * (Npts[1] - Ly); 
    Rz = (2*dot) * (Npts[2] - Lz);

    Rpts[0] = Rx;
    Rpts[1] = Ry;
    Rpts[2] = Rz;

    return Rpts; 
}

double calcSpecular(double Vx, double Vy, double Vz, double * Rpts, double P)
{
    double dot = dotProduct(Vx,Vy,Vz,Rpts);
    double specular = pow(dot,P);
    return specular;
}

int main(int argc, char *argv[])
{
    string outcode1 = "100000";
    string outcode2 = "000001";
    double Ax = 3;
    double Ay = 1;
    double Az = 2;
    int Lx = 1;
    int Ly = 2;
    int Lz = 1;
    double Nx = 1.2;
    double Ny = 1.2;
    double Nz = 1.2;
    double Vx = 0.5;
    double Vy = 0.5;
    double Vz = 0.5;
    double P = 2;

    // Geometric Ops
    int * transPtr = translate(1,1,1,2,3,4);
    printf("(1,1,1) translated by 2, 3, 4: x' = %d, y' = %d, z' = %d \n", transPtr[0], transPtr[1], transPtr[2]);
    double * rotX_ptr = rotateX(1,1,1,45.0);
    printf("(1,1,1) rotated on X by 45 degrees: x' = %f, y' = %f, z' = %f \n", rotX_ptr[0], rotX_ptr[1], rotX_ptr[2]);
    double * rotY_ptr = rotateY(1,1,1,45.0);
    printf("(1,1,1) rotated on Y by 45 degrees: x' = %f, y' = %f, z' = %f \n", rotY_ptr[0], rotY_ptr[1], rotY_ptr[2]);
    double * rotZ_ptr = rotateZ(1,1,1,45.0);
    printf("(1,1,1) rotated on Z by 45 degrees: x' = %f, y' = %f, z' = %f \n", rotZ_ptr[0], rotZ_ptr[1], rotZ_ptr[2]);
    int * scalePtr = scale(2,2,2,2,3,4);
    printf("(2,2,2) scaled by 2, 3, 4: x' = %d, y' = %d, z' = %d \n", scalePtr[0], scalePtr[1], scalePtr[2]);
    double * projectPtr = project(1,1,1);
    printf("(1,1,1) projected on z = d plane: Xp = %f, Yp = %f \n", projectPtr[0], projectPtr[1]);

    // Clipping Ops
    string code = outcode(0.8,1.8,-0.8);
    printf("(0.8,1.8,-0.8) outcode: %s \n", code.c_str());
    bool accept = acceptOutcode(outcode1, outcode2);
    printf("Accept line with outcodes %s and %s?(no = 0, yes = 1): %d \n", outcode1.c_str(), outcode2.c_str(), accept);
    bool reject = rejectOutcode(outcode1, outcode2);
    printf("Reject line with outcodes %s and %s?(no = 0, yes = 1): %d \n", outcode1.c_str(), outcode2.c_str(), reject);
    bool clip = clipOutcode(accept, reject);
    printf("Clip line with outcodes %s and %s?(no = 0, yes = 1): %d \n", outcode1.c_str(), outcode2.c_str(), clip);

    // Shading Ops
    double * normalPtr = normalize(Ax,Ay,Az);
    printf("Ax = 1, Ay = 2, Az = 3 normalized: Bx = %f, By = %f, Bz = %f \n", normalPtr[0], normalPtr[1], normalPtr[2]);
    double dotScalar = dotProduct(Ax,Ay,Az,normalPtr);
    printf("Dot Product of Ax,Ay,Az and Bx,By,Bz = %f \n", dotScalar);
    double * crossPtr = crossProduct(Ax,Ay,Az,normalPtr);
    printf("Cross Product of Ax,Ay,Az and Bx,By,Bz: Cx = %f, Cy = %f, Cz = %f \n", crossPtr[0], crossPtr[1], crossPtr[2]);
    double diffuseTerm = calcDiffuse(Lx,Ly,Lz,Nx,Ny,Nz);
    printf("Diffuse Term = %f \n", diffuseTerm);
    double * idealPtr = calcIdealRef(Lx,Ly,Lz,Nx,Ny,Nz);
    printf("Ideal Reflector of Lx,Ly,Lz and Nx,Ny,Nz: Rx = %f, Ry = %f, Rz = %f \n", idealPtr[0], idealPtr[1], idealPtr[2]);
    double specularTerm = calcSpecular(Vx,Vy,Vz,idealPtr,P);
    printf("Specular Term = %f \n", specularTerm);
}
