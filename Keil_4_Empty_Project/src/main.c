#include "CMSIS/stm32f10x.h"
#include <stdint.h>
#include <math.h>

#define PI 3.14159265358979f
#define GYROMEASERROR PI*(5.f/180.f)
#define BETA sqrt(3.f/4.f)*GYROMEASERROR    

typedef struct quaternion
{
	float w;
	float x;
	float y;
	float z;
} Quaternion;

typedef struct vector
{
	float x;
	float y;
	float z;
} Vector;

Vector quatToEulerAngle(Quaternion q)
{
	Vector ret;
	ret.x = atan2(2.f*q.x*q.y - 2.f*q.w*q.z, 2.f*q.w*q.w + 2.f*q.x*q.x - 1);
	ret.y = -asin(2.f*q.x*q.z + 2.f*q.w*q.y);
	ret.z = atan2(2.f*q.y*q.z - 2.f*q.w*q.x, 2.f*q.w*q.w + 2.f*q.z*q.z - 1);
	return ret;
}

Quaternion scale(Quaternion q, float val)
{
	q.w = val*q.w;
	q.x = val*q.x;
	q.y = val*q.y;
	q.z = val*q.z;
	return q;
}

Quaternion normalize(Quaternion q)
{
	float norm = sqrt(q.w*q.w + q.x*q.x + q.y*q.y +q.z*q.z);
	return scale(q, 1.f/norm);
}

Quaternion inverse(Quaternion q)
{
	q.x = -q.x;
	q.y = -q.y;
	q.z = -q.z;
	return normalize(q);
}

Quaternion mul(Quaternion q1, Quaternion q2)
{
	Quaternion ret;
	ret.w = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
	ret.x = q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
	ret.y = q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x;
	ret.z = q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w;
	return ret;
}

Quaternion summ(Quaternion q1, Quaternion q2)
{
	q1.w += q2.w;
	q1.x += q2.x;
	q1.y += q2.y;
	q1.z += q2.z;
	return q1;
}

Quaternion VecToQuat(Vector v)
{
	return (Quaternion){0.f, v.x, v.y, v.z};
}

Quaternion updateFilterIterator(Quaternion q, Vector a, Vector w, float deltaT)	// q - значение квартерниона на предыдущей итерации, avec - вектор ускорения, wvec -вектор скорости, deltaT - дельта скорости
{
	Quaternion qa = normalize(VecToQuat(a));
	Vector f = {2.f*(q.x*q.z - q.w*q.y) - qa.x,  
							2.f*(q.w*q.x + q.y*q.z) - qa.y,
							2.f*(0.5f - q.x*q.x - q.y*q.y) - qa.z};  // целевая ф-ия    https://habrahabr.ru/post/255661/
	float J[4][3];	// Якобиан
	J[0][0] = -2.f*q.y;
	J[0][1] = 2.f*q.x;
	J[0][2] = 0.f;
	J[1][0] = 2.f*q.z;
	J[1][1] = 2.f*q.w;
	J[1][2] = -4.f*q.x;
	J[2][0] = -2.f*q.w;
	J[2][1] = 2.f*q.z;
	J[2][2] = -4.f*q.y;
	J[3][0] = 2.f*q.x;
	J[3][1] = 2.f*q.y;
	J[3][2] = 0.f;

	Quaternion quatHatDot;		// направление ошибки
	quatHatDot.w = J[0][0]*f.x + J[0][1]*f.y + J[0][2]*f.z;
	quatHatDot.x = J[1][0]*f.x + J[1][1]*f.y + J[1][2]*f.z;
	quatHatDot.y = J[2][0]*f.x + J[2][1]*f.y + J[2][2]*f.z;
	quatHatDot.z = J[3][0]*f.x + J[3][1]*f.y + J[3][2]*f.z;
	quatHatDot = normalize(quatHatDot);
	
	Quaternion quatDotOmega = normalize(mul(scale(q, 0.5f), normalize(VecToQuat(w))));  // скорость изменения ориентации
	q = normalize(summ(q, scale(summ(quatDotOmega, scale(quatHatDot, -BETA)), deltaT)));   // итоговый квартернион
	
	return q;
}

int main(void)       		
{
	Quaternion q = {1.f, 0.f, 0.f, 0.f};									
	Vector a = {1.f, 1.f, 1.f};
	Vector w = {1.f, 1.f, 1.f};     
	Vector ang;
	while(1)
	{
		q = updateFilterIterator(q, a, w, 0.001);	
		ang = quatToEulerAngle(q);
	}
	
}
