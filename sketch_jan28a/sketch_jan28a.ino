#include <Wire.h>
#include <stdint.h>
#define PI 3.14159265358979f
#define G  9.80665f
#define GYROMEASERROR PI*(5.f/180.f)
#define BETA 1.5f*sqrt(3.f/4.f)*GYROMEASERROR
#define ACSADRESS   0x53
#define A_POWER_CTL 0x2D  //Power Control Register
#define A_DATA_FORMAT   0x31
#define A_POINT_DATA 0x32
#define GIROADRESS 0x68
#define G_DLPF_FS 0x16 // full-scale
#define G_PWR_MGM 0x3E //
#define G_POINT_DATA 0x1D
#define g_offx 120
#define g_offy 20
#define g_offz 93

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
  //ret.x = atan2(2.f*q.x*q.y - 2.f*q.w*q.z, 2.f*q.w*q.w + 2.f*q.x*q.x - 1.f);
  //ret.y = asin(2.f*q.x*q.z + 2.f*q.w*q.y);
  //ret.z = atan2(2.f*q.y*q.z - 2.f*q.w*q.x, 2.f*q.w*q.w + 2.f*q.z*q.z - 1.f);
  ret.x = atan2(2.f*(q.w*q.z + q.x*q.y), q.x*q.x + q.w*q.w - q.z*q.z - q.y*q.y);
  ret.y = asin(2.f*(q.x*q.z - q.w*q.y));
  ret.z = atan2(2.f*(q.w*q.x + q.y*q.z), q.z*q.z - q.y*q.y - q.x*q.x + q.w*q.w);
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

Quaternion updateFilterIterator(Quaternion q, Vector a, Vector w, float deltaT) // q - значение квартерниона на предыдущей итерации, avec - вектор ускорения, wvec -вектор скорости, deltaT - дельта скорости
{
  Quaternion qa = normalize(VecToQuat(a));
  Vector f = {2.f*(q.x*q.z - q.w*q.y) - qa.x,  
              2.f*(q.w*q.x + q.y*q.z) - qa.y,
              2.f*(0.5f - q.x*q.x - q.y*q.y) - qa.z};  // целевая ф-ия    https://habrahabr.ru/post/255661/
  float J[4][3];  // Якобиан
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

  Quaternion quatHatDot;    // направление ошибки
  quatHatDot.w = J[0][0]*f.x + J[0][1]*f.y + J[0][2]*f.z;
  quatHatDot.x = J[1][0]*f.x + J[1][1]*f.y + J[1][2]*f.z;
  quatHatDot.y = J[2][0]*f.x + J[2][1]*f.y + J[2][2]*f.z;
  quatHatDot.z = J[3][0]*f.x + J[3][1]*f.y + J[3][2]*f.z;
  quatHatDot = normalize(quatHatDot);

  
  Quaternion quatDotOmega = mul(scale(q, 0.5f), (VecToQuat(w)));  // скорость изменения ориентации
  q.w = q.w + (quatDotOmega.w - BETA*quatHatDot.w)*deltaT;
  q.x = q.x + (quatDotOmega.x - BETA*quatHatDot.x)*deltaT;
  q.y = q.y + (quatDotOmega.y - BETA*quatHatDot.y)*deltaT;
  q.z = q.z + (quatDotOmega.z - BETA*quatHatDot.z)*deltaT;
  q = normalize(q);
  //q = normalize(summ(q, scale(summ(quatDotOmega, scale(quatHatDot, -BETA)), deltaT)));   // итоговый квартернион
  return q;
}

void writeTo(uint8_t device, uint8_t registerAdress, uint8_t value)
{
   Wire.beginTransmission(device);
   Wire.write(registerAdress);
   Wire.write(value);
   Wire.endTransmission();  
}

void readFrom(uint8_t device, uint8_t registerAdress, uint8_t* buff, uint8_t num)
{
  Wire.beginTransmission(device);
  Wire.write(registerAdress); 
  Wire.endTransmission();

  Wire.beginTransmission(device);
  Wire.requestFrom(device, num);
  uint8_t i = 0;
  while(Wire.available())    
  { 
    buff[i] = Wire.read(); 
    i++;
  }
  Wire.endTransmission();   
}

void initGyro()
{
  writeTo(GIROADRESS, G_PWR_MGM, 0x00);
  writeTo(GIROADRESS, G_DLPF_FS, 0x1E);
}

void initAcsel()
{
  writeTo(ACSADRESS, A_POWER_CTL, 0x08);
  writeTo(ACSADRESS, A_DATA_FORMAT, 0x01);  
}

void getGyroData(int32_t* result)
{
  uint8_t buff[6];
  readFrom(GIROADRESS, G_POINT_DATA, buff, 6);
  result[0] = (int32_t)((buff[0] << 8) | buff[1]);
  result[1] = (int32_t)((buff[2] << 8) | buff[3]);
  result[2] = (int32_t)((buff[4] << 8) | buff[5]);  
}

void getAcselData(int32_t* result)
{
  uint8_t buff[6];
  readFrom(ACSADRESS, A_POINT_DATA, buff, 6);
  result[0] = (int32_t)((buff[1] << 8) | buff[0]);
  result[1] = (int32_t)((buff[3] << 8) | buff[2]);
  result[2] = (int32_t)((buff[5] << 8) | buff[4]);  
}


void setup() 
{
  Serial.begin(9600);
  Wire.begin();
  initGyro(); 
  initAcsel(); 
}

Quaternion q = {1.f, 0.f, 0.f, 0.f};
Quaternion quatRotationYaw;
Quaternion quatRotationPitch;
Quaternion quatRotationRoll;
Quaternion quatRotation;
Quaternion aGlobal;
Quaternion VGlobal = {0.f, 0.f, 0.f, 0.f};
Quaternion SGlobal = {0.f, 0.f, 0.f, 0.f};

void loop() 
{
  float wx,wy,wz;
  float ax,ay,az;
  int32_t gyroData[3];
  int32_t acselData[3];
  getGyroData(gyroData);
  getAcselData(acselData);
  wx = PI/180.f*(gyroData[0] / 14.375 + 7.4101);
  wy = PI/180.f*(gyroData[1] / 14.375 + 1.74);
  wz = PI/180.f*(gyroData[2] / 14.375 + 0.42);
  ax = acselData[0] / 128.0;
  ay = acselData[1] / 128.0;
  az = acselData[2] / 128.0; 

  q = updateFilterIterator(q, (Vector){ax, ay, az}, (Vector){wx, wy, wz}, 0.013f);
  Vector temp = quatToEulerAngle(q);
  quatRotationYaw = {cos(temp.x/2), 0.f, 0.f, sin(temp.x/2)};
  quatRotationPitch = {cos(temp.y/2), 0.f, sin(temp.y/2), 0.f};
  quatRotationRoll = {cos(temp.z/2), -sin(temp.z/2), 0.f, 0.f};
  quatRotation = mul(quatRotationRoll, mul(quatRotationPitch, quatRotationYaw));
  //quatRotation = mul(quatRotationRoll, quatRotationYaw);

  aGlobal = mul(inverse(quatRotation), mul((Quaternion){0.f, ax, ay, az}, quatRotation));
  aGlobal = (Quaternion){0.f, aGlobal.x, aGlobal.y, aGlobal.z - 1.f};
  SGlobal = summ(summ(SGlobal, scale(VGlobal, 0.013f)), scale(aGlobal, G*0.0000845f));
  VGlobal = summ(VGlobal, scale(aGlobal, G*0.013f));
  
  //Serial.print("\n yaw=");
  //Serial.print( temp.x );
  //Serial.print("\n pitch=");
  //Serial.println( 180*temp.z/PI );
  //Serial.println( aGlobal.z - 1.0 );
  Serial.println( aGlobal.x);
  //Serial.println( aGlobal.y);
  //Serial.print("\n roll=");
  //Serial.print( temp.z );
  delay(10);
}
