#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <vector>
#include <map>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include "TextureMap.h"
#include "ModelTriangle.h"
#include <algorithm>
#include <string>
#include <iostream>
#include <cmath>
#include "RayTriangleIntersection.h"


using namespace std;


#define WIDTH 320//960//320
#define HEIGHT 240//720//240
#define PI 3.1415926


float buffer[HEIGHT][WIDTH];
void zBuffer(){
    for(auto & y : buffer)
        for(float & x : y)
            x = INT32_MIN;
}
glm::vec3 sphereCameraPosition (0.0, 0.7, 5.2);
glm::vec3 sphereLightPosition (-0.3, 1.3, 1.5);
glm::vec3 boxCameraPosition (0.0, 0.0, 4.0);
glm::mat3 cameraOrientation (1, 0, 0, 0, 1, 0, 0, 0, 1);
glm::vec3 lightPosition (0,0.5,0.5);
glm::vec3 cameraPositionChange(0.0, 0.0, 4.0);
glm::vec3 LP(0, 0.5, 0.5);
TextureMap jupiter("Jupiter.ppm");
TextureMap logoTexture("logo-texture.ppm");
glm::vec3 sphereCenterPoint;
float sphereRadius;

float scalar = 2.0f*HEIGHT/3;
float rayTraceScalar = 50.0f;
float scalingFactor = 0.35;
float focalLength = 2.0;
float ambientStrength = 0.4;

//void draw(DrawingWindow &window) {
//	window.clearPixels();
//	for (size_t y = 0; y < window.height; y++) {
//		for (size_t x = 0; x < window.width; x++) {
//			float red = rand() % 256;
//			float green = 0.0;
//			float blue = 0.0;
//			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
//			window.setPixelColour(x, y, colour);
//		}
//	}
//}

std::vector<float> interpolateSingleFloats(float from, float to, int N){
    vector<float> v;
	float next = from;
	auto a = float((to - from) / float (N - 1));
	for (int i = 0; i < N; i++) {
        v.push_back(next);
        next = next + a;
        //cout<<next<<endl;
	}
	return v;
}

void greyScaleInterpolation(DrawingWindow &window) {
	std::vector<float> result = interpolateSingleFloats(255,0,int(window.width));
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = result[x];
			float green = result[x];
			float blue = result[x];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

std::vector<glm::vec2> interpolateTwoElementValues(glm::vec2 from,glm::vec2 to, int N){
	vector<glm::vec2> v;
	glm::vec2 next = from;
	glm::vec2 a = (to - from) / (float(N- 1));
	 for (int i = 0; i < N; i++) {
		v.push_back(next);
		next = next + a;
	 }
	return v;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from,glm::vec3 to, int N){
	vector<glm::vec3> v;
	glm::vec3 next = from;
	glm::vec3 a = (to - from) / (float(N- 1)); 
	 for (int i = 0; i < N; i++) {
		v.push_back(next);
		next = next + a;
	 }
	return v;
}

void dimensionalColourInterpolation(DrawingWindow &window){
	window.clearPixels();
	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
    glm::vec3 bottomRight(0, 255, 0);    // green 
    glm::vec3 bottomLeft(255, 255, 0);   // yellow
	std::vector<glm::vec3> leftMost;
	std::vector<glm::vec3> rightMost;
	std::vector<glm::vec3> eachRow;
	leftMost = interpolateThreeElementValues(topLeft, bottomLeft, int(window.height));
	rightMost = interpolateThreeElementValues(topRight, bottomRight, int(window.height));
	for (size_t y = 0; y < window.height; y++) {
		eachRow = interpolateThreeElementValues(leftMost[y],rightMost[y],int(window.width));
		for (size_t x = 0; x < window.width; x++) {
			uint32_t colour = (255 << 24) + (int(eachRow[x][0]) << 16) + (int(eachRow[x][1]) << 8) + int(eachRow[x][2]);
			window.setPixelColour(x, y, colour);
		}
	}
}

void lineDraw(DrawingWindow &window, CanvasPoint from, CanvasPoint to, const Colour& c){
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
    float zDiff = to.depth - from.depth;
	float numberOfSpace = max(abs(xDiff),abs(yDiff));
	float xStepSize = xDiff/numberOfSpace;
	float yStepSize = yDiff/numberOfSpace;
    float depthStepSize = zDiff/numberOfSpace;
	int red = c.red;
	int green = c.green;
	int blue = c.blue;
	uint32_t color = (255 << 24) + (red << 16) + (green << 8) + blue;
	for (int i =0; float(i)<numberOfSpace; i++ ){
		float x = from.x + (xStepSize*float(i));
		float y = from.y + (yStepSize*float(i));
        float z = from.depth + (depthStepSize*float(i));
		if(int(round(y))<HEIGHT && int(round(y))>=0 && int(round(x))<WIDTH && int(round(x))>=0){
			if(buffer[int(round(y))][int(round(x))] <= z){
				buffer[int(round(y))][int(round(x))] = z;
				window.setPixelColour(int(round(x)), int(round(y)), color);
			}	
		}
	}
}

//draws "stroked" (unfilled) triangles
CanvasTriangle strokedTriangles(){
	CanvasTriangle triangle;
	for (int i =0; i<3; i++){
		auto random_x = float(rand()%320);
		auto random_y = float(rand()%240);
		if (random_x == random_y){
			random_x = float(rand()%320);
		}
		triangle.vertices[i] = {random_x, random_y};
	}
	return triangle;
}


void drawTriangle(DrawingWindow &window,CanvasTriangle triangle){
    lineDraw(window,triangle.v0(),triangle.v1(),{255,255,255});
	lineDraw(window,triangle.v1(),triangle.v2(),{255,255,255});
	lineDraw(window,triangle.v2(),triangle.v0(),{255,255,255});
}

CanvasTriangle arrangeTriangle(CanvasTriangle triangle){
	if(triangle.v0().y > triangle.v2().y) swap(triangle.v0(), triangle.v2());
	if(triangle.v0().y>triangle.v1().y) swap(triangle.v0(), triangle.v1());
	if(triangle.v1().y > triangle.v2().y) swap(triangle.v1(), triangle.v2());
	return triangle;
}


void filledTriangle(DrawingWindow &window,CanvasTriangle triangle, const Colour& c){
	CanvasTriangle arranged = arrangeTriangle(triangle);
	CanvasPoint top = arranged.v0();
	CanvasPoint bottom = arranged.v2();
	CanvasPoint mid = arranged.v1();	
	for (float y = top.y; y<=mid.y;y++){
	 	float x_left = top.x-((top.x-bottom.x)*(y-top.y)/(bottom.y-top.y));
	 	float x_right = top.x+((y-top.y)*(mid.x-top.x)/(mid.y-top.y));
         //get buffer
	 	lineDraw(window, CanvasPoint(x_left,y), CanvasPoint(x_right,y), c);
	}
	for (float y=mid.y; y<bottom.y;y++) {
        //cout << ((y-mid.y)*(mid.x-bottom.x)/(bottom.y-mid.y)) << endl;
        float x_left = top.x - ((top.x - bottom.x) * (y - top.y) / (bottom.y - top.y));
        float x_right = mid.x - ((y - mid.y) * (mid.x - bottom.x) / (bottom.y - mid.y));
        lineDraw(window, CanvasPoint(x_left, y), CanvasPoint(x_right, y), c);
        // 	for(int x=x_left; x<x_right; x++) window.setPixelColour(x,y,color);
    }
}
unsigned int setColour(const Colour& c){
    int red = c.red;
    int green = c.green;
    int blue = c.blue;
    uint32_t color = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
    return color;
}


unsigned int textureMapperColor(DrawingWindow &window, CanvasTriangle triangle, CanvasPoint point,glm::mat3x3 affineT, TextureMap textureMap) {
	CanvasPoint A = triangle.v0();
	CanvasPoint B = triangle.v1();
	CanvasPoint C = triangle.v2();
	float alpha = (-(point.x-B.x)*(C.y-B.y)+(point.y-B.y)*(C.x-B.x))/(-(A.x-B.x)*(C.y-B.y)+(A.y-B.y)*(C.x-B.x));
	float beta = (-(point.x-C.x)*(A.y-C.y)+(point.y-C.y)*(A.x-C.x))/(-(B.x-C.x)*(A.y-C.y)+(B.y-C.y)*(A.x-C.x));
	float gamma = 1-alpha-beta;
	CanvasPoint texturePoint;
	texturePoint.x = alpha*A.texturePoint.x+beta*B.texturePoint.x+gamma*C.texturePoint.x;
	texturePoint.y = alpha*A.texturePoint.y+beta*B.texturePoint.y+gamma*C.texturePoint.y;
	int index = int(texturePoint.y) * int(textureMap.width) + int(texturePoint.x);
	uint32_t colour = textureMap.pixels[index - 1];
	return colour;
	// 求point的重力坐标
	// 对应到texture的三角里面
	// 知道对应坐标 定位pixel colour
	// return 颜色
}


glm::mat3x3 affine(CanvasTriangle t){
	CanvasTriangle triangle = arrangeTriangle(t);
	float x0,y0,x1,y1,x2,y2,u0,v0,u1,v1,u2,v2;
	x0 = triangle.v0().x;
	y0 = triangle.v0().y;
	x1 = triangle.v1().x;
	y1 = triangle.v1().y;
	x2 = triangle.v2().x;
	y2 = triangle.v2().y;
	u0 = triangle.v0().texturePoint.x;
	v0 = triangle.v0().texturePoint.y;
	u1 = triangle.v1().texturePoint.x;
	v1 = triangle.v1().texturePoint.y;
	u2 = triangle.v2().texturePoint.x;
	v2 = triangle.v2().texturePoint.y;
	glm::mat3x3 texture = {{u0,u1,u2},{v0,v1,v2},{1,1,1}};
	glm::mat3x3 oTriangle = {{x0,x1,x2},{y0,y1,y2},{1,1,1}};
	glm::mat3x3 iTriangle = glm::inverse(oTriangle);
	glm::mat3x3 affineT = iTriangle * texture;
	//std::cout<<glm::to_string(affineT)<<std::endl;
	return affineT;
}





void affineTransformation(DrawingWindow &window, CanvasTriangle triangle, CanvasPoint point,glm::mat3x3 affineT, TextureMap textureMap){
	//inverse of triangle
	TexturePoint p;
	p.x = round(affineT[0][0]*point.x+affineT[0][1]*point.y+affineT[0][2]);
	p.y = round(affineT[1][0]*point.x+affineT[1][1]*point.y+affineT[1][2]);
	//std::cout << p.x << p.y << std::endl;
	//std::cout<<glm::to_string(affineT)<<std::endl;
    int num = (int(p.y)*(int(textureMap.width))+int(p.x));
	uint32_t c = textureMap.pixels[num-1];
	window.setPixelColour(int(round(point.x)),int(round((point.y))),c);
	//return c; 
}

void textureMapper(DrawingWindow &window, CanvasTriangle arranged,glm::mat3x3 affineT, const TextureMap& textureMap){
	CanvasPoint top = arranged.v0();
	CanvasPoint bottom = arranged.v2();
	CanvasPoint mid = arranged.v1();
	for (float y = top.y; y<=mid.y;y++){
		float x_left = top.x-((top.x-bottom.x)*(y-top.y)/(bottom.y-top.y));
		float x_right = top.x+((y-top.y)*(mid.x-top.x)/(mid.y-top.y));
		if (x_left > x_right) std::swap(x_left,x_right);
		for(float x=x_left; x<=x_right; x++) {
			affineTransformation(window, arranged, CanvasPoint(x, y), affineT, textureMap);
			//window.setPixelColour(x,y,color);
		}
	}
	for (float y=mid.y; y<=bottom.y;y++){
		float x_left = top.x-((top.x-bottom.x)*(y-top.y)/(bottom.y-top.y));
		float x_right = mid.x-((y-mid.y)*(mid.x-bottom.x)/(bottom.y-mid.y));
		if (x_left > x_right) std::swap(x_left,x_right);
		for(float x=x_left; x<=x_right; x++) {
			affineTransformation(window, arranged, CanvasPoint(x, y), affineT, textureMap);
			//window.setPixelColour(x,y,color);
		}
	}
}

//The type map returns a matching vector with name on the left and value on the right, so you can use the type map to match the name of a colour to get rgb
std::map<string,Colour> readMaterialFile(const string& fileName){
 ifstream MyReadFile(fileName);
 string fileText;
 std::map<string,Colour> palette;
 Colour c;
 while(getline(MyReadFile,fileText)){
    std::vector<string> splitDelimiter = split(fileText,' ');
    //std::cout<<splitDelimiter[0]<<std::endl;
    if(splitDelimiter[0]=="newmtl"){
		c.name = splitDelimiter[1];
    }
	if(splitDelimiter[0]=="Kd"){
		c.red = int(round(std::stof(splitDelimiter[1]) * 255));
		c.green = int(round(std::stof(splitDelimiter[2]) * 255));
		c.blue = int(round(std::stof(splitDelimiter[3]) * 255));
		palette[c.name]=c;
		// std::cout<<c<<std::endl;
	} 
 }
 MyReadFile.close();
 return palette;
}


//把世界坐标变成相机坐标
glm::vec3 getCameraPoint(glm::vec3 vertexPoint,glm::vec3 cameraPosition){
	vertexPoint.x = vertexPoint.x-cameraPosition.x;
	vertexPoint.y = vertexPoint.y-cameraPosition.y;
	vertexPoint.z = vertexPoint.z-cameraPosition.z;
	return vertexPoint;
}

CanvasPoint getCanvasIntersectionPoint (glm::vec3 cameraPosition, glm::vec3 vertexPosition, float FL, glm::mat3x3 orientation){
	CanvasPoint modelVertex;
	int imagePlaneScaling = HEIGHT *2/3;
	glm::vec3 c = getCameraPoint(vertexPosition,cameraPosition);
	c = c * orientation;
	modelVertex.x = round(float(-imagePlaneScaling) * FL * (c.x / c.z) + float(WIDTH) / 2);
	modelVertex.y = round(float(imagePlaneScaling) * FL * (c.y / c.z) + float(HEIGHT) / 2);
    modelVertex.depth =abs(1/c.z);
	return modelVertex;
}

void pointCloud(DrawingWindow &window, glm::vec3 cameraPosition, float FL, const std::vector<ModelTriangle>& triangle){
	for(auto & i : triangle){
		for(auto v : i.vertices){
				uint32_t color = (255 << 24) + (int(255) << 16) + (int(255) << 8) + int(255);
		    CanvasPoint p=getCanvasIntersectionPoint(cameraPosition, v, FL, cameraOrientation);
			//std::cout<<p.x<<' '<<p.y<<std::endl;
			window.setPixelColour(int(round(p.x)),int(round(p.y)),color);
		}
 	}
}

void wireFrameRender(DrawingWindow &window, glm::vec3 cameraPosition, const std::vector<ModelTriangle>& triangle, float FL, glm::mat3x3 orientation){
	for(ModelTriangle t : triangle){
		CanvasPoint v0=getCanvasIntersectionPoint(cameraPosition, t.vertices[0], FL, orientation);
		CanvasPoint v1=getCanvasIntersectionPoint(cameraPosition, t.vertices[1], FL, orientation);
		CanvasPoint v2=getCanvasIntersectionPoint(cameraPosition, t.vertices[2], FL, orientation);
		CanvasTriangle ct = CanvasTriangle(v0,v1,v2);
		drawTriangle(window,ct);
	}
}

void setBufferTriangle(DrawingWindow &windows,CanvasPoint top,CanvasPoint v1,CanvasPoint v2, const Colour& c){
    //一个三角形被分成两个平底三角形，判断是往top画还是往bottom画
    float d = abs(v1.y-top.y);
    CanvasPoint left=top;
    CanvasPoint right=top;
    float s=1.0f;
    if(top.y - v2.y > 0) s=-1.0f;
    for (int i =1; float(i)<=d; i++){
        float radio = float(i)/abs(d);
        left.x=top.x+((v1.x-top.x)*(radio));
        right.x = top.x+((v2.x-top.x)*(radio));
        left.y += s;
        right.y += s;
        left.depth = top.depth + ((v1.depth - top.depth)*radio);
        right.depth = top.depth + ((v2.depth - top.depth)*radio);
        lineDraw(windows, left, right, c);
    }
}

void occlusion(DrawingWindow &window,CanvasTriangle triangle, const Colour& colour){
    CanvasTriangle arranged = arrangeTriangle(triangle);
    CanvasPoint top = arranged.v0();
    CanvasPoint bottom = arranged.v2();
    CanvasPoint mid = arranged.v1();
    double radio = (bottom.x - top.x)/(bottom.y-top.y);
    float extra_x;
    float extra_y;
    float extraDepth;
    CanvasPoint extraPoint;
    if(radio>=0) extra_x = float(top.x+((mid.y-top.y)*abs(radio)));
    else extra_x = float(top.x-((mid.y-top.y)*abs(radio)));
    extra_y =mid.y;
    float alpha = (-(extra_x-triangle.v1().x)*(triangle.v2().y-triangle.v1().y)+(extra_y-triangle.v1().y)*(triangle.v2().x-triangle.v1().x))/(-(triangle.v0().x-triangle.v1().x)*(triangle.v2().y-triangle.v1().y)+(triangle.v0().y-triangle.v1().y)*(triangle.v2().x-triangle.v1().x));
    float beta = (-(extra_x-triangle.v2().x)*(triangle.v0().y-triangle.v2().y)+(extra_y-triangle.v2().y)*(triangle.v0().x-triangle.v2().x))/(-(triangle.v1().x-triangle.v2().x)*(triangle.v0().y-triangle.v2().y)+(triangle.v1().y-triangle.v2().y)*(triangle.v0().x-triangle.v2().x));
    float gamma = 1 - alpha - beta;
    extraDepth = alpha * triangle.v0().depth + beta * triangle.v1().depth + gamma * triangle.v2().depth;
    extraPoint = CanvasPoint(extra_x, extra_y, extraDepth);
    setBufferTriangle(window,top,mid,extraPoint,colour);
    setBufferTriangle(window,bottom,mid,extraPoint,colour);
}

void rasterisedRender(DrawingWindow &window, const std::vector<ModelTriangle>& triangle, glm::vec3 cameraPosition, float FL, glm::mat3x3 oritation){
	for(ModelTriangle t : triangle){
		CanvasPoint v0=getCanvasIntersectionPoint(cameraPosition, t.vertices[0], FL, oritation);
		CanvasPoint v1=getCanvasIntersectionPoint(cameraPosition, t.vertices[1], FL, oritation);
		CanvasPoint v2=getCanvasIntersectionPoint(cameraPosition, t.vertices[2], FL, oritation);
		CanvasTriangle ct = CanvasTriangle(v0,v1,v2);
		//filledTriangle(window,ct,t.colour);
        occlusion(window,ct,t.colour);
    }
}
//void oribit(DrawingWindow &window,glm::mat3x3 x,std::vector<ModelTriangle> l,glm::vec3 boxCameraPosition,float focalLength,glm::mat3x3 cameraOrientation){
//    boxCameraPosition = x*boxCameraPosition;
//    cameraOrientation = x*cameraOrientation;
//    rasterisedRender(window,l,boxCameraPosition,focalLength,cameraOrientation);
//}

glm::mat3 lookAt(glm::vec3 cameraPosition){
	glm::vec3 vertical (0,1,0);
	glm::vec3 center (0,0,0);
	glm::vec3 forward = glm::normalize(cameraPosition-center);
	glm::vec3 right = glm::cross(vertical,forward);
	glm::vec3 up = glm::cross(forward,right);
	glm::mat3 newOritation = glm::mat3(right.x, right.y, right.z, up.x, up.y, up.z, forward.x, forward.y, forward.z);
	return newOritation;
}

// Scan the screen from the far left to the right, the camera is placed in world coordinates, so convert the screen to world coordinates
//glm::vec3 getWorldPoint(CanvasPoint p, float FL, float S, glm::vec3 cameraPosition,glm::mat3 orientation){
//	glm::vec3 point;
//	point.z = cameraPosition.z - FL;
//	point.x =(point.z*(p.x-float(WIDTH)/2))/(FL * S);
//	point.y =(point.z*(p.y-float(HEIGHT)/2))/(FL * (-S));
//    //point = point * orientation;
//    point = {point.x+cameraPosition.x,point.y+cameraPosition.y,point.z+(-cameraPosition.z)};
//	return point;
//}

glm::vec3 getWorldPoint(CanvasPoint p, float FL, float S, glm::vec3 cameraPosition,glm::mat3 orientation){
    glm::vec3 point;
//    if(n){
//        point.z = cameraPosition.z - FL;
//        point.x =(point.z*(p.x-float(WIDTH)/2))/(FL * S);
//        point.y =(point.z*(p.y-float(HEIGHT)/2))/(FL * (-S));
//        //point = point * orientation;
//        point = {point.x+cameraPosition.x,point.y+cameraPosition.y,point.z+(-cameraPosition.z)};
   // }else{
        point.z = cameraPosition.z - FL;
        point.x =(p.x-float(WIDTH)/2)/S;
        point.y =(p.y-float(HEIGHT)/2)/-S;
        //point = point * orientation;
        point = {point.x+cameraPosition.x,point.y+cameraPosition.y,point.z};
    //}


    return point;
}

glm::vec3 getSphereWorldPoint(CanvasPoint p, float FL, float S, glm::vec3 cameraPosition,glm::mat3 orientation){
    glm::vec3 point;
    point.z = cameraPosition.z - FL;
    point.x =(point.z*(p.x-float(WIDTH)/2))/(FL * S);
    point.y =(point.z*(p.y-float(HEIGHT)/2))/(FL * (-S));
    point = point * orientation;
    point = {point.x+cameraPosition.x,point.y+cameraPosition.y,point.z+(-cameraPosition.z)};
    return point;
}


RayTriangleIntersection getClosetsIntersection(glm::vec3 cameraPosition,glm::vec3 rayDirection,const std::vector<ModelTriangle>& triangle){
	RayTriangleIntersection r;
	float tMin = FLT_MAX;
	int index=0;
	for(ModelTriangle i : triangle){
		glm::vec3 e0 = i.vertices[1] - i.vertices[0];
		glm::vec3 e1 = i.vertices[2] - i.vertices[0];
		glm::vec3 SPVector = cameraPosition - i.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float t = possibleSolution.x;
		float u = possibleSolution.y;
		float v = possibleSolution.z;
		if (t<tMin && t>=0.000001 && u>=0.000001 && u<=1 && v>=0.000001 && v<=1 && (u + v) <= 1) {
			tMin = t;
			// Update RayTriangleIntersection to the nearest points
			r.distanceFromCamera=tMin;
           // cout<<r.distanceFromCamera<<endl;
			r.intersectionPoint = cameraPosition+tMin*rayDirection;
			r.intersectedTriangle=i;
			r.triangleIndex=index;
		}
		index++;
	}
	return r;
}
/*
void shadow(DrawingWindow &window,float FL,float S,glm::vec3 boxCameraPosition,std::vector<ModelTriangle> triangle,glm::vec3 LP){
	for(int y=0; y<HEIGHT; y++){
		for(int x=0; x<WIDTH; x++){
			glm::vec3 worldPoint = getWorldPoint(CanvasPoint(x,y),FL,S,boxCameraPosition);
			glm::vec3 rayDirection = worldPoint-boxCameraPosition;
			RayTriangleIntersection r = getClosetsIntersection(boxCameraPosition,rayDirection,triangle);
			glm::vec3 lightDirection = r.intersectionPoint-LP;
			RayTriangleIntersection light = getClosetsIntersection(LP,lightDirection,triangle);
			//To determine whether the point seen by the camera is the same as the point illuminated by the light, if it is the same it means that it is not in the shadow.
			if(r.triangleIndex==light.triangleIndex){
				uint32_t color = setColour(r.intersectedTriangle.colour);
				window.setPixelColour(x,y,color);
			}else window.setPixelColour(x,y,setColour(Colour(0,0,0)));
		}
	}
}
*/

vector<glm::vec3> multipleLight(glm::vec3 LP){
    vector<glm::vec3> ml;
    float r = 16.0f;
    float step = 0.05f;
    for(int i=int(-r); i<int(r); i++){
        for(int j=int(-r); j<int(r); j++){
            if(fabs(i)+fabs(j)<=r){
                glm::vec3 light (LP.x+float(i)*step,LP.y,LP.z+float(j)*step);
                ml.push_back(light);
            }
        }
    }
    return ml;
}

vector<glm::vec3> multipleLightPoint = multipleLight(LP);

float softShadow(const RayTriangleIntersection& closetIntersection,const std::vector<ModelTriangle>& triangle,glm::vec3 LP){
    float number = 1.0f;
    for(glm::vec3 l : multipleLightPoint ) {
        glm::vec3 mLightDirection = glm::normalize(closetIntersection.intersectionPoint - l);
        RayTriangleIntersection mLight = getClosetsIntersection(l, mLightDirection, triangle);
        if(mLight.triangleIndex == closetIntersection.triangleIndex) number++;
    }
    number = number/float(multipleLightPoint.size());
    return number;
}


glm::mat3x3 modelTriangleAffine(ModelTriangle triangle,const TextureMap& file){
    float x0,y0,x1,y1,x2,y2,z0,z1,z2,u0,v0,u1,v1,u2,v2;
    x0 = triangle.vertices[0].x;
    y0 = triangle.vertices[0].y;
    z0 = triangle.vertices[0].z;
    x1 = triangle.vertices[1].x;
    y1 = triangle.vertices[1].y;
    z1 = triangle.vertices[1].z;
    x2 = triangle.vertices[2].x;
    y2 = triangle.vertices[2].y;
    z2 = triangle.vertices[2].z;
    u0 = triangle.texturePoints[0].x * float (file.width);
    v0 = triangle.texturePoints[0].y * float (file.height);
    u1 = triangle.texturePoints[1].x * float (file.width);
    v1 = triangle.texturePoints[1].y * float (file.height);
    u2 = triangle.texturePoints[2].x * float (file.width);
    v2 = triangle.texturePoints[2].y * float (file.height);

    glm::mat3x3 texture = {{u0,u1,u2},{v0,v1,v2},{1,1,1}};
    glm::mat3x3 oTriangle = {{x0,x1,x2},{y0,y1,y2},{z0,z1,z2}};
    glm::mat3x3 iTriangle = glm::inverse(oTriangle);
    glm::mat3x3 affineT = iTriangle * texture ;
    //std::cout<<glm::to_string(affineT)<<std::endl;
    return affineT;
}

uint32_t textureModelTriangle(glm::vec3 texturePoint, const TextureMap& file){
    int index = int(texturePoint.y) * int(file.width) + int(texturePoint.x);
    uint32_t c = file.pixels[index-1];
    return  c;
}

Colour textureMapFloor(const TextureMap& file, const RayTriangleIntersection& closetIntersection){
    glm::mat3x3 affine = modelTriangleAffine(closetIntersection.intersectedTriangle,file);
    glm::vec3 texturePoint = closetIntersection.intersectionPoint * affine;
    Colour c;
    uint32_t color = textureModelTriangle(texturePoint,file);
    c.blue = int((color) & 0xff);
    c.green = int((color >> 8) & 0xff);
    c.red = int((color >> 16) & 0xff);
    return c;
}

TextureMap file = TextureMap("chessboard.ppm");

float fresnel(glm::vec3 incident, glm::vec3 normal, float refractiveIndex) {
    float cosi = glm::dot(incident, normal);
    float etai = 1;
    float etat = refractiveIndex;
    if (cosi > 0){
        std::swap(etai, etat);
    }
    float eta = etai / etat;
    float sint = eta * sqrtf(std::max(0.f, 1 - cosi * cosi));
    if (sint >= 1){
        return 1;
    } else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (Rs * Rs + Rp * Rp) * 0.5f;
    }
    return 0;
}

glm::vec3 vertexNormal(const std::vector<ModelTriangle>& triangle,glm::vec3 point){
    glm::vec3 vertexN;
    float count = 0;
    glm::vec3 verticeSum;
    for(ModelTriangle i : triangle) {
        if(i.vertices[0] == point || i.vertices[1] == point || i.vertices[2] == point) {
            verticeSum+=i.normal;
            count++;
        }
    }
    vertexN = verticeSum * (1/count);
    return glm::normalize(vertexN);
}

Colour drawPhoneSphereInBox(RayTriangleIntersection closetIntersection, const std::vector<ModelTriangle>& triangle, glm::vec3 cameraPosition, glm::vec3 LP, Colour c){
    float specularExponent = 256.0f;
    float sourceStrength = 8;
    glm::vec3 v0;
    glm::vec3 v1;
    glm::vec3 v2;

    v0 = closetIntersection.intersectedTriangle.vertices[0];
    v1 = closetIntersection.intersectedTriangle.vertices[1];
    v2 = closetIntersection.intersectedTriangle.vertices[2];

    glm::vec3 lightDirection = glm::normalize(closetIntersection.intersectionPoint - LP);
    float lightDistance = glm::length(LP - closetIntersection.intersectionPoint);
    float PL = sourceStrength / float(( 4 *PI * lightDistance * lightDistance));

    glm::vec3 point = closetIntersection.intersectionPoint;
    glm::vec3 n0 = vertexNormal(triangle, closetIntersection.intersectedTriangle.vertices[0]);
    glm::vec3 n1 = vertexNormal(triangle,closetIntersection.intersectedTriangle.vertices[1]);
    glm::vec3 n2 = vertexNormal(triangle,closetIntersection.intersectedTriangle.vertices[2]);

    float alpha = (-(point.x-v1.x)*(v2.y-v1.y)+(point.y-v1.y)*(v2.x-v1.x))/(-(v0.x-v1.x)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.x-v1.x));
    float beta = (-(point.x-v2.x)*(v0.y-v2.y)+(point.y-v2.y)*(v0.x-v2.x))/(-(v1.x-v2.x)*(v0.y-v2.y)+(v1.y-v2.y)*(v0.x-v2.x));
    float gamma = 1-alpha-beta;

    glm::vec3 normal = n0*alpha + n1*beta + gamma*n2;

    float lightAngle = max(0.0f,glm::dot(normal, lightDirection));
    glm::vec3 view = glm::normalize(cameraPosition - closetIntersection.intersectionPoint);
    glm::vec3 reflection = lightDirection - 2.0f*(normal * glm::dot(lightDirection, normal));
    float specularLight = glm::pow(glm::dot(view,reflection),specularExponent);
    float brightness = PL * lightAngle + specularLight + ambientStrength;

    float red = std::min(255.0f, c.red * brightness);
    float green = std::min(255.0f, c.green * brightness);
    float blue = std::min(255.0f, c.blue * brightness);
    //判断点是否在圆内
    if(closetIntersection.intersectionPoint == glm::vec3(0, 0, 0)) {
        red = 0;
        green = 0;
        blue = 0;
    }
    return Colour(int(round(red)), (int(round(green))), int(round(blue)));
}

Colour shootRay(glm::vec3 cameraPosition, glm::vec3 rayDirection, const std::vector<ModelTriangle>& triangle,glm::vec3 LP,float sourceStrength, int index) {
    if(index >= 5) {return {255, 255, 255};}
    float specularExponent = 256.0f;
    RayTriangleIntersection closetIntersection = getClosetsIntersection(cameraPosition, rayDirection, triangle);
    glm::vec3 lightDirection = glm::normalize(closetIntersection.intersectionPoint - LP);
    RayTriangleIntersection light = getClosetsIntersection(LP, lightDirection, triangle);
    //model表面的点到灯的方向向量
    float lightDistance = glm::length(LP - closetIntersection.intersectionPoint);
    float PL = sourceStrength / float((4 * PI * lightDistance * lightDistance));
    float lightAngle = max(0.0f, glm::dot(closetIntersection.intersectedTriangle.normal, lightDirection));
    glm::vec3 view = glm::normalize(closetIntersection.intersectionPoint - cameraPosition);
    glm::vec3 reflection = lightDirection - 2.0f * (closetIntersection.intersectedTriangle.normal * glm::dot(lightDirection,closetIntersection.intersectedTriangle.normal));
    float specularLight = glm::pow(glm::dot(view, reflection), specularExponent);
    float lightIntensity = PL * lightAngle + specularLight + ambientStrength;
    float num = softShadow(closetIntersection, triangle, LP);
    Colour c = closetIntersection.intersectedTriangle.colour;
    if (closetIntersection.intersectedTriangle.colour.name == "Yellow") {
        glm::vec3 mirrorReflection = glm::normalize(rayDirection-2.0f * closetIntersection.intersectedTriangle.normal * glm::dot(rayDirection,closetIntersection.intersectedTriangle.normal));
        return shootRay(closetIntersection.intersectionPoint,mirrorReflection,triangle,LP,sourceStrength, index+1);
    } else if(closetIntersection.intersectedTriangle.colour.name=="Cobbles"){
        //closetIntersection.intersectedTriangle.colour.name=="Cobbles"
        c = textureMapFloor(file,closetIntersection);
        //cout<<1<<endl;
    } else if (closetIntersection.intersectedTriangle.colour.name == "Blue") {
        // glass = reflection + refraction
        // reflection
        glm::vec3 glassReflection = glm::normalize(rayDirection - 2.0f * closetIntersection.intersectedTriangle.normal * glm::dot(rayDirection,closetIntersection.intersectedTriangle.normal));
        glm::vec3 updatePoint = closetIntersection.intersectionPoint;
        Colour reflectionColour = shootRay(closetIntersection.intersectionPoint, glassReflection, triangle, LP, sourceStrength, index + 1);
        // refraction
        // Calculate the refraction direction
        float air = 1.0f;
        float glass = 1.3f;
        glm::vec3 incidentRay = rayDirection;
        glm::vec3 normal = closetIntersection.intersectedTriangle.normal;
        glm::vec3 direction (0, 0, 0);
        float angle = std::fmax(std::fmin(glm::dot(incidentRay, normal), 1.0f), -1.0f);
        float ratio = 0.0f;
        float k = 0.0f;
        if(angle > 0) {
            std::swap(glass, air);
            normal = -normal;
        }
        ratio = air / glass;
        k = 1 - ratio * ratio * (1 - abs(angle) * abs(angle));
        if(k >= 0) {
            direction = glm::normalize(incidentRay * ratio + normal * (ratio * angle - sqrt(k)));
        }
        glm::vec3 adjustPoint = closetIntersection.intersectionPoint;
        if(angle < 0) {
            adjustPoint = closetIntersection.intersectionPoint - closetIntersection.intersectedTriangle.normal * float(0.00001);
        } else {
            adjustPoint = closetIntersection.intersectionPoint + closetIntersection.intersectedTriangle.normal * float(0.00001);
        }
        Colour refractionColour = shootRay(adjustPoint, direction, triangle, LP, sourceStrength, index + 1);
        // mix colour
        // fresnel
        float reflectiveConstant = fresnel(rayDirection, closetIntersection.intersectedTriangle.normal, 1.5f);
        float refractiveConstant = 1-reflectiveConstant;
        Colour mixColour;
        mixColour.red = int(reflectiveConstant * float(reflectionColour.red) + refractiveConstant * float(refractionColour.red));
        mixColour.green = int(reflectiveConstant * float(reflectionColour.green) + refractiveConstant * float(refractionColour.green));
        mixColour.blue = int(reflectiveConstant * float(reflectionColour.blue) + refractiveConstant * float(refractionColour.blue));
        return mixColour;
    } else if (closetIntersection.intersectedTriangle.colour.name == "SkyBlue") {
        glm::vec3 spherePoint = (closetIntersection.intersectionPoint - sphereCenterPoint) / sphereRadius;
        float PHI = atan2(spherePoint.z, spherePoint.x);
        float THETA = asin(spherePoint.y);
        float u = 1-(PHI + 3.1415) / (2*3.1415);
        float v = (THETA + 3.1415/2) / 3.1415;
        int pointX = int((u)* jupiter.width);
        int pointY = int((1 - v)*jupiter.height - 0.001f);
        uint32_t colour32 = jupiter.pixels[pointX + pointY * jupiter.width];
        int blue = (colour32) & 0xff;
        int green = (colour32 >> 8) & 0xff;
        int red = (colour32 >> 16) & 0xff;
        c = Colour(red, green, blue);
        c = drawPhoneSphereInBox(closetIntersection, triangle, cameraPosition, LP, c);
        return c;
    } else if (closetIntersection.intersectedTriangle.colour.name == "Carrot") {
        c = textureMapFloor(logoTexture,closetIntersection);
        std::array<TexturePoint, 3> triangleTexturePoints = closetIntersection.intersectedTriangle.texturePoints;

        glm::vec3 A = closetIntersection.intersectedTriangle.vertices[0];
        glm::vec3 B = closetIntersection.intersectedTriangle.vertices[1] + 0.0001f;
        glm::vec3 C = closetIntersection.intersectedTriangle.vertices[2] + 0.0001f;
        glm::vec2 textureA = glm::vec2(triangleTexturePoints[0].x, triangleTexturePoints[0].y);
        glm::vec2 textureB = glm::vec2(triangleTexturePoints[1].x, triangleTexturePoints[1].y);
        glm::vec2 textureC = glm::vec2(triangleTexturePoints[2].x, triangleTexturePoints[2].y);
        glm::vec3 normal = closetIntersection.intersectedTriangle.normal;
        //std::cout << glm::to_string(normal) << std::endl;

        glm::vec3 deltaPos1 = B - A;
        glm::vec3 deltaPos2 = C - A;
        glm::vec2 deltaUV1 = textureB - textureA;
        glm::vec2 deltaUV2 = textureC - textureA;

        float r = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV1.y * deltaUV2.x);
        glm::vec3 tangent = (deltaPos1 * deltaUV2.y   - deltaPos2 * deltaUV1.y)*r;
        tangent = glm::normalize(tangent - (glm::dot(tangent, normal) * normal));
        glm::vec3 bitangent = (deltaPos2 * deltaUV1.x   - deltaPos1 * deltaUV2.x)*r;
        bitangent = glm::normalize(bitangent - (glm::dot(bitangent, normal) * normal) - (glm::dot(bitangent, tangent) * tangent));
        normal = glm::normalize(normal);
        glm::mat3x3 tangentMatrix = {tangent, bitangent, normal};
        //std::cout << glm::to_string(tangentMatrix) << std::endl;

        //float pointGray = c.red * 0.3 + c.green * 0.59 + c.blue * 0.11;
        glm::vec3 pixelNormal = glm::vec3(c.red/(float)255, c.green/(float)255, c.blue/(float)255);
        glm::vec3 textureNormal = 2.0f * pixelNormal - 1.0f;
        glm::vec3 worldNormal = textureNormal * tangentMatrix;
        //std::cout << glm::to_string(worldNormal) << std::endl;
        

        lightDistance = glm::length(LP - closetIntersection.intersectionPoint);
        PL = 16 / float((4 * PI * lightDistance * lightDistance));
        lightAngle = glm::dot(worldNormal, lightDirection);
        lightAngle = glm::clamp<float>(lightAngle, 0.0, 1.0);
        view = glm::normalize(closetIntersection.intersectionPoint - cameraPosition);
        reflection = lightDirection - 2.0f * (worldNormal * glm::dot(lightDirection,worldNormal));
        specularLight = glm::pow(fmax(0.0f, glm::dot(view, reflection)), 256);
        lightIntensity = (PL * lightAngle + specularLight + ambientStrength);

    }
    float red = std::max(0.0f, std::min(255.0f, float((c.red)) * lightIntensity));
    float green = std::max(0.0f, std::min(255.0f, float((c.green)) * lightIntensity));
    float blue = std::max(0.0f, std::min(255.0f, float((c.blue)) * lightIntensity));
    auto shadowValue = glm::clamp<float>((0.6f/ambientStrength*num),0,1);
//    auto shadowValue = ambientStrength;

    if (closetIntersection.triangleIndex != light.triangleIndex && closetIntersection.intersectedTriangle.colour.name != "Yellow") {
        if(light.intersectedTriangle.colour.name == "Blue") {
            shadowValue = 0.8;
        }
        red = red * shadowValue;
        green = green * shadowValue;
        blue = blue * shadowValue;
    }

    return {int(round(red)),int(round(green)),int(round(blue))};
}

void ray(DrawingWindow &window, float FL, float S, glm::vec3 cameraPosition, const std::vector<ModelTriangle>& triangle, glm::vec3 LP,glm::mat3x3 orientation){
    float sourceStrength = 16;
	for(int y=0; y<HEIGHT; y++){
		for(int x=0; x<WIDTH; x++){
			glm::vec3 worldPoint = (getWorldPoint(CanvasPoint(float(x),float(y)), FL, 2.0f*HEIGHT/3, cameraPosition, orientation));
            glm::vec3 rayDirection = glm::normalize(worldPoint-cameraPosition)*glm::inverse(orientation);
            Colour c = shootRay(cameraPosition,rayDirection,triangle,LP,sourceStrength, 1);
			uint32_t color = (255 << 24) + (int(round(c.red)) << 16) + (int(round(c.green)) << 8) + int(round(c.blue));
			window.setPixelColour(x,y,color);
		}
	}
}

void drawFlatSphere(DrawingWindow &window, float FL, float S, glm::vec3 cameraPosition, const std::vector<ModelTriangle>& triangle, glm::vec3 LP,glm::mat3x3 orientation){
	float sourceStrength = 8;
	float specularExponent = 256.0f;
	for(int y=0; y<HEIGHT; y++){
		for(int x=0; x<WIDTH; x++){
			glm::vec3 worldPoint = getSphereWorldPoint(CanvasPoint(float(x),float(y)), FL, S, cameraPosition, orientation);
			glm::vec3 rayDirection = glm::normalize(worldPoint-cameraPosition);
			RayTriangleIntersection r = getClosetsIntersection(cameraPosition,rayDirection,triangle);
			glm::vec3 lightDirection = glm::normalize(r.intersectionPoint - LP);
			RayTriangleIntersection light = getClosetsIntersection(LP, lightDirection, triangle);
			//model表面的点到灯的方向向量
			float lightDistance = glm::length(LP - r.intersectionPoint);
			//cout<<lightDistance<<endl;
			float PL = sourceStrength / float( 4 *PI * lightDistance * lightDistance);
			float lightAngle = max(0.0f,glm::dot(r.intersectedTriangle.normal,lightDirection));
			//if(r.intersectedTriangle.normal.x != 0) cout<<r.intersectedTriangle.normal.x<<r.intersectedTriangle.normal.y<<r.intersectedTriangle.normal.z<<endl;
			glm::vec3 view = glm::normalize(cameraPosition-r.intersectionPoint);
			glm::vec3 reflection = lightDirection - 2.0f*(r.intersectedTriangle.normal*glm::dot(lightDirection,r.intersectedTriangle.normal));
			float specularLight = glm::pow(glm::dot(view,reflection),specularExponent);
			//cout<<glm::dot(view,reflection)<<endl;
			float lightIntensity = PL*lightAngle+ specularLight+ambientStrength;
			//cout<<lightDistance<<endl;
			float red = std::min(255.0f,(255)*lightIntensity);
			float green = std::min(255.0f,(0)*lightIntensity);
			float blue = std::min(255.0f,(0)*lightIntensity);
			//判断点是否在圆内
			if(r.intersectionPoint == glm::vec3(0, 0, 0)) {
                red = 0;
                green = 0;
                blue = 0;
            }
			uint32_t color = (255 << 24) + (int(round(red)) << 16) + (int(round(green)) << 8) + int(round(blue));
			window.setPixelColour(x,y,color);
		}
	}
}

void drawGouraudSphere(DrawingWindow &window, float FL, float S, glm::vec3 cameraPosition, const std::vector<ModelTriangle>& triangle, glm::vec3 LP,glm::mat3x3 orientation){
    float sourceStrength = 8;
    float specularExponent = 256;
    for(int y=0; y<HEIGHT; y++){
        for(int x=0; x<WIDTH; x++){
            glm::vec3 v0;
            glm::vec3 v1;
            glm::vec3 v2;

            glm::vec3 worldPoint = getSphereWorldPoint(CanvasPoint(float(x),float (y)), FL, S, cameraPosition, orientation);
            glm::vec3 rayDirection = glm::normalize(worldPoint-cameraPosition);
            RayTriangleIntersection closetIntersection = getClosetsIntersection(cameraPosition, rayDirection, triangle);

            v0 = closetIntersection.intersectedTriangle.vertices[0];
            v1 = closetIntersection.intersectedTriangle.vertices[1];
            v2 = closetIntersection.intersectedTriangle.vertices[2];

            glm::vec3 lightDirection0 = glm::normalize(v0 - LP);
            glm::vec3 lightDirection1 = glm::normalize(v1 - LP);
            glm::vec3 lightDirection2 = glm::normalize(v2 - LP);

            float lightDistance0 = glm::length(LP - v0);
            float lightDistance1 = glm::length(LP - v1);
            float lightDistance2 = glm::length(LP - v2);

            float PL0 = sourceStrength / float(( 4 *PI * lightDistance0 * lightDistance0));
            float PL1 = sourceStrength / float(( 4 *PI * lightDistance1 * lightDistance1));
            float PL2 = sourceStrength / float(( 4 *PI * lightDistance2 * lightDistance2));

            glm::vec3 point = closetIntersection.intersectionPoint;
            glm::vec3 n0 = vertexNormal(triangle, closetIntersection.intersectedTriangle.vertices[0]);
            glm::vec3 n1 = vertexNormal(triangle,closetIntersection.intersectedTriangle.vertices[1]);
            glm::vec3 n2 = vertexNormal(triangle,closetIntersection.intersectedTriangle.vertices[2]);
            //重心坐标就是做线性插值
            float alpha = (-(point.x-v1.x)*(v2.y-v1.y)+(point.y-v1.y)*(v2.x-v1.x))/(-(v0.x-v1.x)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.x-v1.x));
            float beta = (-(point.x-v2.x)*(v0.y-v2.y)+(point.y-v2.y)*(v0.x-v2.x))/(-(v1.x-v2.x)*(v0.y-v2.y)+(v1.y-v2.y)*(v0.x-v2.x));
            float gamma = 1-alpha-beta;

            float lightAngle0 = max(0.0f,glm::dot(n0, lightDirection0));
            float lightAngle1 = max(0.0f,glm::dot(n1, lightDirection1));
            float lightAngle2 = max(0.0f,glm::dot(n2, lightDirection2));

            glm::vec3 view0 = glm::normalize(cameraPosition - v0);
            glm::vec3 view1 = glm::normalize(cameraPosition - v1);
            glm::vec3 view2 = glm::normalize(cameraPosition - v2);

            glm::vec3 reflection0 = lightDirection0 - 2.0f*(n0 * glm::dot(lightDirection0, n0));
            glm::vec3 reflection1 = lightDirection1 - 2.0f*(n1 * glm::dot(lightDirection1, n1));
            glm::vec3 reflection2 = lightDirection2 - 2.0f*(n2 * glm::dot(lightDirection2, n2));

            float specularLight0 = glm::pow(glm::dot(view0,reflection0),specularExponent);
            float specularLight1 = glm::pow(glm::dot(view1,reflection1),specularExponent);
            float specularLight2 = glm::pow(glm::dot(view2,reflection2),specularExponent);

            float lightIntensity0 = PL0*lightAngle0+ specularLight0+ambientStrength;
            float lightIntensity1 = PL1*lightAngle1+ specularLight1+ambientStrength;
            float lightIntensity2 = PL2*lightAngle2+ specularLight2+ambientStrength;

            float lightIntensity = lightIntensity0*alpha + lightIntensity1*beta + lightIntensity2*gamma;

            //cout<<lightDistance<<endl;
            float red = std::min(255.0f,(255)*lightIntensity);
            float green = std::min(255.0f,(0)*lightIntensity);
            float blue = std::min(255.0f,(0)*lightIntensity);
            //判断点是否在圆内
            if(closetIntersection.intersectionPoint == glm::vec3(0, 0, 0)) {
                red = 0;
                green = 0;
                blue = 0;
            }
            uint32_t color = (255 << 24) + (int(round(red)) << 16) + (int(round(green)) << 8) + int(round(blue));
            window.setPixelColour(x,y,color);
        }
    }
}

void drawGouraudBox(DrawingWindow &window, float FL, float S, glm::vec3 cameraPosition, const std::vector<ModelTriangle>& triangle, glm::vec3 LP,glm::mat3x3 orientation){
    float sourceStrength = 8;
    float specularExponent = 256;
    for(int y=0; y<HEIGHT; y++){
        for(int x=0; x<WIDTH; x++){
            glm::vec3 v0;
            glm::vec3 v1;
            glm::vec3 v2;

            glm::vec3 worldPoint = getSphereWorldPoint(CanvasPoint(float(x),float (y)), FL, S, cameraPosition, orientation);
            glm::vec3 rayDirection = glm::normalize(worldPoint-cameraPosition);
            RayTriangleIntersection closetIntersection = getClosetsIntersection(cameraPosition, rayDirection, triangle);
            if(closetIntersection.intersectionPoint == glm::vec3{0,0,0}){
                uint32_t color = (255 << 24) + (0 << 16) + (0 << 8) + 0;
                window.setPixelColour(x,y,color);
                continue;
            }

            v0 = closetIntersection.intersectedTriangle.vertices[0];
            v1 = closetIntersection.intersectedTriangle.vertices[1];
            v2 = closetIntersection.intersectedTriangle.vertices[2];

            glm::vec3 point = closetIntersection.intersectionPoint;
            float alpha;
            float beta;
            float gamma;

            if(((v0.x-v1.x)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.x-v1.x))==0){
                alpha = (-(point.z-v1.z)*(v2.y-v1.y)+(point.y-v1.y)*(v2.z-v1.z))/(-(v0.z-v1.z)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.z-v1.z));
                beta = (-(point.z-v2.z)*(v0.y-v2.y)+(point.y-v2.y)*(v0.z-v2.z))/(-(v1.z-v2.z)*(v0.y-v2.y)+(v1.y-v2.y)*(v0.z-v2.z));
                gamma = 1-alpha-beta;
            }else if(v1.y==v2.y||v0.y==v1.y){
                alpha = (-(point.x-v1.x)*(v2.z-v1.z)+(point.z-v1.z)*(v2.x-v1.x))/(-(v0.x-v1.x)*(v2.z-v1.z)+(v0.z-v1.z)*(v2.x-v1.x));
                beta = (-(point.x-v2.x)*(v0.z-v2.z)+(point.z-v2.z)*(v0.x-v2.x))/(-(v1.x-v2.x)*(v0.z-v2.z)+(v1.z-v2.z)*(v0.x-v2.x));
                gamma = 1-alpha-beta;
            }else{
                alpha = (-(point.x-v1.x)*(v2.y-v1.y)+(point.y-v1.y)*(v2.x-v1.x))/(-(v0.x-v1.x)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.x-v1.x));
                beta = (-(point.x-v2.x)*(v0.y-v2.y)+(point.y-v2.y)*(v0.x-v2.x))/(-(v1.x-v2.x)*(v0.y-v2.y)+(v1.y-v2.y)*(v0.x-v2.x));
                gamma = 1-alpha-beta;
            }
            glm::vec3 lightDirection0 = glm::normalize(v0 - LP);
            glm::vec3 lightDirection1 = glm::normalize(v1 - LP);
            glm::vec3 lightDirection2 = glm::normalize(v2 - LP);

            glm::vec3 lightDirection = glm::normalize(closetIntersection.intersectionPoint - LP);


            float lightDistance0 = glm::length(LP - v0);
            float lightDistance1 = glm::length(LP - v1);
            float lightDistance2 = glm::length(LP - v2);

            RayTriangleIntersection light = getClosetsIntersection(LP, lightDirection, triangle);

            float PL0 = sourceStrength / float(( 4 *PI * lightDistance0 * lightDistance0));
            float PL1 = sourceStrength / float(( 4 *PI * lightDistance1 * lightDistance1));
            float PL2 = sourceStrength / float(( 4 *PI * lightDistance2 * lightDistance2));

            glm::vec3 n0 = vertexNormal(triangle, closetIntersection.intersectedTriangle.vertices[0]);
            glm::vec3 n1 = vertexNormal(triangle,closetIntersection.intersectedTriangle.vertices[1]);
            glm::vec3 n2 = vertexNormal(triangle,closetIntersection.intersectedTriangle.vertices[2]);

            float lightAngle0 = max(0.0f,glm::dot(n0, lightDirection0));
            float lightAngle1 = max(0.0f,glm::dot(n1, lightDirection1));
            float lightAngle2 = max(0.0f,glm::dot(n2, lightDirection2));

            glm::vec3 view0 = glm::normalize(cameraPosition - v0);
            glm::vec3 view1 = glm::normalize(cameraPosition - v1);
            glm::vec3 view2 = glm::normalize(cameraPosition - v2);

            glm::vec3 reflection0 = lightDirection0 - 2.0f*(n0 * glm::dot(lightDirection0, n0));
            glm::vec3 reflection1 = lightDirection1 - 2.0f*(n1 * glm::dot(lightDirection1, n1));
            glm::vec3 reflection2 = lightDirection2 - 2.0f*(n2 * glm::dot(lightDirection2, n2));

            float specularLight0 = glm::pow(glm::dot(view0,reflection0),specularExponent);
            float specularLight1 = glm::pow(glm::dot(view1,reflection1),specularExponent);
            float specularLight2 = glm::pow(glm::dot(view2,reflection2),specularExponent);

            float lightIntensity0 = PL0*lightAngle0+ specularLight0+ambientStrength;
            float lightIntensity1 = PL1*lightAngle1+ specularLight1+ambientStrength;
            float lightIntensity2 = PL2*lightAngle2+ specularLight2+ambientStrength;

            float lightIntensity = lightIntensity0*alpha + lightIntensity1*beta + lightIntensity2*gamma;
            //cout<<lightDistance<<endl;
            float red = std::max(0.0f,std::min(255.0f, float ((closetIntersection.intersectedTriangle.colour.red)) * lightIntensity));
            float green = std::max(0.0f,std::min(255.0f, float ((closetIntersection.intersectedTriangle.colour.green)) * lightIntensity));
            float blue = std::max(0.0f,std::min(255.0f, float ((closetIntersection.intersectedTriangle.colour.blue)) * lightIntensity));

            if(closetIntersection.triangleIndex != light.triangleIndex){
                red = red*ambientStrength;
                green = green*ambientStrength;
                blue = blue*ambientStrength;
            }
            uint32_t color = (255 << 24) + (int(round(red)) << 16) + (int(round(green)) << 8) + int(round(blue));
            window.setPixelColour(x,y,color);
        }
    }
}


void drawPhoneSphere(DrawingWindow &window, float FL, float S, glm::vec3 cameraPosition, const std::vector<ModelTriangle>& triangle, glm::vec3 LP,glm::mat3 orientation){
    float sourceStrength = 8;
    float specularExponent = 256;
    for(int y=0; y<HEIGHT; y++){
        for(int x=0; x<WIDTH; x++){
            glm::vec3 v0;
            glm::vec3 v1;
            glm::vec3 v2;

            glm::vec3 worldPoint = getSphereWorldPoint(CanvasPoint(float(x),float(y)), FL, S, cameraPosition, orientation);
            glm::vec3 rayDirection = glm::normalize(worldPoint-cameraPosition);
            RayTriangleIntersection closetIntersection = getClosetsIntersection(cameraPosition, rayDirection, triangle);

            v0 = closetIntersection.intersectedTriangle.vertices[0];
            v1 = closetIntersection.intersectedTriangle.vertices[1];
            v2 = closetIntersection.intersectedTriangle.vertices[2];

            glm::vec3 lightDirection = glm::normalize(closetIntersection.intersectionPoint - LP);
            RayTriangleIntersection light = getClosetsIntersection(LP, lightDirection, triangle);
            float lightDistance = glm::length(LP - closetIntersection.intersectionPoint);
            float PL = sourceStrength / float(( 4 *PI * lightDistance * lightDistance));

            glm::vec3 point = closetIntersection.intersectionPoint;
            glm::vec3 n0 = vertexNormal(triangle, closetIntersection.intersectedTriangle.vertices[0]);
            glm::vec3 n1 = vertexNormal(triangle,closetIntersection.intersectedTriangle.vertices[1]);
            glm::vec3 n2 = vertexNormal(triangle,closetIntersection.intersectedTriangle.vertices[2]);

            float alpha = (-(point.x-v1.x)*(v2.y-v1.y)+(point.y-v1.y)*(v2.x-v1.x))/(-(v0.x-v1.x)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.x-v1.x));
            float beta = (-(point.x-v2.x)*(v0.y-v2.y)+(point.y-v2.y)*(v0.x-v2.x))/(-(v1.x-v2.x)*(v0.y-v2.y)+(v1.y-v2.y)*(v0.x-v2.x));
            float gamma = 1-alpha-beta;

            glm::vec3 normal = n0*alpha + n1*beta + gamma*n2;

            float lightAngle = max(0.0f,glm::dot(normal, lightDirection));
            glm::vec3 view = glm::normalize(cameraPosition - closetIntersection.intersectionPoint);
            glm::vec3 reflection = lightDirection - 2.0f*(normal * glm::dot(lightDirection, normal));
            float specularLight = glm::pow(glm::dot(view,reflection),specularExponent);
            float brightness = PL * lightAngle + specularLight + ambientStrength;

            float red = std::min(255.0f, (255) * brightness);
            float green = std::min(255.0f, (0) * brightness);
            float blue = std::min(255.0f, (0) * brightness);
            //判断点是否在圆内
            if(closetIntersection.intersectionPoint == glm::vec3(0, 0, 0)) {
                red = 0;
                green = 0;
                blue = 0;
            }
            uint32_t color = (255 << 24) + (int(round(red)) << 16) + (int(round(green)) << 8) + int(round(blue));
            window.setPixelColour(x,y,color);
        }
    }
}

void drawPhoneBox(DrawingWindow &window, float FL, float S, glm::vec3 cameraPosition, const std::vector<ModelTriangle>& triangle, glm::vec3 LP,glm::mat3 orientation){
    float sourceStrength = 8;
    float specularExponent = 256;
    for(int y=0; y<HEIGHT; y++){
        for(int x=0; x<WIDTH; x++){
            glm::vec3 v0;
            glm::vec3 v1;
            glm::vec3 v2;

            glm::vec3 worldPoint = getSphereWorldPoint(CanvasPoint(float(x),float(y)), FL, S, cameraPosition, orientation);
            glm::vec3 rayDirection = glm::normalize(worldPoint-cameraPosition);
            RayTriangleIntersection closetIntersection = getClosetsIntersection(cameraPosition, rayDirection, triangle);

            if(closetIntersection.intersectionPoint == glm::vec3{0,0,0}){
                uint32_t color = (255 << 24) + (0 << 16) + (0 << 8) + 0;
                window.setPixelColour(x,y,color);
                continue;
            }

            v0 = closetIntersection.intersectedTriangle.vertices[0];
            v1 = closetIntersection.intersectedTriangle.vertices[1];
            v2 = closetIntersection.intersectedTriangle.vertices[2];

            glm::vec3 lightDirection = glm::normalize(closetIntersection.intersectionPoint - LP);
            RayTriangleIntersection light = getClosetsIntersection(LP, lightDirection, triangle);
            float lightDistance = glm::length(LP - closetIntersection.intersectionPoint);
            float PL = sourceStrength / float(( 4 *PI * lightDistance * lightDistance));

            glm::vec3 point = closetIntersection.intersectionPoint;
            glm::vec3 n0 = vertexNormal(triangle, closetIntersection.intersectedTriangle.vertices[0]);
            glm::vec3 n1 = vertexNormal(triangle,closetIntersection.intersectedTriangle.vertices[1]);
            glm::vec3 n2 = vertexNormal(triangle,closetIntersection.intersectedTriangle.vertices[2]);

            float alpha;
            float beta;
            float gamma;

            if(((v0.x-v1.x)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.x-v1.x))==0){
                alpha = (-(point.z-v1.z)*(v2.y-v1.y)+(point.y-v1.y)*(v2.z-v1.z))/(-(v0.z-v1.z)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.z-v1.z));
                beta = (-(point.z-v2.z)*(v0.y-v2.y)+(point.y-v2.y)*(v0.z-v2.z))/(-(v1.z-v2.z)*(v0.y-v2.y)+(v1.y-v2.y)*(v0.z-v2.z));
                gamma = 1-alpha-beta;
            }else{
                alpha = (-(point.x-v1.x)*(v2.y-v1.y)+(point.y-v1.y)*(v2.x-v1.x))/(-(v0.x-v1.x)*(v2.y-v1.y)+(v0.y-v1.y)*(v2.x-v1.x));
                beta = (-(point.x-v2.x)*(v0.y-v2.y)+(point.y-v2.y)*(v0.x-v2.x))/(-(v1.x-v2.x)*(v0.y-v2.y)+(v1.y-v2.y)*(v0.x-v2.x));
                gamma = 1-alpha-beta;
            }
            glm::vec3 normal = n0*alpha + n1*beta + gamma*n2;
            float lightAngle = max(0.0f,glm::dot(normal, lightDirection));
            glm::vec3 view = glm::normalize(cameraPosition - closetIntersection.intersectionPoint);
            glm::vec3 reflection = lightDirection - 2.0f*(normal * glm::dot(lightDirection, normal));
            float specularLight = glm::pow(glm::dot(view,reflection),specularExponent);

            float brightness = PL * lightAngle + specularLight + ambientStrength;

            float red = std::max(0.0f,std::min(255.0f, float ((closetIntersection.intersectedTriangle.colour.red)) * brightness));
            float green = std::max(0.0f,std::min(255.0f, float ((closetIntersection.intersectedTriangle.colour.green)) * brightness));
            float blue = std::max(0.0f,std::min(255.0f, float ((closetIntersection.intersectedTriangle.colour.blue)) * brightness));

            if(closetIntersection.triangleIndex != light.triangleIndex){
                red = red*ambientStrength;
                green = green*ambientStrength;
                blue = blue*ambientStrength;
            }
            uint32_t color = (255 << 24) + (int(round(red)) << 16) + (int(round(green)) << 8) + int(round(blue));
            window.setPixelColour(x,y,color);
        }
    }
}

void get_xyz(double x1,double y1,double z1,
             double x2,double y2,double z2,
             double x3,double y3,double z3,
             double x4,double y4,double z4)//空间四点确定球心坐标(克莱姆法则)
{
    double x, y, z;
    double a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,d,d1,d2,d3;
    a11=2*(x2-x1); a12=2*(y2-y1); a13=2*(z2-z1);
    a21=2*(x3-x2); a22=2*(y3-y2); a23=2*(z3-z2);
    a31=2*(x4-x3); a32=2*(y4-y3); a33=2*(z4-z3);
    b1=x2*x2-x1*x1+y2*y2-y1*y1+z2*z2-z1*z1;
    b2=x3*x3-x2*x2+y3*y3-y2*y2+z3*z3-z2*z2;
    b3=x4*x4-x3*x3+y4*y4-y3*y3+z4*z4-z3*z3;
    d=a11*a22*a33+a12*a23*a31+a13*a21*a32-a11*a23*a32-a12*a21*a33-a13*a22*a31;
    d1=b1*a22*a33+a12*a23*b3+a13*b2*a32-b1*a23*a32-a12*b2*a33-a13*a22*b3;
    d2=a11*b2*a33+b1*a23*a31+a13*a21*b3-a11*a23*b3-b1*a21*a33-a13*b2*a31;
    d3=a11*a22*b3+a12*b2*a31+b1*a21*a32-a11*b2*a32-a12*a21*b3-b1*a22*a31;
    x=d1/d;
    y=d2/d;
    z=d3/d;
    sphereCenterPoint = glm::vec3(x, y, z);
}

std::vector<ModelTriangle> loadNewOBJFile(const std::vector<string>& objFilenames, float SF, bool t){
    std::vector<ModelTriangle> triangle;
    std::vector<glm::vec3> allSpherePoint;
    for(const string& objFilename : objFilenames){
        string fileText;
        std::vector<glm::vec3> points;
        std::vector<TexturePoint> vt;
        std::string colourName;
        std::map<string,Colour> palette;
        if(objFilename == "textured-cornell-box.obj"){
            palette = readMaterialFile("textured-cornell-box.mtl");
        }

        ifstream MyReadFile(objFilename);
        while (getline(MyReadFile, fileText)){
            std::vector<std::string> splitDelimiter = split(fileText,' ');
            if(splitDelimiter[0]=="v"){
                float x = std::stof(splitDelimiter[1]);
                float y = std::stof(splitDelimiter[2]);
                float z = std::stof(splitDelimiter[3]);
                glm::vec3 p = {x * SF, y * SF, z * SF};
                //cout<<p.x<<endl;
                points.push_back(p);
            }
            if(splitDelimiter[0]=="f"){
                //std::cout<<s[3]<<std::endl;
                if(objFilename == "textured-cornell-box.obj") {
                    int fx = std::stoi(splitDelimiter[1]);
                    int fy = std::stoi(splitDelimiter[2]);
                    int fz = std::stoi(splitDelimiter[3]);
                    std::array<glm::vec3,3> face = {points[fx-1],points[fy-1],points[fz-1]};
                    ModelTriangle m = ModelTriangle();
                    m.vertices = face;
                    m.colour = palette[colourName];
                    m.colour.name = colourName;
                    if (m.colour.name == "Cobbles") {
                        std::vector<std::string> f0 = split(splitDelimiter[1], '/');
                        std::vector<std::string> f1 = split(splitDelimiter[2], '/');
                        std::vector<std::string> f2 = split(splitDelimiter[3], '/');
                        m.texturePoints = {vt[(std::stoi(f0[1])) - 1], vt[(std::stoi(f1[1])) - 1],
                                           vt[(std::stoi(f2[1])) - 1]};
                    }
                    //cout<<m<<endl;
                    m.normal=glm::normalize(glm::cross(m.vertices[2]-m.vertices[0],m.vertices[1]-m.vertices[0]));
                    // cout<<m.colour.name<<endl;
                    triangle.push_back(m);
                }else if(objFilename == "sphere.obj"){
                    glm::vec3 move (0.4,1.1,-0.2);
                    //glm::vec3 move (99,99,99);
                    int fx = std::stoi(splitDelimiter[1]);
                    int fy = std::stoi(splitDelimiter[2]);
                    int fz = std::stoi(splitDelimiter[3]);
                    std::array<glm::vec3,3> face = {points[fx-1]-move,points[fy-1]-move,points[fz-1]-move};
                    ModelTriangle m = ModelTriangle();
                    m.vertices = face;
                    m.colour = Colour(135,206,235);
                    m.colour.name = "SkyBlue";
                    allSpherePoint.push_back(m.vertices[0]);
                    allSpherePoint.push_back(m.vertices[1]);
                    allSpherePoint.push_back(m.vertices[2]);
                    //cout<<m<<endl;
                    m.normal=glm::normalize(glm::cross(m.vertices[2]-m.vertices[0],m.vertices[1]-m.vertices[0]));
                    // cout<<m.colour.name<<endl;
                    triangle.push_back(m);
                }else if(objFilename == "logo.obj"){
                    glm::vec3 move (-0.2234,-0.1987,0.9957);
                    int fx = std::stoi(splitDelimiter[1]);
                    int fy = std::stoi(splitDelimiter[2]);
                    int fz = std::stoi(splitDelimiter[3]);
                    std::array<glm::vec3,3> face = {points[fx-1]/250.0f-move,points[fy-1]/250.0f-move,points[fz-1]/250.0f-move};
                    ModelTriangle m = ModelTriangle();
                    m.vertices = face;
                    m.colour = Colour(237,145,33);
                    m.colour.name = "Carrot";
                    std::vector<std::string> f0 = split(splitDelimiter[1], '/');
                    std::vector<std::string> f1 = split(splitDelimiter[2], '/');
                    std::vector<std::string> f2 = split(splitDelimiter[3], '/');
                    m.texturePoints = {vt[(std::stoi(f0[1])) - 1], vt[(std::stoi(f1[1])) - 1], vt[(std::stoi(f2[1])) - 1]};
                    //cout<<m<<endl;
                    m.normal=glm::normalize(glm::cross(m.vertices[2]-m.vertices[0],m.vertices[1]-m.vertices[0]));
                    // cout<<m.colour.name<<endl;
                    triangle.push_back(m);
                }

            }
            if(splitDelimiter[0]=="vt"){
                float x = std::stof(splitDelimiter[1]);
                float y = std::stof(splitDelimiter[2]);
                TexturePoint p = {x, y};
                vt.push_back(p);
            }
            if(splitDelimiter[0]=="usemtl"){
                //std::cout<<splitDelimiter[1]<<std::endl;
                colourName = splitDelimiter[1];
            }
        }
        MyReadFile.close();
    }
    get_xyz(allSpherePoint[0].x, allSpherePoint[0].y, allSpherePoint[0].z,
    allSpherePoint[10].x, allSpherePoint[10].y, allSpherePoint[10].z,
    allSpherePoint[20].x, allSpherePoint[20].y, allSpherePoint[20].z,
    allSpherePoint[30].x, allSpherePoint[30].y, allSpherePoint[30].z);
    sphereRadius = glm::distance(sphereCenterPoint, allSpherePoint[0]);
    /* std::cout << glm::to_string(sphereCenterPoint) << std::endl;
    std::cout << sphereRadius << std::endl;
    exit(1); */
    return triangle;
}

std::vector<ModelTriangle> loadOBJFile(const string& objFilename, float SF, bool t){
	string fileText;
	std::vector<ModelTriangle> triangle;
	std::vector<glm::vec3> points;
    std::vector<TexturePoint> vt;
    std::string colourName;
	std::map<string,Colour> palette;
	if(t) palette = readMaterialFile("textured-cornell-box.mtl");
	ifstream MyReadFile(objFilename);
	while (getline(MyReadFile, fileText)){
		std::vector<std::string> splitDelimiter = split(fileText,' ');
		if(splitDelimiter[0]=="v"){
			float x = std::stof(splitDelimiter[1]);
			float y = std::stof(splitDelimiter[2]);
			float z = std::stof(splitDelimiter[3]);
			glm::vec3 p = {x * SF, y * SF, z * SF};
			//cout<<p.x<<endl;
			points.push_back(p);
	 	}
	 	if(splitDelimiter[0]=="f"){
			//std::cout<<s[3]<<std::endl;
            int fx = std::stoi(splitDelimiter[1]);
            int fy = std::stoi(splitDelimiter[2]);
            int fz = std::stoi(splitDelimiter[3]);
	 		std::array<glm::vec3,3> face = {points[fx-1],points[fy-1],points[fz-1]};
	 		ModelTriangle m = ModelTriangle();
	 		m.vertices = face;
			if(t){
                m.colour = palette[colourName];
                m.colour.name = colourName;
                if(m.colour.name == "Cobbles"){
                    std::vector<std::string> f0 = split(splitDelimiter[1],'/');
                    std::vector<std::string> f1 = split(splitDelimiter[2],'/');
                    std::vector<std::string> f2 = split(splitDelimiter[3],'/');
                    m.texturePoints={vt[(std::stoi(f0[1]))-1],vt[(std::stoi(f1[1]))-1],vt[(std::stoi(f2[1]))-1]};
                }
            }
			//cout<<m<<endl;
			m.normal=glm::normalize(glm::cross(m.vertices[2]-m.vertices[0],m.vertices[1]-m.vertices[0]));
           // cout<<m.colour.name<<endl;
	 		triangle.push_back(m);
	 	}
        if(splitDelimiter[0]=="vt"){
            float x = std::stof(splitDelimiter[1]);
            float y = std::stof(splitDelimiter[2]);
            TexturePoint p = {x, y};
            vt.push_back(p);
        }
		if(splitDelimiter[0]=="usemtl"){
			//std::cout<<splitDelimiter[1]<<std::endl;
			colourName = splitDelimiter[1];
		}
	}
	MyReadFile.close();
	return triangle;
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		float move = 0.1;
        zBuffer();
        bool rotation = false;
        bool transform = true;
        std::vector<std::string> objFileNames {"textured-cornell-box.obj", "sphere.obj", "logo.obj"};
		std::vector<ModelTriangle> l= loadNewOBJFile(objFileNames,scalingFactor,true);
		std::vector<ModelTriangle> sphere= loadOBJFile("sphere.obj",scalingFactor,false);
		if (event.key.keysym.sym == SDLK_LEFT) {
			std::cout << "light" << std::endl;
			window.clearPixels();
            lightPosition.z = lightPosition.z + move;
			//rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
            ray(window, focalLength, scalar, boxCameraPosition, l, lightPosition, cameraOrientation);
        }else if (event.key.keysym.sym == SDLK_RIGHT){
			std::cout << "RIGHT" << std::endl;
			window.clearPixels();
            boxCameraPosition.x = boxCameraPosition.x - move;
			//rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
            ray(window, focalLength, scalar, boxCameraPosition, l, lightPosition, cameraOrientation);
        }else if (event.key.keysym.sym == SDLK_UP){
			std::cout << "UP" << std::endl;
			window.clearPixels();
            boxCameraPosition.y = boxCameraPosition.y - move;
			//rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
            ray(window, focalLength, scalar, boxCameraPosition, l, lightPosition, cameraOrientation);
        }else if (event.key.keysym.sym == SDLK_DOWN) {
			window.clearPixels();
			std::cout << "DOWN" << std::endl;
            boxCameraPosition.y = boxCameraPosition.y + move;
			//rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
            ray(window, focalLength, scalar, boxCameraPosition, l, lightPosition, cameraOrientation);

        }else if (event.key.keysym.sym == SDLK_a){
			std::cout << "FORWARD" << std::endl;
			window.clearPixels();
            boxCameraPosition.z = boxCameraPosition.z - move;
			rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
		}else if (event.key.keysym.sym == SDLK_b){
			std::cout << "BACKWARD" << std::endl;
			window.clearPixels();
            boxCameraPosition.z = boxCameraPosition.z + move;
			rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
		}else if (event.key.keysym.sym == SDLK_x){
			std::cout << "ROTATION OF X" << std::endl;
			window.clearPixels();
			glm::mat3x3 x = {{1,0,0},{0,cos(0.1),-sin(0.1)},{0,sin(0.1),cos(0.1)}};
            boxCameraPosition = x * boxCameraPosition;
            cameraOrientation = x * cameraOrientation;
			//rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
            ray(window, focalLength, rayTraceScalar, boxCameraPosition, l, lightPosition, cameraOrientation);
        }else if (event.key.keysym.sym == SDLK_y) {
			std::cout << "ROTATION OF Y" << std::endl;
			window.clearPixels();
			glm::mat3x3 x = {{cos(0.1),0,sin(0.1)},{0,1,0},{-sin(0.1),0,cos(0.1)}};
            boxCameraPosition = x * boxCameraPosition;
            cameraOrientation = x * cameraOrientation;
            cout<<boxCameraPosition.z<<endl;
			//rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
            ray(window, focalLength, rayTraceScalar, boxCameraPosition, l, lightPosition, cameraOrientation);
        }else if(event.key.keysym.sym == SDLK_j){
			std::cout << "RAY TRACE" << std::endl;
			window.clearPixels();
			ray(window, focalLength, rayTraceScalar, boxCameraPosition, l, lightPosition, cameraOrientation);
		}else if(event.key.keysym.sym == SDLK_9){
            std::cout << "RAY TRACE" << std::endl;
            window.clearPixels();
            boxCameraPosition.x = boxCameraPosition.x + move;
            ray(window, focalLength, rayTraceScalar, boxCameraPosition, l, lightPosition, cameraOrientation);
        }else if(event.key.keysym.sym == SDLK_0){
            std::cout << "RAY TRACE" << std::endl;
            window.clearPixels();
            boxCameraPosition.x = boxCameraPosition.x - move;
            ray(window, focalLength, rayTraceScalar, boxCameraPosition, l, lightPosition, cameraOrientation);
        }else if(event.key.keysym.sym == SDLK_h){
			std::cout << "AMBIENT LIGHT ADD" << std::endl;
			window.clearPixels();
			if(ambientStrength>=1) ambientStrength = 0.1;
			ambientStrength += move;
			ray(window, focalLength, scalar, boxCameraPosition, l, lightPosition, cameraOrientation);
		}else if(event.key.keysym.sym == SDLK_1){
			std::cout << "draw Flat sphere" << std::endl;
			window.clearPixels();
            drawFlatSphere(window, focalLength, scalar, sphereCameraPosition, sphere, sphereLightPosition, cameraOrientation);
		}else if (event.key.keysym.sym == SDLK_2){
            std::cout << "draw Goraud sphere" << std::endl;
            window.clearPixels();
            drawGouraudSphere(window, focalLength, scalar, sphereCameraPosition, sphere, sphereLightPosition, cameraOrientation);
        }else if (event.key.keysym.sym == SDLK_3){
            std::cout << "draw Phone sphere" << std::endl;
            window.clearPixels();
            drawPhoneSphere(window, focalLength, scalar, sphereCameraPosition, sphere, sphereLightPosition, cameraOrientation);
        }else if (event.key.keysym.sym == SDLK_4){
            std::cout << "Goraud sphere light add" << std::endl;
            window.clearPixels();
            sphereLightPosition.z = sphereLightPosition.z-move;
            drawGouraudSphere(window, focalLength, scalar, sphereCameraPosition, sphere, sphereLightPosition, cameraOrientation);
        }else if(event.key.keysym.sym == SDLK_5){
            std::cout << "draw ground box" << std::endl;
            window.clearPixels();
            if(ambientStrength>=1) ambientStrength = 0.1;
            ambientStrength += move;
            drawGouraudBox(window, focalLength, rayTraceScalar, boxCameraPosition, l, lightPosition, cameraOrientation);
        }else if(event.key.keysym.sym == SDLK_6){
            std::cout << "draw phone box" << std::endl;
            window.clearPixels();
            if(ambientStrength>=1) ambientStrength = 0.1;
            ambientStrength += move;
            drawPhoneBox(window, focalLength, rayTraceScalar, boxCameraPosition, l, lightPosition, cameraOrientation);
        }else if(event.key.keysym.sym == SDLK_7){
            std::cout << "draw phone box light move" << std::endl;
            window.clearPixels();
            lightPosition.y = lightPosition.y+move;
            drawPhoneBox(window, focalLength, rayTraceScalar, boxCameraPosition, l, lightPosition, cameraOrientation);
        }else if(event.key.keysym.sym == SDLK_d){
			std::cout << "RAY TRACE LEFT" << std::endl;
			window.clearPixels();
            boxCameraPosition.x = boxCameraPosition.x-move;
			ray(window, focalLength, scalar, boxCameraPosition, l, lightPosition, cameraOrientation);
		}else if(event.key.keysym.sym == SDLK_f){
			std::cout << "RAY TRACE RIGHT" << std::endl;
			window.clearPixels();
			boxCameraPosition.x = boxCameraPosition.x+move;
			ray(window, focalLength, 50, boxCameraPosition, l, lightPosition, cameraOrientation);
		}else if(event.key.keysym.sym == SDLK_g){
			std::cout << "RAY TRACE UP" << std::endl;
			window.clearPixels();
			boxCameraPosition.y = boxCameraPosition.y+move;
			ray(window, focalLength, 50, boxCameraPosition, l, lightPosition, cameraOrientation);
		}else if(event.key.keysym.sym == SDLK_h){
			std::cout << "RAY TRACE DOWN" << std::endl;
			window.clearPixels();
            boxCameraPosition.y = boxCameraPosition.y-move;
			ray(window, focalLength, scalar, boxCameraPosition, l, lightPosition, cameraOrientation);
		}else if (event.key.keysym.sym == SDLK_v) {
			CanvasTriangle triangle = strokedTriangles();
			Colour colour = {rand()%256,rand()%256,rand()%256};
		    lineDraw(window,triangle.v0(),triangle.v1(),colour);
			lineDraw(window,triangle.v1(),triangle.v2(),colour);
			lineDraw(window,triangle.v2(),triangle.v0(),colour);
		}else if (event.key.keysym.sym == SDLK_l) {
			CanvasTriangle triangle = strokedTriangles();
			Colour colour = {rand()%256,rand()%256,rand()%256};
			//filledTriangle(window,triangle,colour);
            occlusion(window,triangle,colour);
			// white edge
			drawTriangle(window,triangle);
		}else if (event.key.keysym.sym == SDLK_t){
			window.clearPixels();
			CanvasTriangle canvas = {{160,10},{300,230},{10,150}};
			canvas.v0().texturePoint = {195,5};
			canvas.v1().texturePoint = {395,380};
			canvas.v2().texturePoint = {65,330};
			CanvasTriangle canvasT = arrangeTriangle(canvas);
			TextureMap textureFile = TextureMap("texture.ppm");
			glm::mat3x3 affineT = affine(canvasT);
			//affineTransformation(window, canvasT, CanvasPoint(300, 230), affineT, textureFile);
			textureMapper(window,canvasT,affineT,textureFile);
			//drawTriangle(window,canvas);
		}else if (event.key.keysym.sym == SDLK_p){
			window.clearPixels();
			pointCloud(window, boxCameraPosition, focalLength, l);
		}else if(event.key.keysym.sym == SDLK_w){
			window.clearPixels();
			wireFrameRender(window, boxCameraPosition, l, focalLength, cameraOrientation);
		}else if(event.key.keysym.sym == SDLK_r){
			window.clearPixels();
			rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
            //cout<<zBuffer[3][4]<<'/'<<zBuffer[5][1] <<endl;
        }else if (event.key.keysym.sym == SDLK_c) {
            zBuffer();
			window.clearPixels();
            boxCameraPosition = {0, 0, 4};
            sphereLightPosition = {-0.3, 1.3, 0.1};
            lightPosition = {0,0.5,0.5};
            sphereLightPosition = {-0.3, 1.3, 1.5};
            cameraOrientation ={{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
		}
	}else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);	
	auto height = float(window.height);
	auto width = float(window.width);
	SDL_Event event;
	bool press_X = false;
	bool press_Y = false;
	std::vector<ModelTriangle> l= loadOBJFile("textured-cornell-box.obj",scalingFactor,true);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		if(event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_m) press_X = true;
		if(event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_n) press_Y = true;
		if(press_X){
			cout << "X" << endl;
			if(event.key.keysym.sym == SDLK_s){
			std::cout << "STOP" << std::endl;
                zBuffer();
				window.clearPixels();
				press_X = false;
			}
            zBuffer();
			window.clearPixels();
			glm::mat3x3 x = {{1,0,0},{0,cos(0.01),-sin(0.01)},{0,sin(0.01),cos(0.01)}};;
            boxCameraPosition = x * boxCameraPosition;
            cameraOrientation = x * cameraOrientation;
			//cameraOrientation = lookAt(boxCameraPosition);
			//rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
            ray(window, focalLength, 50, boxCameraPosition, l, lightPosition, cameraOrientation);

        }else if(press_Y){
			cout << "Y" << endl;
			if(event.key.keysym.sym == SDLK_s){
			std::cout << "STOP" << std::endl;
                zBuffer();
				window.clearPixels();
				press_Y = false;
			}
            zBuffer();
			window.clearPixels();
			glm::mat3x3 x = glm::mat3(cos(0.01),0,-sin(0.01),0,1,0,sin(0.01),0,cos(0.01));
            boxCameraPosition = x * boxCameraPosition;
			//cameraOrientation = x*cameraOrientation;
			cameraOrientation = lookAt(boxCameraPosition);
            //rasterisedRender(window, l, boxCameraPosition, focalLength, cameraOrientation);
           ray(window, focalLength, 50, boxCameraPosition, l, lightPosition, cameraOrientation);

        }
		//rasterisedRender(window,l,boxCameraPosition,2.0,cameraOrientation);
		// Need to render thse frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
