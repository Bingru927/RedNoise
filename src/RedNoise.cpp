#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <map>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include "TextureMap.h"
#include "ModelTriangle.h"
#include <algorithm>
#include <string>
#include <iostream>
#include <iostream>
#include <cctype>
#include <cmath>
#include "RayTriangleIntersection.h"


using namespace std;


#define WIDTH 320
#define HEIGHT 240
#define PI 3.1415926


float buffer[HEIGHT][WIDTH];
void zbuffer(){
    for(int y = 0; y < HEIGHT; y++)
        for(int x = 0; x < WIDTH; x++)
            buffer[y][x] = INT32_MIN;
}

glm::vec3 cameraPosition (0.0, -0.0, 4.0);
glm::mat3 cameraOritation (1,0,0,0,1,0,0,0,1);
glm::vec3 lightPosition (0,0.5,0.5);

float scalar = 2*HEIGHT/3;
float scalingFactor = 0.35;
float focalLength = 2.0;
float ambientStrength = 0.1;

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

std::vector<float> interpolateSingleFloats(float from, float to, int N){
    vector<float> v;
	float next = from;
	float a = (to - from) / (N - 1);
	for (int i = 0; i < N; i++) {
        v.push_back(next);
        next = next + a;
        //cout<<next<<endl;
	}
	return v;
}

void greyScaleInterpolation(DrawingWindow &window) {
	std::vector<float> result = interpolateSingleFloats(255,0,window.width);
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
	leftMost = interpolateThreeElementValues(topLeft, bottomLeft, window.height);
	rightMost = interpolateThreeElementValues(topRight, bottomRight, window.height);
	for (size_t y = 0; y < window.height; y++) {
		eachRow = interpolateThreeElementValues(leftMost[y],rightMost[y],window.width);
		for (size_t x = 0; x < window.width; x++) {
			uint32_t colour = (255 << 24) + (int(eachRow[x][0]) << 16) + (int(eachRow[x][1]) << 8) + int(eachRow[x][2]);
			window.setPixelColour(x, y, colour);
		}
	}
}

void lineDraw(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour){
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
    float zDiff = to.depth - from.depth;
	float numberOfSpace = max(abs(xDiff),abs(yDiff));
	float xStepSize = xDiff/numberOfSpace;
	float yStepSize = yDiff/numberOfSpace;
    float depthStepSize = zDiff/numberOfSpace;
	float red = colour.red;
	float green = colour.green;
	float blue = colour.blue;
	uint32_t color = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	for (float i =0.0; i<numberOfSpace; i++ ){
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
        float z = from.depth + (depthStepSize*i);
		if(int(round(y))<HEIGHT && int(round(y))>=0 && int(round(x))<WIDTH && int(round(x))>=0){
			if(buffer[int(round(y))][int(round(x))] <= z){
				buffer[int(round(y))][int(round(x))] = z;
				window.setPixelColour(round(x), round(y), color);
			}	
		}
	}
}

//draws "stroked" (unfilled) triangles
CanvasTriangle strokedTriangles(){
	CanvasTriangle triangle;
	for (int i =0; i<3; i++){
		float random_x = rand()%320;
		float random_y = rand()%240;
		if (random_x == random_y){
			random_x = rand()%320;
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


void filledTriangle(DrawingWindow &window,CanvasTriangle triangle, Colour colour){
	CanvasTriangle arranged = arrangeTriangle(triangle);
	CanvasPoint top = arranged.v0();
	CanvasPoint bottom = arranged.v2();
	CanvasPoint mid = arranged.v1();	
	for (float y = top.y; y<=mid.y;y++){
	 	float x_left = top.x-((top.x-bottom.x)*(y-top.y)/(bottom.y-top.y));
	 	float x_right = top.x+((y-top.y)*(mid.x-top.x)/(mid.y-top.y));
         //get buffer
	 	lineDraw(window,CanvasPoint(x_left,y),CanvasPoint(x_right,y),colour);
	}
	for (float y=mid.y; y<bottom.y;y++) {
        //cout << ((y-mid.y)*(mid.x-bottom.x)/(bottom.y-mid.y)) << endl;
        float x_left = top.x - ((top.x - bottom.x) * (y - top.y) / (bottom.y - top.y));
        float x_right = mid.x - ((y - mid.y) * (mid.x - bottom.x) / (bottom.y - mid.y));
        lineDraw(window, CanvasPoint(x_left, y), CanvasPoint(x_right, y), colour);
        // 	for(int x=x_left; x<x_right; x++) window.setPixelColour(x,y,color);
    }
}
unsigned int setColour(Colour colour){
    float red = colour.red;
    float green = colour.green;
    float blue = colour.blue;
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
	uint32_t c = textureMap.pixels[((p.y)*(int(textureMap.width))+p.x)-1];
	window.setPixelColour(point.x,(point.y),c);
	//return c; 
}

void textureMapper(DrawingWindow &window, CanvasTriangle arranged,glm::mat3x3 affineT, TextureMap textureMap){
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
std::map<string,Colour> readMaterialFile(string MLTfilename){
 ifstream MyReadFile(MLTfilename);
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
		c.red = std::stof(splitDelimiter[1]) * 255;
		c.green = std::stof(splitDelimiter[2]) * 255;
		c.blue = std::stof(splitDelimiter[3]) * 255;
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

CanvasPoint getCanvasIntersectionPoint (glm::vec3 cameraPosition, glm::vec3 vertexPosition,float focalLength,glm::mat3x3 cameraOritation){
	CanvasPoint modelVertex;
	int imagePlaneScaling = HEIGHT *2/3;
	glm::vec3 c = getCameraPoint(vertexPosition,cameraPosition);
	c = c*cameraOritation;
	//cout<<cameraPosition.x<<endl;
	modelVertex.x = round((-imagePlaneScaling) * focalLength * (c.x/c.z) + WIDTH/2);
	modelVertex.y = round(imagePlaneScaling * focalLength * (c.y/c.z) + HEIGHT/2);
    modelVertex.depth =abs(1/c.z);
    //cout<<modelVertex.depth<<endl;
	return modelVertex;
}

void pointCloud(DrawingWindow &window,glm::vec3 cameraPosition,float focalLength,std::vector<ModelTriangle> triangle){
	for(int i=0;i<int(triangle.size());i++){
		for(int j=0;j<int((triangle[i].vertices).size());j++){
			glm::vec3 v= triangle[i].vertices[j];
			uint32_t color = (255 << 24) + (int(255) << 16) + (int(255) << 8) + int(255);
		    CanvasPoint p=getCanvasIntersectionPoint(cameraPosition,v,focalLength,cameraOritation);
			//std::cout<<p.x<<' '<<p.y<<std::endl;
			window.setPixelColour(p.x,p.y,color);
		}
 	}
}

void wireFrameRender(DrawingWindow &window,glm::vec3 cameraPosition,std::vector<ModelTriangle> triangle,float focalLength,glm::mat3x3 cameraOritation){
	for(ModelTriangle t : triangle){
		CanvasPoint v0=getCanvasIntersectionPoint(cameraPosition,t.vertices[0],focalLength,cameraOritation);
		CanvasPoint v1=getCanvasIntersectionPoint(cameraPosition,t.vertices[1],focalLength,cameraOritation);
		CanvasPoint v2=getCanvasIntersectionPoint(cameraPosition,t.vertices[2],focalLength,cameraOritation);
		CanvasTriangle ct = CanvasTriangle(v0,v1,v2);
		drawTriangle(window,ct);
	}
}

void setBufferTriangle(DrawingWindow &windows,CanvasPoint top,CanvasPoint v1,CanvasPoint v2, Colour colour){
    float d = abs(v1.y-top.y);
    CanvasPoint left=top;
    CanvasPoint right=top;
    int s=1;
    if(top.y - v2.y > 0) s=-1;
    for (float i =1; i<=d; i++){
        float radio = i/abs(d);
        left.x=top.x+((v1.x-top.x)*(radio));
        right.x = top.x+((v2.x-top.x)*(radio));
        left.y += s;
        right.y += s;
        left.depth = top.depth + ((v1.depth - top.depth)*radio);
        right.depth = top.depth + ((v2.depth - top.depth)*radio);
        lineDraw(windows,left,right,colour);
    }
}

void occlusion(DrawingWindow &window,CanvasTriangle triangle, Colour colour){
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

void rasterisedRender(DrawingWindow &window,std::vector<ModelTriangle> triangle,glm::vec3 cameraPosition,float focalLength,glm::mat3x3 cameraOritation){
	for(ModelTriangle t : triangle){
		CanvasPoint v0=getCanvasIntersectionPoint(cameraPosition,t.vertices[0],focalLength,cameraOritation);
		CanvasPoint v1=getCanvasIntersectionPoint(cameraPosition,t.vertices[1],focalLength,cameraOritation);
		CanvasPoint v2=getCanvasIntersectionPoint(cameraPosition,t.vertices[2],focalLength,cameraOritation);
		CanvasTriangle ct = CanvasTriangle(v0,v1,v2);
		//filledTriangle(window,ct,t.colour);
        occlusion(window,ct,t.colour);
    }
}
//void oribit(DrawingWindow &window,glm::mat3x3 x,std::vector<ModelTriangle> l,glm::vec3 cameraPosition,float focalLength,glm::mat3x3 cameraOritation){
//    cameraPosition = x*cameraPosition;
//    cameraOritation = x*cameraOritation;
//    rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
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
glm::vec3 getWorldPoint(CanvasPoint p,float focalLength,float scalar,glm::vec3 cameraPosition){
	glm::vec3 point;
	point.z = cameraPosition.z-focalLength;
	point.x =(point.z*(p.x-WIDTH/2))/(focalLength*scalar);
	point.y =(point.z*(p.y-HEIGHT/2))/(focalLength*(-scalar));
	return point;
}


RayTriangleIntersection getClosetsIntersection(glm::vec3 cameraPosition,glm::vec3 rayDirection,std::vector<ModelTriangle> triangle){
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
		if (t<tMin && t>=0 && u>=0 && u<=1 && v>=0 && v<=1 && (u + v) <= 1) {
			tMin = t;
			// Update RayTriangleIntersection to the nearest point
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

void shadow(DrawingWindow &window,float focalLength,float scalar,glm::vec3 cameraPosition,std::vector<ModelTriangle> triangle,glm::vec3 lightPosition){
	for(int y=0; y<HEIGHT; y++){
		for(int x=0; x<WIDTH; x++){
			glm::vec3 worldPoint = getWorldPoint(CanvasPoint(x,y),focalLength,scalar,cameraPosition);
			glm::vec3 rayDirection = worldPoint-cameraPosition;
			RayTriangleIntersection r = getClosetsIntersection(cameraPosition,rayDirection,triangle);
			glm::vec3 lightDirection = r.intersectionPoint-lightPosition;
			RayTriangleIntersection light = getClosetsIntersection(lightPosition,lightDirection,triangle);
			//To determine whether the point seen by the camera is the same as the point illuminated by the light, if it is the same it means that it is not in the shadow.
			if(r.triangleIndex==light.triangleIndex){
				uint32_t color = setColour(r.intersectedTriangle.colour);
				window.setPixelColour(x,y,color);
			}else window.setPixelColour(x,y,setColour(Colour(0,0,0)));
		}
	}
}
void ray(DrawingWindow &window,float focalLength,float scalar,glm::vec3 cameraPosition,std::vector<ModelTriangle> triangle,glm::vec3 lightPosition){
	float sourceStrength = 8;
	int specularExponent = 256;
	for(int y=0; y<HEIGHT; y++){
		for(int x=0; x<WIDTH; x++){
			glm::vec3 worldPoint = getWorldPoint(CanvasPoint(x,y),focalLength,scalar,cameraPosition);
			glm::vec3 rayDirection = glm::normalize(worldPoint-cameraPosition);
            //cout<<rayDirection.x<<rayDirection.y<<endl;
			RayTriangleIntersection r = getClosetsIntersection(cameraPosition,rayDirection,triangle);
			glm::vec3 lightDirection = glm::normalize(r.intersectionPoint-lightPosition);
			RayTriangleIntersection light = getClosetsIntersection(lightPosition,lightDirection,triangle);
			//model表面的点到灯的方向向量
			float lightDistance = glm::length(lightPosition-r.intersectionPoint);
			//cout<<lightDistance<<endl;
			float PL = sourceStrength / ( 4 *PI * lightDistance * lightDistance);
			float lightAngle = max(0.0f,glm::dot(r.intersectedTriangle.normal,lightDirection));
			glm::vec3 view = glm::normalize(cameraPosition-r.intersectionPoint);
			glm::vec3 reflection = lightDirection - 2.0f*(r.intersectedTriangle.normal*glm::dot(lightDirection,r.intersectedTriangle.normal));
			float specularLight = glm::pow(glm::dot(view,reflection), specularExponent);
			float lightIntensity = PL*lightAngle+ specularLight+ambientStrength;
			float red = std::max(0.0f,std::min(255.0f,(r.intersectedTriangle.colour.red)*lightIntensity));
			float green = std::max(0.0f,std::min(255.0f,(r.intersectedTriangle.colour.green)*lightIntensity));
			float blue = std::max(0.0f,std::min(255.0f,(r.intersectedTriangle.colour.blue)*lightIntensity));
			if(r.triangleIndex != light.triangleIndex){
				red = red*ambientStrength;
				green = green*ambientStrength;
				blue = blue*ambientStrength;
			}
			uint32_t color = (255 << 24) + (int(round(red)) << 16) + (int(round(green)) << 8) + int(round(blue));
			window.setPixelColour(x,y,color);
		}
	}
}



void drawSphere(DrawingWindow &window,float focalLength,float scalar,glm::vec3 cameraPosition,std::vector<ModelTriangle> triangle,glm::vec3 lightPosition){
	float sourceStrength = 8;
	int specularExponent = 258;
	for(int y=0; y<HEIGHT; y++){
		for(int x=0; x<WIDTH; x++){
			glm::vec3 worldPoint = getWorldPoint(CanvasPoint(x,y),focalLength,scalar,cameraPosition);
			glm::vec3 rayDirection = glm::normalize(worldPoint-cameraPosition);
			RayTriangleIntersection r = getClosetsIntersection(cameraPosition,rayDirection,triangle);
			glm::vec3 lightDirection = glm::normalize(r.intersectionPoint-lightPosition);
			RayTriangleIntersection light = getClosetsIntersection(lightPosition,lightDirection,triangle);
			//model表面的点到灯的方向向量
			float lightDistance = glm::length(lightPosition-r.intersectionPoint);
			//cout<<lightDistance<<endl;
			float PL = sourceStrength / ( 4 *PI * lightDistance * lightDistance);
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

std::vector<ModelTriangle> loadOBJFile( string objFilename, float scalingFactor,bool t){
	string fileText;
	std::vector<ModelTriangle> triangle;
	std::vector<glm::vec3> points;
	std::string colourName;
	std::map<string,Colour> palette;
	if(t) palette = readMaterialFile("cornell-box.mtl");
	ifstream MyReadFile(objFilename);
	while (getline(MyReadFile, fileText)){
		std::vector<std::string> splitDelimiter = split(fileText,' ');
		if(splitDelimiter[0]=="v"){
			float x = std::stof(splitDelimiter[1]);
			float y = std::stof(splitDelimiter[2]);
			float z = std::stof(splitDelimiter[3]);
			glm::vec3 p = {x*scalingFactor,y*scalingFactor,z*scalingFactor};
			//cout<<p.x<<endl;
			points.push_back(p);
	 	}
	 	if(splitDelimiter[0]=="f"){
			//std::vector<std::string> f = split(fileText,' ');
			//std::cout<<s[3]<<std::endl;
            int fx = std::stoi(splitDelimiter[1]);
            int fy = std::stoi(splitDelimiter[2]);
            int fz = std::stoi(splitDelimiter[3]);
	 		std::array<glm::vec3,3> face = {points[fx-1],points[fy-1],points[fz-1]};
	 		ModelTriangle m = ModelTriangle();
	 		m.vertices = face;
			if(t) m.colour = palette[colourName];
			//cout<<m<<endl;
			m.normal=glm::normalize(glm::cross(m.vertices[2]-m.vertices[0],m.vertices[1]-m.vertices[0]));
	 		triangle.push_back(m);
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
		zbuffer();
		std::vector<ModelTriangle> l= loadOBJFile("cornell-box.obj",scalingFactor,true);
		std::vector<ModelTriangle> sphere= loadOBJFile("sphere.obj",scalingFactor,false);
		if (event.key.keysym.sym == SDLK_LEFT) {
			std::cout << "LEFT" << std::endl;
			window.clearPixels();
			cameraPosition.x = cameraPosition.x + move;
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
		}else if (event.key.keysym.sym == SDLK_RIGHT){
			std::cout << "RIGHT" << std::endl;
			window.clearPixels();
			cameraPosition.x = cameraPosition.x - move;
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
		}else if (event.key.keysym.sym == SDLK_UP){
			std::cout << "UP" << std::endl;
			window.clearPixels();
			cameraPosition.y = cameraPosition.y - move;
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
		}else if (event.key.keysym.sym == SDLK_DOWN) {
			window.clearPixels();
			std::cout << "DOWN" << std::endl;
			cameraPosition.y = cameraPosition.y + move;
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
		}else if (event.key.keysym.sym == SDLK_a){
			std::cout << "FORWARD" << std::endl;
			window.clearPixels();
			cameraPosition.z = cameraPosition.z - move;
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
		}else if (event.key.keysym.sym == SDLK_b){
			std::cout << "BACKWARD" << std::endl;
			window.clearPixels();
			cameraPosition.z = cameraPosition.z + move;
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
		}else if (event.key.keysym.sym == SDLK_x){
			std::cout << "ROTATION OF X" << std::endl;
			window.clearPixels();
			glm::mat3x3 x = {{1,0,0},{0,cos(0.1),-sin(0.1)},{0,sin(0.1),cos(0.1)}};
			cameraPosition = x*cameraPosition;
			cameraOritation = x*cameraOritation;
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
		}else if (event.key.keysym.sym == SDLK_y) {
			std::cout << "ROTATION OF Y" << std::endl;
			window.clearPixels();
			glm::mat3x3 x = {{cos(0.1),0,sin(0.1)},{0,1,0},{-sin(0.1),0,cos(0.1)}};
			cameraPosition = x*cameraPosition;
			cameraOritation = x*cameraOritation;
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
		}else if(event.key.keysym.sym == SDLK_j){
			std::cout << "RAY TRACE" << std::endl;
			window.clearPixels();
			ray(window,focalLength,scalar,cameraPosition,l,lightPosition);
		}else if(event.key.keysym.sym == SDLK_4){
            std::cout << "LIGHT ADD" << std::endl;
            window.clearPixels();
            lightPosition.x = lightPosition.x-move;
            drawSphere(window,focalLength,scalar,cameraPosition,sphere,lightPosition);
        }else if(event.key.keysym.sym == SDLK_1){
			std::cout << "AMBIENT LIGHT ADD" << std::endl;
			window.clearPixels();	
			if(ambientStrength>=1) ambientStrength = 0.1;
			ambientStrength += move;
			ray(window,focalLength,scalar,cameraPosition,l,lightPosition);
		}else if(event.key.keysym.sym == SDLK_2){
			std::cout << "sphere" << std::endl;
			window.clearPixels();
			drawSphere(window,focalLength,scalar,cameraPosition,sphere,lightPosition);
		}else if (event.key.keysym.sym == SDLK_3){
            std::cout << "FORWARD" << std::endl;
            window.clearPixels();
           // cameraPosition.y = cameraPosition.y - move;
            cameraPosition.z = cameraPosition.z - move;
            ambientStrength += move;
            drawSphere(window,focalLength,scalar,cameraPosition,sphere,lightPosition);
        }else if (event.key.keysym.sym == SDLK_5){
            std::cout << "FORWARD" << std::endl;
            window.clearPixels();
           // cameraPosition.y = cameraPosition.y + move;
            cameraPosition.z = cameraPosition.z + move;
            ambientStrength += move;
            drawSphere(window,focalLength,scalar,cameraPosition,sphere,lightPosition);
        }else if(event.key.keysym.sym == SDLK_d){
			std::cout << "RAY TRACE LEFT" << std::endl;
			window.clearPixels();
			lightPosition.x = lightPosition.x-move;
			ray(window,focalLength,scalar,cameraPosition,l,lightPosition);
		}else if(event.key.keysym.sym == SDLK_f){
			std::cout << "RAY TRACE RIGHT" << std::endl;
			window.clearPixels();
			lightPosition.x = lightPosition.x+move;
			ray(window,focalLength,scalar,cameraPosition,l,lightPosition);
		}else if(event.key.keysym.sym == SDLK_g){
			std::cout << "RAY TRACE UP" << std::endl;
			window.clearPixels();
			lightPosition.y = lightPosition.y+move;
			ray(window,focalLength,scalar,cameraPosition,l,lightPosition);
		}else if(event.key.keysym.sym == SDLK_h){
			std::cout << "RAY TRACE DOWN" << std::endl;
			window.clearPixels();
			lightPosition.y = lightPosition.y-move;
			ray(window,focalLength,scalar,cameraPosition,l,lightPosition);
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
			TextureMap file = TextureMap("texture.ppm");
			glm::mat3x3 affineT = affine(canvasT);
			//affineTransformation(window, canvasT, CanvasPoint(300, 230), affineT, file);
			textureMapper(window,canvasT,affineT,file);
			//drawTriangle(window,canvas);
		}else if (event.key.keysym.sym == SDLK_p){
			window.clearPixels();
			pointCloud(window,cameraPosition,focalLength,l);
		}else if(event.key.keysym.sym == SDLK_w){
			window.clearPixels();
			wireFrameRender(window,cameraPosition,l,focalLength,cameraOritation);
		}else if(event.key.keysym.sym == SDLK_r){
			window.clearPixels();
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
            //cout<<zBuffer[3][4]<<'/'<<zBuffer[5][1] <<endl;
        }else if (event.key.keysym.sym == SDLK_c) {
            zbuffer();
			window.clearPixels();
			cameraPosition = {0,0,4};
			cameraOritation ={{1,0,0},{0,1,0},{0,0,1}};
		}
	}else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);	
	float height = window.height;
	float width = window.width;
	SDL_Event event;
	bool press_X = false;
	bool press_Y = false;
	std::vector<ModelTriangle> l= loadOBJFile("cornell-box.obj",scalingFactor,true);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		if(event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_m) press_X = true;
		if(event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_n) press_Y = true;
		if(press_X){
			cout << "X" << endl;
			if(event.key.keysym.sym == SDLK_s){
			std::cout << "STOP" << std::endl;
				zbuffer();
				window.clearPixels();
				press_X = false;
			}
			zbuffer();
			window.clearPixels();
			glm::mat3x3 x = {{1,0,0},{0,cos(0.01),-sin(0.01)},{0,sin(0.01),cos(0.01)}};;
			cameraPosition = x*cameraPosition;
			cameraOritation = x*cameraOritation;
			//cameraOritation = lookAt(cameraPosition);
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
		}else if(press_Y){
			cout << "Y" << endl;
			if(event.key.keysym.sym == SDLK_s){
			std::cout << "STOP" << std::endl;
				zbuffer();
				window.clearPixels();
				press_Y = false;
			}
			zbuffer();
			window.clearPixels();
			glm::mat3x3 x = glm::mat3(cos(0.01),0,-sin(0.01),0,1,0,sin(0.01),0,cos(0.01));
			cameraPosition = x*cameraPosition;
			//cameraOritation = x*cameraOritation;
			cameraOritation = lookAt(cameraPosition);
			rasterisedRender(window,l,cameraPosition,focalLength,cameraOritation);
		}
		//rasterisedRender(window,l,cameraPosition,2.0,cameraOritation);
		// Need to render thse frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
