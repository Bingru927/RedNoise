#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include "TextureMap.h"
#include <stdint.h>



#define WIDTH 320
#define HEIGHT 240
using namespace std;

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
	 }
	 return v;
}

void greyScaleInterpolation(DrawingWindow &window) {
	window.clearPixels();
	std::vector<float> result = interpolateSingleFloats(255,0,window.width);
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = result[x];
			float green = result[x];
			float blue = result[x];
			//i第一个是灰度 相当于是32bit分成8个8个，然后是灰度，红色，绿色，和蓝色，然后加在一起就是一个颜色
			//32个0和1，然后分成4部分，每部分8个部分，0-255 一共256个值
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
	//算出左边的y列和右边的y列，然后作为参数传过去，那计算的就是左到右每一行，每一行的颜色，然后再把这个作为红绿蓝，得到每一个点的color。
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
	//step size
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSpace = max(abs(xDiff),abs(yDiff));
	float xStepSize = xDiff/numberOfSpace;
	float yStepSize = yDiff/numberOfSpace;
	//color
	float red = colour.red;
	float green = colour.green;
	float blue = colour.blue;
	uint32_t color = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	//set color
	for (float i =0.0; i<numberOfSpace; i++ ){
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		window.setPixelColour(round(x), round(y), color);
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


//arrange the triangle:
CanvasTriangle arrangeTriangle(CanvasTriangle triangle){
	if(triangle.v0().y > triangle.v2().y) swap(triangle.v0(), triangle.v2());
	if(triangle.v0().y>triangle.v1().y) swap(triangle.v0(), triangle.v1());
	if(triangle.v1().y > triangle.v2().y) swap(triangle.v1(), triangle.v2());
	return triangle;
}


unsigned int textureMappercolor(DrawingWindow &window, CanvasTriangle triangle, CanvasPoint point, TextureMap textureMap) {
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

glm::mat3x3 affine(CanvasTriangle triangle){
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
	return affineT;

}

unsigned int affineTransformation(DrawingWindow &window, CanvasTriangle triangle, CanvasPoint point, TextureMap textureMap){
	//inverse of triangle
	TexturePoint p;
	glm::mat3x3 canvansP = {{point.x,0,0},{point.y,0,0},{1,0,0}};
	glm::mat3x3 affineT = affine(triangle);
	p.x = round((affineT*canvansP)[0][0]);
	p.y = round((affineT*canvansP)[1][0]);
	//std::cout<<glm::to_string(affineT)<<std::endl;
	uint32_t c = textureMap.pixels[((p.y)-1)*(int(textureMap.width))+p.x];
	return c;
}

void textureMapper(DrawingWindow &window, CanvasTriangle triangle, TextureMap textureMap){
	CanvasTriangle arranged = arrangeTriangle(triangle);
	CanvasPoint top = arranged.v0();
	CanvasPoint bottom = arranged.v2();
	CanvasPoint mid = arranged.v1();
	CanvasPoint ne; 
	ne.x= top.x+((mid.y-top.y)/(bottom.y-top.y))*(bottom.x-top.x);
	ne.y = mid.y;
	vector<glm::vec2> v_left = interpolateTwoElementValues({ne.x,ne.y},{top.x,top.y},(mid.y-top.y));
	vector<glm::vec2> v_right = interpolateTwoElementValues({mid.x,mid.y},{top.x,top.y},(mid.y-top.y));
	// for (int i=0;i<(abs(mid.x-ne.x));i++){
	// 	vector<glm::vec2> v = interpolateTwoElementValues({v_left[i].x,v_left[i].y},{v_right[i].x,v_right[i].y},abs(v_right[i].x-v_left[i].x));
	// 	for (int j=0;j<(abs(v_right[i].x)-v_left[i].x);i++){
	// 		uint32_t c = affineTransformation(window,triangle,{v[j].x,v[j].y},textureMap);
	// 		window.setPixelColour(v[j].x,v[j].y,c);
	// 	}
	// }
	for (int y = top.y; y<mid.y;y++){
		float x_left = top.x-((top.x-bottom.x)*(y-top.y)/(bottom.y-top.y));
		float x_right = top.x+((y-top.y)*(mid.x-top.x)/(mid.y-top.y));
		if (x_left > x_right) std::swap(x_left,x_right);
		for(int x=x_left; x<x_right; x++) {
			uint32_t c = textureMappercolor(window,arranged,CanvasPoint(x,y),textureMap);
			//uint32_t c = affineTransformation(window,arranged,CanvasPoint(x,y),textureMap);
			window.setPixelColour(x,y,c);
		}
	}
		for (int y=mid.y; y<bottom.y;y++){
		float x_left = top.x-((top.x-bottom.x)*(y-top.y)/(bottom.y-top.y));
		float x_right = mid.x-((y-mid.y)*(mid.x-bottom.x)/(bottom.y-mid.y));
		if (x_left > x_right) std::swap(x_left,x_right);
		for(int x=x_left; x<x_right; x++){
			uint32_t c = textureMappercolor(window,arranged,CanvasPoint(x,y),textureMap);
			//uint32_t c = affineTransformation(window,arranged,CanvasPoint(x,y),textureMap);

			window.setPixelColour(x,y,c);
		}
	}
}



void filledtriangle(DrawingWindow &window,CanvasTriangle triangle, Colour colour){
	CanvasTriangle arranged = arrangeTriangle(triangle);
	CanvasPoint top = arranged.v0();
	CanvasPoint bottom = arranged.v2();
	CanvasPoint mid = arranged.v1();
	//set color
	float red = colour.red;
	float green = colour.green;
	float blue = colour.blue;
	uint32_t color = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	//points fromup to down
	for (int y = top.y; y<mid.y;y++){
		float x_left = top.x-((top.x-bottom.x)*(y-top.y)/(bottom.y-top.y));
		float x_right = top.x+((y-top.y)*(mid.x-top.x)/(mid.y-top.y));
		if (x_left > x_right) std::swap(x_left,x_right);
		for(int x=x_left; x<x_right; x++) window.setPixelColour(x,y,color);
	}
	for (int y=mid.y; y<bottom.y;y++){
		float x_left = top.x-((top.x-bottom.x)*(y-top.y)/(bottom.y-top.y));
		float x_right = mid.x-((y-mid.y)*(mid.x-bottom.x)/(bottom.y-mid.y));
		if (x_left > x_right) std::swap(x_left,x_right);
		for(int x=x_left; x<x_right; x++) window.setPixelColour(x,y,color);
	}
}


void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) {std::cout << "DOWN" << std::endl;}
		else if (event.key.keysym.sym == SDLK_u) {
			CanvasTriangle triangle = strokedTriangles();
			Colour colour = {rand()%256,rand()%256,rand()%256};
		    lineDraw(window,triangle.v0(),triangle.v1(),colour);
			lineDraw(window,triangle.v1(),triangle.v2(),colour);
			lineDraw(window,triangle.v2(),triangle.v0(),colour);
		}
		else if (event.key.keysym.sym == SDLK_f) {
			CanvasTriangle triangle = strokedTriangles();
			Colour colour = {rand()%256,rand()%256,rand()%256};
			filledtriangle(window,triangle,colour);
			// white edge
			lineDraw(window,triangle.v0(),triangle.v1(),{255,255,255});
			lineDraw(window,triangle.v1(),triangle.v2(),{255,255,255});
			lineDraw(window,triangle.v2(),triangle.v0(),{255,255,255});
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	//week02 task 2
	// std::vector<float> result1;
	// result1 = interpolateSingleFloats(2.2, 8.5, 7);
	// for(size_t i=0; i<result1.size(); i++) std::cout << result1[i] << " ";
	// std::cout << std::endl;
	//week02 task 4
	// std::vector<glm::vec3> result2;
	// glm::vec3 from(1.0, 4.0, 9.2);
	// glm::vec3 to(4.0, 1.0, 9.8);
    //result2 = interpolateThreeElementValues(from,to, 4);
	//for (int i = 0; i < (int)result2.size(); i++) std::cout<<glm::to_string(result2[i])<<std::endl;
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	float height = window.height;
	float width = window.width;
	SDL_Event event;
	char c;
	CanvasTriangle canvas = {{160,10},{300,230},{10,150}};
	canvas.v0().texturePoint = {195,5};
	canvas.v1().texturePoint = {395,380};
	canvas.v2().texturePoint = {65,330};

	TextureMap file = TextureMap("texture.ppm");
		textureMapper(window,canvas,file);
		lineDraw(window,canvas.v0(),canvas.v1(),{255,255,255});
		lineDraw(window,canvas.v1(),canvas.v2(),{255,255,255});
 	    lineDraw(window,canvas.v2(),canvas.v0(),{255,255,255});

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) {
			handleEvent(event, window);
		}
		
		//std::cout<<file<<std::endl;
		//draw(window);
		//week02 task 3
		//greyScaleInterpolation(window);
		//dimensionalColourInterpolation(window);
		/* week03 task 2
		lineDraw(window,{0,0},{width/2,height/2},{0, 255, 0});
		lineDraw(window,{width/2,height/2},{width/2,height},{255,0,0});
		lineDraw(window,{width/3,height/2},{2*width/3,height/2},{0, 255, 0});
		lineDraw(window,{2*width/3,height/2},{width/3,height/2},{0, 255, 0});
		lineDraw(window,{width/2,height/3},{width/2,2*height/3},{0, 0, 255});*/
		// Need to render thse frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
