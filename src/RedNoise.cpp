#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#define WIDTH 320
#define HEIGHT 240
using namespace std;

/*void draw(DrawingWindow &window) {
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
}*/

std::vector<float> interpolateSingleFloats(float f, float to, int N){
	 vector<float> v;
	 float next = f;
	 float a = (to - f) / (N - 1);
	 for (int i = 0; i < N; i++) {
		v.push_back(next);
		next = next + a;
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

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from,glm::vec3 to, int N){
	vector<glm::vec3> v;
	vector<glm::vec3> next = from;
	 vector<glm::vec3> a = (from - to) / (N - 1);
	 for (int i = 0; i < N; i++) {
		v.push_back(next);
		next = next + a;
	 }
	 return v;

}


void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	std::vector<glm::vec3> result;
    result = interpolateSingleFloats((1.0, 4.0, 9.2), (4.0, 1.0, 9.8), 4);
    for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
	std::cout << std::endl;
	//DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		greyScaleInterpolation(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
	
}
