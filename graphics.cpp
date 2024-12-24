#define pb push_back

// Include standard headers
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <GLFW/glfw3.h>
GLFWwindow* window;

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
using namespace std;

#include <shader.hpp>

typedef long long ll;

const ll inf = 1000000000000000009;
const float earth_radius = 6371.0;
const int KEY_N = 512;
bool key_toggle[KEY_N];
int window_width = 768;
int window_height = 768;
float window_ratio = 1.0*window_width/window_height;
int lod_mode = -1;
float zoom, zoom3d;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if(key >= KEY_N) return;
	if(action == GLFW_PRESS) {
		key_toggle[key] = !key_toggle[key];
		
		if(key == GLFW_KEY_0)
			lod_mode = -1;
		else if(key == GLFW_KEY_1)
			lod_mode = 0;
		else if(key == GLFW_KEY_2)
			lod_mode = 1;
		else if(key == GLFW_KEY_3)
			lod_mode = 2;
		else if(key == GLFW_KEY_4)
			lod_mode = 3;
		else if(key == GLFW_KEY_5)
			lod_mode = 4;
	}
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
	zoom *= pow(1.125, yoffset);
	zoom3d *= pow(1.0+1.0/16.0, yoffset);
	if(zoom3d < 0.5) zoom3d = 0.5;
	if(zoom3d > 16.0) zoom3d = 16.0;
}

void window_size_callback(GLFWwindow* window, int width, int height) {
	window_width = width;
	window_height = height;
	window_ratio = 1.0*window_width/window_height;
	glViewport(0, 0, width, height);
}

GLuint index_buffer[5];
GLuint model_buffer;
int index_buffer_size[5];

const int SRTM_SIZE = 1201;
GLuint shader2D, shader3D;

ll triangle_count = 0;

class Tile {
public:
	
	bool loaded;
	int lon, lat;
	GLuint vertex_buffer, normal_buffer;
	
	Tile() {
		
		glGenBuffers(1, &vertex_buffer);
		glGenBuffers(1, &normal_buffer);
		loaded = false;
	}
	
	void load_from_file(string filename) {
	
		loaded = false;
		ifstream file(filename, ios::in|ios::binary);
		if(!file) return;
		unsigned char buffer[2];
		GLfloat* height = new GLfloat[SRTM_SIZE*SRTM_SIZE];
		for (int i = 0; i < SRTM_SIZE; ++i)
		{
		    for (int j = 0; j < SRTM_SIZE; ++j) 
		    {
		        if(!file.read( reinterpret_cast<char*>(buffer), sizeof(buffer) )) return;
		        height[i*SRTM_SIZE+j] = (buffer[0] << 8) | buffer[1];
		    }
		}
		GLfloat* normal_data = new GLfloat[3*SRTM_SIZE*SRTM_SIZE];
		for(int i = 1; i < SRTM_SIZE-1; ++i) {
			for(int j = 1; j < SRTM_SIZE-1; ++j) {
			
				float lon_this = lon+1.0*j/1200;
				float lat_this = lat+1.0*(1200-i)/1200;
				glm::vec3 pos_this = glm::vec3(sin(glm::radians(lon_this))*cos(glm::radians(lat_this)), sin(glm::radians(lat_this)), cos(glm::radians(lon_this))*cos(glm::radians(lat_this)));
				
				glm::vec3 east = normalize(cross(glm::vec3(0, 1, 0), pos_this));
				glm::vec3 north = cross(pos_this, east);
				
				pos_this *= earth_radius+height[i*SRTM_SIZE+j]/1000;
				
				// North
				float lon_north = lon+1.0*j/1200;
				float lat_north = lat+1.0*(1200-(i-1))/1200;
				glm::vec3 pos_north = (earth_radius+height[(i-1)*SRTM_SIZE+j]/1000)*glm::vec3(sin(glm::radians(lon_north))*cos(glm::radians(lat_north)), sin(glm::radians(lat_north)), cos(glm::radians(lon_north))*cos(glm::radians(lat_north)));
				
				glm::vec3 norm_north = cross(east, normalize(pos_north-pos_this));
				
				// South
				float lon_south = lon+1.0*j/1200;
				float lat_south = lat+1.0*(1200-(i+1))/1200;
				glm::vec3 pos_south = (earth_radius+height[(i+1)*SRTM_SIZE+j]/1000)*glm::vec3(sin(glm::radians(lon_south))*cos(glm::radians(lat_south)), sin(glm::radians(lat_south)), cos(glm::radians(lon_south))*cos(glm::radians(lat_south)));
				
				glm::vec3 norm_south = cross(-east, normalize(pos_south-pos_this));
				
				// East
				float lon_east = lon+1.0*(j+1)/1200;
				float lat_east = lat+1.0*(1200-i)/1200;
				glm::vec3 pos_east = (earth_radius+height[i*SRTM_SIZE+j+1]/1000)*glm::vec3(sin(glm::radians(lon_east))*cos(glm::radians(lat_east)), sin(glm::radians(lat_east)), cos(glm::radians(lon_east))*cos(glm::radians(lat_east)));
				
				glm::vec3 norm_east = cross(-north, normalize(pos_east-pos_this));
				
				// West
				float lon_west = lon+1.0*(j-1)/1200;
				float lat_west = lat+1.0*(1200-i)/1200;
				glm::vec3 pos_west = (earth_radius+height[i*SRTM_SIZE+j-1]/1000)*glm::vec3(sin(glm::radians(lon_west))*cos(glm::radians(lat_west)), sin(glm::radians(lat_west)), cos(glm::radians(lon_west))*cos(glm::radians(lat_west)));
				
				glm::vec3 norm_west = cross(north, normalize(pos_west-pos_this));
				
				
				// Normal vector
				glm::vec3 normal = normalize(norm_north+norm_south+norm_east+norm_west);
				//normal = east;
				normal_data[3*(i*SRTM_SIZE+j)+0] = normal.x;
				normal_data[3*(i*SRTM_SIZE+j)+1] = normal.y;
				normal_data[3*(i*SRTM_SIZE+j)+2] = normal.z;
			}
		}
		glBindBuffer(GL_ARRAY_BUFFER, normal_buffer);
		glBufferData(GL_ARRAY_BUFFER, 3*SRTM_SIZE*SRTM_SIZE*sizeof(GLfloat), normal_data, GL_STATIC_DRAW);
		delete[] normal_data;
		
		glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
		glBufferData(GL_ARRAY_BUFFER, SRTM_SIZE*SRTM_SIZE*sizeof(GLfloat), height, GL_STATIC_DRAW);
		delete[] height;
		loaded = true;
	}
	
	void draw2D(glm::vec3 pos, int lod) {
	
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
		glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 0, (void*)0);
		
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, model_buffer);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
		
		glUniform3f(glGetUniformLocation(shader2D, "pos"), pos.x, pos.y, pos.z);
		glUniform2f(glGetUniformLocation(shader2D, "tilePos"), lon, lat);
		glUniform1f(glGetUniformLocation(shader2D, "windowRatio"), window_ratio);
		
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer[lod]);
		glDrawElements(GL_TRIANGLES, index_buffer_size[lod], GL_UNSIGNED_INT, (void*)0);
		triangle_count += index_buffer_size[lod]/3;
		
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
	}
	
	void draw3D(glm::mat4 VP, glm::vec3 light, int lod) {
	
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
		glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 0, (void*)0);
		
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, model_buffer);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
		
		glEnableVertexAttribArray(2);
		glBindBuffer(GL_ARRAY_BUFFER, normal_buffer);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		
		glUniformMatrix4fv(glGetUniformLocation(shader3D, "VP"), 1, GL_FALSE, &VP[0][0]);
		glUniform2f(glGetUniformLocation(shader3D, "tilePos"), lon, lat);
		glUniform3f(glGetUniformLocation(shader3D, "light"), light.x, light.y, light.z);
		
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer[lod]);
		glDrawElements(GL_TRIANGLES, index_buffer_size[lod], GL_UNSIGNED_INT, (void*)0);
		triangle_count += index_buffer_size[lod]/3;
		
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(2);
	}
};

int main(int argc, const char *argv[])
{
	int min_longitude = -180;
	int max_longitude = 180;
	int min_latitude = -90;
	int max_latitude = 90;

	int argptr = 1;
	if(argptr >= argc) {
		cout << "Usage: zadanie6 [data catalog name]\n";
		return 0;
	}
	argptr++;
	float init_height, init_lon, init_lat;
	bool start_set = false;
	while(argptr < argc) {
		string s = argv[argptr];
		if(s == "-lon") {
			if(argptr+2 < argc) {
				min_longitude = atoi(argv[argptr+1]);
				max_longitude = atoi(argv[argptr+2]);
				argptr += 3;
			}
			else {
				cout << "Usage: -lon n m, where n, m range from -180 to 180\n";
				return 0;
			}
		}
		else if(s == "-lat") {
			if(argptr+2 < argc) {
				min_latitude = atoi(argv[argptr+1]);
				max_latitude = atoi(argv[argptr+2]);
				argptr += 3;
			}
			else {
				cout << "Usage: -lat n m, where n, m range from -90 to 90\n";
				return 0;
			}
		}
		else if(s == "-start") {
			if(argptr+3 < argc) {
				init_lon = atof(argv[argptr+1]);
				init_lat = atof(argv[argptr+2]);
				init_height = atof(argv[argptr+3])/1000.0;
				argptr += 4;
				start_set = true;
			}
			else {
				cout << "Usage: -start lon lat h\n";
				return 0;
			}
		}
		else {
			cout << "Error\n";
			return 0;
		}
	}
	cout << "longitude: [" << min_longitude << ", " << max_longitude << "]\n";
	cout << "latitude: [" << min_latitude << ", " << max_latitude << "]\n";
	glm::vec3 position2D = glm::vec3((min_longitude+max_longitude)/2.0, (min_latitude+max_latitude)/2.0, 1);
	float height = 1;
	if(start_set) {
		height = init_height;
		position2D = glm::vec3(init_lon, init_lat, 1);
	}

	// Initialise GLFW
	if( !glfwInit() )
	{
		fprintf( stderr, "Failed to initialize GLFW\n" );
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(window_width, window_height, "Zadanie 6", NULL, NULL);
	if( window == NULL ){
		fprintf( stderr, "Failed to open GLFW window\n" );
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}

	// Ensure we can capture the escape key being pressed below
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    glfwSetKeyCallback(window, key_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetWindowSizeCallback(window, window_size_callback);
    // Hide the mouse and enable unlimited mouvement
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    
    // Set the mouse at the center of the screen
    glfwPollEvents();
    glfwSetCursorPos(window, window_width/2, window_height/2);
    glfwSwapInterval(0);

	glClearColor(0.529f, 0.808f, 0.922f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);
	
	GLfloat* model_data = new GLfloat[2*SRTM_SIZE*SRTM_SIZE];
	for (int i = 0; i < SRTM_SIZE; ++i)
	{
		for (int j = 0; j < SRTM_SIZE; ++j) 
		{
			model_data[2*(i*SRTM_SIZE+j)+0] = j/1200.0;
			model_data[2*(i*SRTM_SIZE+j)+1] = (1200-i)/1200.0;
		}
	}
	glGenBuffers(1, &model_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, model_buffer);
	glBufferData(GL_ARRAY_BUFFER, 2*SRTM_SIZE*SRTM_SIZE*sizeof(GLfloat), &model_data[0], GL_STATIC_DRAW);
	delete[] model_data;
	
	// LOD = 0
	int* index_data = new int[6*1200*1200]; 
	int k = 0;
	for(int i = 0; i < 1200; ++i) {
		for(int j = 0; j < 1200; ++j) {
		
			index_data[k++] = SRTM_SIZE*i+j;
			index_data[k++] = SRTM_SIZE*(i+1)+j;
			index_data[k++] = SRTM_SIZE*i+j+1;
			
			index_data[k++] = SRTM_SIZE*(i+1)+j;
			index_data[k++] = SRTM_SIZE*i+j+1;
			index_data[k++] = SRTM_SIZE*(i+1)+j+1;
		}
	}
	glGenBuffers(1, &index_buffer[0]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer[0]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, k*sizeof(int), &index_data[0], GL_STATIC_DRAW);
	index_buffer_size[0] = k;
	
	// LOD = 1
	k = 0;
	for(int i = 0; i < 1200; i += 2) {
		for(int j = 0; j < 1200; j += 2) {
		
			index_data[k++] = SRTM_SIZE*i+j;
			index_data[k++] = SRTM_SIZE*(i+2)+j;
			index_data[k++] = SRTM_SIZE*i+j+2;
			
			index_data[k++] = SRTM_SIZE*(i+2)+j;
			index_data[k++] = SRTM_SIZE*i+j+2;
			index_data[k++] = SRTM_SIZE*(i+2)+j+2;
		}
	}
	glGenBuffers(1, &index_buffer[1]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer[1]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, k*sizeof(int), &index_data[0], GL_STATIC_DRAW);
	index_buffer_size[1] = k;
	
	// LOD = 2
	k = 0;
	for(int i = 0; i < 1200; i += 4) {
		for(int j = 0; j < 1200; j += 4) {
		
			index_data[k++] = SRTM_SIZE*i+j;
			index_data[k++] = SRTM_SIZE*(i+4)+j;
			index_data[k++] = SRTM_SIZE*i+j+4;
			
			index_data[k++] = SRTM_SIZE*(i+4)+j;
			index_data[k++] = SRTM_SIZE*i+j+4;
			index_data[k++] = SRTM_SIZE*(i+4)+j+4;
		}
	}
	glGenBuffers(1, &index_buffer[2]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer[2]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, k*sizeof(int), &index_data[0], GL_STATIC_DRAW);
	index_buffer_size[2] = k;
	
	// LOD = 3
	k = 0;
	for(int i = 0; i < 1200; i += 8) {
		for(int j = 0; j < 1200; j += 8) {
		
			index_data[k++] = SRTM_SIZE*i+j;
			index_data[k++] = SRTM_SIZE*(i+8)+j;
			index_data[k++] = SRTM_SIZE*i+j+8;
			
			index_data[k++] = SRTM_SIZE*(i+8)+j;
			index_data[k++] = SRTM_SIZE*i+j+8;
			index_data[k++] = SRTM_SIZE*(i+8)+j+8;
		}
	}
	glGenBuffers(1, &index_buffer[3]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer[3]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, k*sizeof(int), &index_data[0], GL_STATIC_DRAW);
	index_buffer_size[3] = k;
	
	// LOD = 4
	k = 0;
	for(int i = 0; i < 1200; i += 16) {
		for(int j = 0; j < 1200; j += 16) {
		
			index_data[k++] = SRTM_SIZE*i+j;
			index_data[k++] = SRTM_SIZE*(i+16)+j;
			index_data[k++] = SRTM_SIZE*i+j+16;
			
			index_data[k++] = SRTM_SIZE*(i+16)+j;
			index_data[k++] = SRTM_SIZE*i+j+16;
			index_data[k++] = SRTM_SIZE*(i+16)+j+16;
		}
	}
	glGenBuffers(1, &index_buffer[4]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer[4]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, k*sizeof(int), &index_data[0], GL_STATIC_DRAW);
	index_buffer_size[4] = k;
	
	delete[] index_data;
	
	// Tiles
	shader2D = LoadShaders("shaders/2d.vertexshader", "shaders/2d.fragmentshader");
	shader3D = LoadShaders("shaders/3d.vertexshader", "shaders/3d.fragmentshader");
	string catalog = argv[1];
	vector<Tile> tiles(1);
	int tile_count = 0;
	
	for(int i = min_longitude; i < max_longitude; ++i) {
		for(int j = min_latitude; j < max_latitude; ++j) {
			string filename = catalog;
			if(j < 0) {
				filename += "S";
				if(-j < 10) filename += "0";
				filename += to_string(-j);
			}
			else {
				filename += "N";
				if(j < 10) filename += "0";
				filename += to_string(j);
			}
			if(i < 0) {
				filename += "W";
				if(-i < 10) filename += "0";
				if(-i < 100) filename += "0";
				filename += to_string(-i);
			}
			else {
				filename += "E";
				if(i < 10) filename += "0";
				if(i < 100) filename += "0";
				filename += to_string(i);
			}
			filename += ".hgt";
			tiles[tile_count].lon = i;
			tiles[tile_count].lat = j;
			tiles[tile_count].load_from_file(filename);
			if(tiles[tile_count].loaded) {
				tiles.pb(Tile());
				tile_count++;
			}
		}
	}
	cout << "Loaded " << tile_count << " tiles\n";
	
	// Sphere
	GLuint sphere_shader = LoadShaders("shaders/sphere.vertexshader", "shaders/sphere.fragmentshader");
	GLfloat* sphere_data = new GLfloat[3*(361*180+179*361)];
	k = 0;
	for(int i = 0; i < 180; ++i) {
		for(int j = 0; j <= 360; ++j) {
			float sin_lat = sin(glm::radians(1.0*j));
			float cos_lat = cos(glm::radians(1.0*j));
			float sin_lon = sin(glm::radians(1.0*i));
			float cos_lon = cos(glm::radians(1.0*i));
			sphere_data[k++] = sin_lon*cos_lat;
			sphere_data[k++] = sin_lat;
			sphere_data[k++] = cos_lon*cos_lat;
		}
	}
	for(int i = -89; i <= 89; ++i) {
		for(int j = 0; j <= 360; ++j) {
			float sin_lat = sin(glm::radians(1.0*i));
			float cos_lat = cos(glm::radians(1.0*i));
			float sin_lon = sin(glm::radians(1.0*j));
			float cos_lon = cos(glm::radians(1.0*j));
			sphere_data[k++] = sin_lon*cos_lat;
			sphere_data[k++] = sin_lat;
			sphere_data[k++] = cos_lon*cos_lat;
		}
	}
	GLuint sphere_buffer;
	glGenBuffers(1, &sphere_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, sphere_buffer);
	glBufferData(GL_ARRAY_BUFFER, k*sizeof(GLfloat), sphere_data, GL_STATIC_DRAW);
	delete[] sphere_data;
	
	// Camera
	zoom = min(2.0/(max_latitude-min_latitude),
			   window_ratio/cos(glm::radians(position2D.y))/((max_longitude-min_longitude)/2.0));
	position2D.z = zoom;
	
	const float init_fov = 45.0;
	float mouse_speed = 0.001;
	float hor_angle = 0.0;
	float ver_angle = 0.0;
	float near_clipping = 0.064;
	float far_clipping = 512.0;
	float sin_lat = sin(glm::radians(position2D.y));
	float cos_lat = cos(glm::radians(position2D.y));
	float sin_lon = sin(glm::radians(position2D.x));
	float cos_lon = cos(glm::radians(position2D.x));
	glm::vec3 position3D = glm::vec3(sin_lon*cos_lat, sin_lat, cos_lon*cos_lat);
	zoom3d = 1;

	// For speed computation
	double lastTime = glfwGetTime();
	double lastTimeFrames = lastTime;
	int nbFrames = 0;
	bool running = true;
	int perf_lod = 0;
	//ll perf_triangles = inf;

	do{

		// Measure speed
		double currentTime = glfwGetTime();
		float deltaTime = float(currentTime - lastTime);
		lastTime = currentTime;
		nbFrames++;
		if ( currentTime - lastTimeFrames >= 1.0 ) {
			// printf and reset
			printf("FPS: %d\n", nbFrames);
			printf("avg triangles per frame: %lld\n", triangle_count/nbFrames);
			if(nbFrames <= 10) {
				//perf_triangles = triangle_count;
				if(perf_lod < 4) perf_lod++;
			}
			else if(perf_lod > 0) {
				//perf_triangles = inf;
				perf_lod--;
			}
			nbFrames = 0;
			triangle_count = 0;
			lastTimeFrames += 1.0;
		}

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		if(key_toggle[GLFW_KEY_D])
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		else
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			
		if(key_toggle[GLFW_KEY_TAB]) {
		
			if(glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
				ver_angle -= deltaTime/zoom3d;
			if(glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
				ver_angle += deltaTime/zoom3d;
			if(glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
				hor_angle += deltaTime/zoom3d;
			if(glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
				hor_angle -= deltaTime/zoom3d;
				
			if(ver_angle <= -M_PI/2) ver_angle = -M_PI/2+0.01;
			if(ver_angle >= M_PI/2) ver_angle = M_PI/2-0.01;
			
			if(glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
				height *= pow(1.125, deltaTime);
			if(glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
				height /= pow(1.125, deltaTime);
				
			if(height < 0.016) height = 0.016;
			
			// Camera
			glm::vec3 east = normalize(cross(glm::vec3(0, 1, 0), position3D));
			glm::vec3 north = cross(position3D, east);
			glm::vec3 right = cos(hor_angle)*east+sin(hor_angle)*north;
			glm::vec3 dir = -sin(hor_angle)*east+cos(hor_angle)*north;
			glm::vec3 look_dir = cos(ver_angle)*dir+sin(ver_angle)*position3D;
			
			if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
				position3D = normalize((earth_radius+height)*position3D+height*deltaTime*dir);
			if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
				position3D = normalize((earth_radius+height)*position3D-height*deltaTime*dir);
				
			float sin_lat = position3D.y;
			float cos_lat = sqrt(1-sin_lat*sin_lat);
			float sin_lon = position3D.x/cos_lat;
			float cos_lon = position3D.z/cos_lat;
			float lat = asin(sin_lat);
			float lon = asin(sin_lon);
			if(cos_lon < 0 && sin_lon > 0) lon = M_PI-lon;
			if(cos_lon < 0 && sin_lon < 0) lon = -M_PI-lon;
			position2D = glm::vec3(glm::degrees(lon), glm::degrees(lat), zoom);
			
			// Light
			glm::vec3 out = glm::cross(east, glm::vec3(0, 1, 0));
			glm::vec3 light = glm::normalize(out+east);
			
			// Projection matrix
			glm::mat4 ProjectionMatrix = glm::perspective(glm::radians(init_fov/zoom3d), window_ratio, near_clipping, far_clipping);
			glm::mat4 ViewMatrix = lookAt((earth_radius+height)*position3D, (earth_radius+height)*position3D+look_dir, position3D);
			glm::mat4 VP = ProjectionMatrix * ViewMatrix;
			
			int lod = perf_lod;
			if(0 <= lod_mode && lod_mode <= 4) {
				lod = lod_mode;
			}
			
			glUseProgram(sphere_shader);
			
			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, sphere_buffer);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
			
			glUniformMatrix4fv(glGetUniformLocation(sphere_shader, "VP"), 1, GL_FALSE, &VP[0][0]);
			glDrawArrays(GL_LINE_STRIP, 0, 180*361+179*361);
			
			glDisableVertexAttribArray(0);
			
			glUseProgram(shader3D);
			
			for(int i = 0; i < tile_count; ++i) {
				
				float dx = min(abs(tiles[i].lon-position2D.x), abs(tiles[i].lon+1-position2D.x));
				float dy = min(abs(tiles[i].lat-position2D.y), abs(tiles[i].lat+1-position2D.y));
				float d = max(dx, dy);
				
				if(d <= 2) tiles[i].draw3D(VP, light, lod);
			}
		}
		else {
			if(glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
				position2D.y -= deltaTime/zoom;
			if(glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
				position2D.y += deltaTime/zoom;
			if(glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
				position2D.x -= deltaTime/zoom/cos(glm::radians(position2D.y));
			if(glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
				position2D.x += deltaTime/zoom/cos(glm::radians(position2D.y));
				
			if(position2D.y > 90.0) position2D.y = 90.0;
			if(position2D.y < -90.0) position2D.y = -90.0;
			if(position2D.x > 180.0) position2D.y = 180.0;
			if(position2D.x < -180.0) position2D.y = -180.0;
			
			float sin_lat = sin(glm::radians(position2D.y));
			float cos_lat = cos(glm::radians(position2D.y));
			float sin_lon = sin(glm::radians(position2D.x));
			float cos_lon = cos(glm::radians(position2D.x));
			position3D = glm::vec3(sin_lon*cos_lat, sin_lat, cos_lon*cos_lat);
				
			position2D.z = zoom;
				
			glUseProgram(shader2D);
			
			int lod = perf_lod;
			if(0 <= lod_mode && lod_mode <= 4) {
				lod = lod_mode;
			}
			
			for(int i = 0; i < tile_count; ++i) {
				float min_x = position2D.z*cos(glm::radians(position2D.y))*(tiles[i].lon-position2D.x)/window_ratio;
				float min_y = position2D.z*(tiles[i].lat-position2D.y);
				
				float max_x = position2D.z*cos(glm::radians(position2D.y))*(1+tiles[i].lon-position2D.x)/window_ratio;
				float max_y = position2D.z*(1+tiles[i].lat-position2D.y);
			
				if(max_x > -1 && max_y > -1 && min_x < 1 && min_y < 1)
					tiles[i].draw2D(position2D, lod);
			}
		}

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
		   glfwWindowShouldClose(window) == 0 && running);

	// Cleanup VBO and shader
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}
