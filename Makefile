default: graphics

%: %.cpp
	g++ -O3 -I. $< -o $@ shader.cpp -lGLEW  -lGL -lglfw 
clean:
	rm a.out *.o *~ graphics
