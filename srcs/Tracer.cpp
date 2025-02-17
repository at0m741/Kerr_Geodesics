#include <Geodesics.h>
#include <GL/glu.h>
#include <GLFW/glfw3.h>
#define WIDTH 1920
#define HEIGHT 1080
#define R 200.0  

extern double (*geodesic_points)[5];
extern int num_points;
extern double a;
/* float rotation_angle = 0.0f; */
int displayed_points = 0;



void draw_blackhole() {
    glPointSize(1.0f);
    glBegin(GL_POINTS);

    for (int i = 0; i < displayed_points; i++) {
        double x = geodesic_points[i][0];
        double y = geodesic_points[i][1];
        double z = geodesic_points[i][2];

        double r = sqrt(x * x + y * y);
        double r_horizon = 2.0;

        float t = (float)i / num_points;  
        glColor3f(1.0f - t, t, 0.5f);  

        glVertex3d(x, y, z);
    }

    glEnd();
}



void process_input(GLFWwindow *window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, 1);
    }
}

void generate_blackhole_image() {
    if (!glfwInit()) {
        printf("Erreur: Impossible d'initialiser GLFW\n");
        return;
    }

    GLFWwindow *window = glfwCreateWindow(WIDTH, HEIGHT, "Black Hole Visualization", NULL, NULL);
    if (!window) {
        printf("Erreur: Impossible de créer la fenêtre OpenGL\n");
        glfwTerminate();
        return;
    }

    glfwMakeContextCurrent(window);
    glEnable(GL_DEPTH_TEST);

    while (!glfwWindowShouldClose(window)) {
        process_input(window);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
		glFrustum(-0.5, 0.5, -0.5, 0.5, 0.1, 2000.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        gluLookAt(80, 80, 80, 
                  0, 0, 0,      
                  0, 1, 0);    

        draw_blackhole();

        if (displayed_points < num_points) {
            displayed_points += 65000;
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
}
