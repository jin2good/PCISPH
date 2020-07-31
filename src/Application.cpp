#include "Global.h"

#include "renderer.h"
#include "shader.h"
#include "VertexBuffer.h"
#include "VertexBufferLayout.h"
#include "IndexBuffer.h"
#include "VertexArray.h"
#include "Camera.h"
#include "ParticleSystem.h"
#include "Parameter.h"


void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);

/* Settings */
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

/* Camera */
Camera camera(glm::vec3(0.0f, 0.0f, 5.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;
bool pause = true;
int leftmouse = GLFW_RELEASE;

/* timing */
double deltaTime = 0.0f;	// time between current frame and last frame
double lastFrame = 0.0f;

int main(void)
{
    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */

    window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "SPH Simulation", NULL, NULL);
    if (!window)    
    {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    leftmouse = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    //glfwSwapInterval(1);


    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    if (glewInit() != GLEW_OK) {
        std::cout << "ERROR: failed to initialize GLEW" << std::endl;
    }

    {
        // Global Opengl Configuration
        glEnable(GL_DEPTH_TEST);

        /* Box */
        float boxpos[24] = {
            0.0f,   0.0f,  0.0f,
            PARTICLE_INITIAL_BOUNDARY_X, 0.0f,                        0.0f,
            PARTICLE_INITIAL_BOUNDARY_X, 0.0f,                        PARTICLE_INITIAL_BOUNDARY_Z,
            0.0f,                        0.0f,                        PARTICLE_INITIAL_BOUNDARY_Z,
            0.0f,                        PARTICLE_INITIAL_BOUNDARY_Y, 0.0f,
            PARTICLE_INITIAL_BOUNDARY_X, PARTICLE_INITIAL_BOUNDARY_Y, 0.0f,
            PARTICLE_INITIAL_BOUNDARY_X, PARTICLE_INITIAL_BOUNDARY_Y, PARTICLE_INITIAL_BOUNDARY_Z,
            0.0f,                        PARTICLE_INITIAL_BOUNDARY_Y, PARTICLE_INITIAL_BOUNDARY_Z,
        };
        unsigned int boxind[30] = {
            0, 2, 1,
            0, 3, 2,
            0, 1, 5,
            0, 5, 4,
            1, 2, 6,
            1, 6, 5,
            3, 6, 2,
            3, 7, 6,
            0, 4, 7,
            0, 7, 3,
        };
        VertexArray box;
        VertexBuffer box_vb(boxpos, 24 * sizeof(float));
        VertexBufferLayout box_layout;
        box_layout.Push<float>(3);
        box.AddBuffer(box_vb, box_layout);

        IndexBuffer box_ind(boxind, 30);

        /* particle system */
        ParticleSystem* ps = new ParticleSystem(GRANULAR);
        float* particle_position = ps->GetParticlePositionArray();


        VertexArray p_va;
        

        VertexBuffer p_vb(particle_position, (SHOW_PARTICLE_BOX ? ps->particle_count : PARTICLE_COUNT) * 3 * sizeof(float));
        VertexBufferLayout p_layout;
        p_layout.Push<float>(3);
        p_va.AddBuffer(p_vb, p_layout);

        /* Bound Shader*/
        Shader shader("shaders//vertex.shader", "shaders//fragment.shader");

        Renderer renderer;
        
            while (!glfwWindowShouldClose(window))
            {
                renderer.Clear();
                // per-frame time logic
                double currentFrame = glfwGetTime();
                deltaTime = currentFrame - lastFrame;
                lastFrame = currentFrame;
             
                // input
                if(ENABLE_KEY_INPUT)
                    processInput(window);
                
                /* Activate Shader*/
                shader.Bind();

                // pass projection matrix to shader (note that in this case it could change every frame)
                glm::mat4 projectionMatrix = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
                shader.SetUniformMat4f("projMat", projectionMatrix);

                // camera/view transformation
                glm::mat4 viewMatrix = camera.GetViewMatrix();
                shader.SetUniformMat4f("viewMat", viewMatrix);

                // for each object
                glm::mat4 worldMatrix = glm::mat4(1.0f); // for now
                shader.SetUniformMat4f("worldMat", worldMatrix);

                /* Update */
                if (!pause) {
                    double updatetime = dt;
                    if(ENABLE_DEBUG_MODE)
                        std::cout << "time taken: " << deltaTime << std::endl;
                    ps->Update(updatetime);
                }
                    

                float *ptr = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);
                if (ptr) {
                    float* pos = ps->GetParticlePositionArray();
                    
                    for (int i = 0; i < (SHOW_PARTICLE_BOX ? ps->particle_count : PARTICLE_COUNT) * 3; ++i) {
                        *ptr = *pos; ++pos; ++ptr;
                    }
                    glUnmapBuffer(GL_ARRAY_BUFFER);
                }                
                
                /* Draw Box */
                if (SHOW_BOX) {
                    glLineWidth(2.0f);
                    glColor3f(1.0f, 0.0f, 0.0f);
                    shader.SetUniform4f("u_Color", 1.0f, 1.0f, 1.0f, 1.0f);
                    glBegin(GL_LINES);
                        

                    glVertex3f(0.0f, 0.0f, 0.0f);
                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, 0.0f, 0.0f);

                    glVertex3f(0.0f, 0.0f, 0.0f);
                    glVertex3f(0.0f, PARTICLE_INITIAL_BOUNDARY_Y, 0.0f);

                    glVertex3f(0.0f, 0.0f, 0.0f);
                    glVertex3f(0.0f, 0.0f, PARTICLE_INITIAL_BOUNDARY_Z);

                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, 0.0f, 0.0f);
                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, PARTICLE_INITIAL_BOUNDARY_Y, 0.0f);

                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, PARTICLE_INITIAL_BOUNDARY_Y, 0.0f);
                    glVertex3f(0.0f, PARTICLE_INITIAL_BOUNDARY_Y, 0.0f);

                    glVertex3f(0.0f, PARTICLE_INITIAL_BOUNDARY_Y, PARTICLE_INITIAL_BOUNDARY_Z);
                    glVertex3f(0.0f, 0.0f, PARTICLE_INITIAL_BOUNDARY_Z);

                    glVertex3f(0.0f, PARTICLE_INITIAL_BOUNDARY_Y, PARTICLE_INITIAL_BOUNDARY_Z);
                    glVertex3f(0.0f, PARTICLE_INITIAL_BOUNDARY_Y, 0.0f);

                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, 0.0f, 0.0f);
                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, 0.0f, PARTICLE_INITIAL_BOUNDARY_Z);

                    glVertex3f(0.0f, 0.0f, PARTICLE_INITIAL_BOUNDARY_Z);
                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, 0.0f, PARTICLE_INITIAL_BOUNDARY_Z);

                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, PARTICLE_INITIAL_BOUNDARY_Y, 0.0f);
                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, PARTICLE_INITIAL_BOUNDARY_Y, PARTICLE_INITIAL_BOUNDARY_Z);

                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, PARTICLE_INITIAL_BOUNDARY_Y, PARTICLE_INITIAL_BOUNDARY_Z);
                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, 0.0f, PARTICLE_INITIAL_BOUNDARY_Z);

                    glVertex3f(0.0f, PARTICLE_INITIAL_BOUNDARY_Y, PARTICLE_INITIAL_BOUNDARY_Z);
                    glVertex3f(PARTICLE_INITIAL_BOUNDARY_X, PARTICLE_INITIAL_BOUNDARY_Y, PARTICLE_INITIAL_BOUNDARY_Z);

                    glEnd();
                }
                
                /* Draw Points */
                shader.SetUniform4f("u_Color", 0.6f, 0.6f, 0.2f, 1.0f);
                renderer.DrawPoints(p_va, shader, (SHOW_PARTICLE_BOX ? ps->particle_count : PARTICLE_COUNT));
                
                ptr = ps->GetParticlePositionArray();

                /* Swap front and back buffers */
                glfwSwapBuffers(window);

                /* Poll for and process events */
                if (ENABLE_KEY_INPUT)
                    glfwPollEvents();
            }
    }
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
        pause = !pause;
}


// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}


// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = (float)xpos;
        lastY = (float)ypos;
        firstMouse = false;
    }
    
    float xoffset = (float)xpos - lastX;
    float yoffset = lastY - (float)ypos; // reversed since y-coordinates go from bottom to top

    lastX = (float)xpos;
    lastY = (float)ypos;
    
    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}


