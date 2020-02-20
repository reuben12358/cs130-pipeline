#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];
    std::fill(state.image_color, state.image_color + (width * height), make_pixel(0,0,0));
    std::fill(state.image_depth, state.image_depth + (width * height), 2);
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    //std::cout<<"TODO: implement rendering."<<std::endl;
    for (int i = 0; i < state.num_vertices; i += 3) {
        data_geometry* pass[3];
        for (int j = 0; j < 3; ++j) {
            data_geometry* temp_g = new data_geometry();
            temp_g -> data = state.vertex_data + (i + j)*state.floats_per_vertex;
            data_vertex temp_v;
            temp_v.data = temp_g -> data;
            state.vertex_shader(temp_v, *temp_g, state.uniform_data);
            pass[j] = temp_g; 
        }
        rasterize_triangle(state, pass);
    }
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, data_geometry* in[3])
{
    // std::cout<<"TODO: implement rasterization"<<std::endl;
    for (int i = 0; i < 3; ++i) {
        data_vertex vertex;
        data_geometry* k = in[i];
        vertex.data = in[i] -> data;
        state.vertex_shader(vertex, *k, state.uniform_data);
        for (int j = 0; j < 3; ++j) {
            in[i] -> gl_Position[j] = in[i] -> gl_Position[j] / in[i] -> gl_Position[3];
        }
        in[i] -> gl_Position[0] = (in[i] -> gl_Position[0] + 1) * (state.image_width / 2);
        in[i] -> gl_Position[1] = (in[i] -> gl_Position[1] + 1) * (state.image_height / 2);
        int temp = in[i]->gl_Position[1] * state.image_width + in[i]->gl_Position[0];
        state.image_color[temp] = make_pixel(255, 255, 255);
    }

    float area = ((in[0]->gl_Position[0] * (in[1]->gl_Position[1] - in[2]->gl_Position[1])) + 
		        (in[1]->gl_Position[0] * (in[2]->gl_Position[1] - in[0]->gl_Position[1])) + 
		        (in[2]->gl_Position[0] * (in[0]->gl_Position[1] - in[1]->gl_Position[1]))) / 2;

    int minX = std::min(std::min(in[0]->gl_Position[0], in[1]->gl_Position[0]), in[2]->gl_Position[0]);
    int maxX = std::max(std::max(in[0]->gl_Position[0], in[1]->gl_Position[0]), in[2]->gl_Position[0]);
    int minY = std::min(std::min(in[0]->gl_Position[1], in[1]->gl_Position[1]), in[2]->gl_Position[1]);
    int maxY = std::max(std::max(in[0]->gl_Position[1], in[1]->gl_Position[1]), in[2]->gl_Position[1]);

    for (int i = minX; i < maxX; ++i) {
        for (int j = minY; j < maxY; ++j) {
            float b = ((in[0]->gl_Position[0] * (j - in[2]->gl_Position[1])) + 
          			  (i * (in[2]->gl_Position[1] - in[0]->gl_Position[1])) + 
			          (in[2]->gl_Position[0] * (in[0]->gl_Position[1] - j))) / 2;

            float y = ((in[0]->gl_Position[0] * (in[1]->gl_Position[1] - j)) + 
			          (in[1]->gl_Position[0] * (j - in[0]->gl_Position[1])) + 
          			  (i * (in[0]->gl_Position[1] - in[1]->gl_Position[1]))) / 2;
            double beta = b/area;
            double gamma = y/area;
            double alpha = 1 - beta - gamma;
            if (alpha >= 0 && beta >= 0 && gamma >= 0 && (alpha + beta + gamma) <= 1.001) {
                state.image_color[(j * state.image_width) + i] = make_pixel(255, 255, 255);
            }
        }
    }
}

