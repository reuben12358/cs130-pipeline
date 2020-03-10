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
    if (type == render_type::triangle) {
        for (int i = 0; i < state.num_vertices; i = i + 3) {
            const data_geometry* pass[3];
            for (int j = 0; j < 3; ++j) {
                data_geometry* dg = new data_geometry();
                dg->data = state.vertex_data + (i + j)*state.floats_per_vertex;
                data_vertex dv;
                dv.data = dg->data;
                state.vertex_shader(dv, *dg, state.uniform_data);
                pass[j] = dg;
            }
            clip_triangle(state, pass, 0);
            for (int j = 0; j < 3; ++j) {
                delete pass[j]; 
            }
        }
    }
    else if (type == render_type::indexed) {
        for (int i = 0; i < 3 * state.num_triangles; i = i + 3) {
            const data_geometry* pass[3];
            for (int j = 0; j < 3; ++j) {
                data_geometry* dg = new data_geometry();
                dg->data = state.vertex_data + state.index_data[i + j] * state.floats_per_vertex;
                data_vertex dv;
                dv.data = dg->data;
                state.vertex_shader(dv, *dg, state.uniform_data);
                pass[j] = dg;
            }
            clip_triangle(state, pass, 0);
            for (int j = 0; j < 3; ++j) {
                delete pass[j]; 
            }
        }
    }
    else if (type == render_type::fan) { //First vertex is same for all triangles
        for (int i = 0; i < state.num_vertices; ++i) {
            const data_geometry* pass[3];
            for (int j = 0; j < 3; ++j) {
                data_geometry* dg = new data_geometry();
                if (j == 0) {
                    dg->data = state.vertex_data;
                }
                else {
                    dg->data = state.vertex_data + (i + j) * state.floats_per_vertex;
                }
                data_vertex dv;
                dv.data = dg->data;
                state.vertex_shader(dv, *dg, state.uniform_data);
                pass[j] = dg;
            }
            clip_triangle(state, pass, 0);
            for (int j = 0; j < 3; ++j) {
                delete pass[j]; 
            }
        }
    }
    else if (type == render_type::strip) { //Two vertices are shared between triangles
        for (int i = 0; i < (state.num_vertices - 2); i++) {
            const data_geometry* pass[3];
            for (int j = 0; j < 3; j++) {
                data_geometry* dg = new data_geometry();
                dg->data = state.vertex_data + (i + j) * state.floats_per_vertex;
                data_vertex dv;
                dv.data = dg->data;
                state.vertex_shader(dv, *dg, state.uniform_data);
                pass[j] = dg;
            }
            clip_triangle(state, pass, 0);
            for (int j = 0; j < 3; j++) {
                delete pass[j]; 
            }
        }
    }
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6) {
        rasterize_triangle(state, in);
        return;
    }

    float cnt = 0, ind = 0, val = 1;
    ivec3 inVec = {0, 0, 0};
    if (face == 0) { 		// x = 1
        ind = 0;
        val = 1;
    } 
    else if (face == 1) { // x = -1
        ind = 0;
        val = -1;
    } 
    else if (face == 2) { // y = 1
        ind = 1;
        val = 1;
    } 
    else if (face == 3) { // y = -1
        ind = 1;
        val = -1;
    } 
    else if (face == 4) { // z = 1
        ind = 2;
        val = 1;
    } 
    else if (face == 5) { // z = -1
        ind = 2;
        val = -1;
    }

    //If vertex is inside, set it's inVec value to 1 and increment a counter (inVec stores which vertex is inside, and the counter tells how many vertices are inside)
    for (int i = 0; i < 3; i++) {
        if (val < 0 && in[i]->gl_Position[ind] >= (val * in[i]->gl_Position[3])) {
            cnt++;
            inVec[i] = 1;
        } 
        else if (val > 0 && in[i]->gl_Position[ind] <= (val * in[i]->gl_Position[3])) {
            cnt++;
            inVec[i] = 1;
        }
  }

    if (cnt == 0) {
        return;
    }
    else if (cnt == 1) {
        const data_geometry* A = 0;
        const data_geometry* B = 0;
        const data_geometry* C = 0;

        for (int i = 0; i < 3; i++) {
            if (inVec[i] == 1) {
                A = in[i];
                B = in[(i + 1) % 3];
                C = in[(i + 2) % 3];
                break;
            }
        }

        float wA = A->gl_Position[3], wB = B->gl_Position[3], wC = C->gl_Position[3];
        float aV = A->gl_Position[ind], bV = B->gl_Position[ind], cV = C->gl_Position[ind];
        float alphaB = ((val * wB) - bV) / ((aV - (val * wA)) + ((val * wB) - bV));
        float alphaC = ((val * wC) - cV) / ((aV - (val * wA)) + ((val * wC) - cV));

        data_geometry* AB = new data_geometry();
        data_geometry* AC = new data_geometry();
    
        AB->data = new float[state.floats_per_vertex];
        AC->data = new float[state.floats_per_vertex];
        for (int i = 0; i < state.floats_per_vertex; i++) {
            if (state.interp_rules[i] == interp_type::flat) {
                AB->data[i] = A->data[i];
                AC->data[i] = A->data[i];
            } 
            else if (state.interp_rules[i] == interp_type::smooth) {
                AB->data[i] = (alphaB * A->data[i]) + ((1 - alphaB) * B->data[i]);
                AC->data[i] = (alphaC * A->data[i]) + ((1 - alphaC) * C->data[i]);
            } 
            else if (state.interp_rules[i] == interp_type::noperspective) {
                float k = (alphaB * A->gl_Position[3]) + ((1 - alphaB) * B->gl_Position[3]);
                float alpha = (alphaB * A->gl_Position[3]) / k;
                AB->data[i] = (alpha * A->data[i]) + ((1 - alpha) * B->data[i]);
                k = (alphaC * A->gl_Position[3]) + ((1 - alphaC) * C->gl_Position[3]);
                alpha = (alphaC * A->gl_Position[3]) / k;
                AC->data[i] = (alpha * A->data[i]) + ((1 - alpha) * C->data[i]);
            }
        }
        for (int i = 0; i < 4; i++) {
            AB->gl_Position[i] = (alphaB * A->gl_Position[i]) + ((1 - alphaB) * B->gl_Position[i]);
            AC->gl_Position[i] = (alphaC * A->gl_Position[i]) + ((1 - alphaC) * C->gl_Position[i]);
        }

        const data_geometry* pass[3];
        pass[0] = A;
        pass[1] = AB;
        pass[2] = AC;
        clip_triangle(state,pass,face+1);
        delete AB;
        delete AC;
        return;
    }
    else if (cnt == 2) {
        const data_geometry* A = 0;
        const data_geometry* B = 0;
        const data_geometry* C = 0;

        for (int i = 0; i < 3; i++) {
            if (inVec[i] == 0) {
                A = in[i];
                B = in[(i + 1) % 3];
                C = in[(i + 2) % 3];
                break;
            }
        }

        float wA = A->gl_Position[3], wB = B->gl_Position[3], wC = C->gl_Position[3];
        float aV = A->gl_Position[ind], bV = B->gl_Position[ind], cV = C->gl_Position[ind];
        float alphaB = ((val * wB) - bV) / ((aV - (val * wA)) + ((val * wB) - bV));
        float alphaC = ((val * wC) - cV) / ((aV - (val * wA)) + ((val * wC) - cV));

        data_geometry* AB = new data_geometry();
        data_geometry* AC = new data_geometry();
    
        AB->data = new float[state.floats_per_vertex];
        AC->data = new float[state.floats_per_vertex];
        for (int i = 0; i < state.floats_per_vertex; i++) {
            if (state.interp_rules[i] == interp_type::flat) {
                AB->data[i] = A->data[i];
                AC->data[i] = A->data[i];
            } 
            else if (state.interp_rules[i] == interp_type::smooth) {
                AB->data[i] = (alphaB * A->data[i]) + ((1 - alphaB) * B->data[i]);
                AC->data[i] = (alphaC * A->data[i]) + ((1 - alphaC) * C->data[i]);
            } 
            else if (state.interp_rules[i] == interp_type::noperspective) {
                float k = (alphaB * A->gl_Position[3]) + ((1 - alphaB) * B->gl_Position[3]);
                float alpha = (alphaB * A->gl_Position[3]) / k;
                AB->data[i] = (alpha * A->data[i]) + ((1 - alpha) * B->data[i]);
                k = (alphaC * A->gl_Position[3]) + ((1 - alphaC) * C->gl_Position[3]);
                alpha = (alphaC * A->gl_Position[3]) / k;
                AC->data[i] = (alpha * A->data[i]) + ((1 - alpha) * C->data[i]);
            }
        }
        for (int i = 0; i < 4; i++) {
            AB->gl_Position[i] = (alphaB * A->gl_Position[i]) + ((1 - alphaB) * B->gl_Position[i]);
            AC->gl_Position[i] = (alphaC * A->gl_Position[i]) + ((1 - alphaC) * C->gl_Position[i]);
        }
        const data_geometry* pass[3];
        pass[0] = AB;
        pass[1] = B;
        pass[2] = C;
        clip_triangle(state,pass,face+1);
        pass[0] = C;
        pass[1] = AC;
        pass[2] = AB;
        clip_triangle(state,pass,face+1);
        delete AB;
        delete AC;
        return;
    }    
    else if (cnt == 3) {
        clip_triangle(state,in,face+1);
        return;
    } 
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    // std::cout<<"TODO: implement rasterization"<<std::endl;
    vec3 a,b,c;

    /*ax*/ a[0] = (((in[0] -> gl_Position[0] / in[0] -> gl_Position[3]) + 1) * (state.image_width/2)) - 0.5;
    /*ay*/ a[1] = (((in[0] -> gl_Position[1] / in[0] -> gl_Position[3]) + 1) * (state.image_height/2)) - 0.5;
    /*az*/ a[2] = in[0] -> gl_Position[2] / in[0] -> gl_Position[3];
    /*bx*/ b[0] = (((in[1] -> gl_Position[0] / in[1] -> gl_Position[3]) + 1) * (state.image_width/2)) - 0.5;
    /*by*/ b[1] = (((in[1] -> gl_Position[1] / in[1] -> gl_Position[3]) + 1) * (state.image_height/2)) - 0.5;
    /*bz*/ b[2] = in[1] -> gl_Position[2] / in[1] -> gl_Position[3];
    /*cx*/ c[0] = (((in[2] -> gl_Position[0] / in[2] -> gl_Position[3]) + 1) * (state.image_width/2)) - 0.5;
    /*cy*/ c[1] = (((in[2] -> gl_Position[1] / in[2] -> gl_Position[3]) + 1) * (state.image_height/2)) - 0.5;
    /*cz*/ c[2] = in[2] -> gl_Position[2] / in[2] -> gl_Position[3];

    float area = ((a[0] * (b[1] - c[1])) + (b[0] * (c[1] - a[1])) + (c[0] * (a[1] - b[1]))) / 2;

    int minX = std::min(std::min(a[0], b[0]), c[0]);
    int maxX = std::max(std::max(a[0], b[0]), c[0]);
    int minY = std::min(std::min(a[1], b[1]), c[1]);
    int maxY = std::max(std::max(a[1], b[1]), c[1]);

    for (int i = minX; i <= maxX; ++i) {
        for (int j = minY; j <= maxY; ++j) {
            float bt = ((a[0] * ((float)j - c[1])) + 
          			  ((float)i * (c[1] - a[1])) + 
			          (c[0] * (a[1] - (float)j))) / 2;

            float y = ((a[0] * (b[1] - (float)j)) + 
			          (b[0] * ((float)j - a[1])) + 
          			  ((float)i * (a[1] - b[1]))) / 2;

            float al = (((float)i * (b[1] - c[1])) +
                      (b[0] * (c[1] - (float)j)) +
                      (c[0] * ((float)j - b[1]))) / 2;

            double beta = bt/area;
            double gamma = y/area;
            double alpha = al/area;
            if (alpha >= 0 && beta >= 0 && gamma >= 0 && (alpha + beta + gamma) <= 1.001) {
                // implement color for test 6 onwards
                float color[state.floats_per_vertex];
                // calculate co-ords using method from shroeders lecture
                for (int i = 0; i < state.floats_per_vertex; ++i) {
                    if (state.interp_rules[i] == interp_type::flat) {
                        color[i] = in[0] -> data[i];
                    }
                    else if (state.interp_rules[i] == interp_type::noperspective) {
                        color[i] = (alpha * in[0] -> data[i]) + (beta * in[1] -> data[i]) + (gamma * in[2] -> data[i]);
                    }
                    else if (state.interp_rules[i] == interp_type::smooth) {
                        float k = alpha / in[0] -> gl_Position[3] + beta / in[1] -> gl_Position[3] + gamma / in[2] -> gl_Position[3];
                        color[i] = ((alpha * in[0] -> data[i] / in[0] -> gl_Position[3]) 
                        + (beta * in[1] -> data[i] / in[1] -> gl_Position[3])
                        + (gamma * in[2] -> data[i] / in[2] -> gl_Position[3])) / k;
                    }
                    else { //invalid
                        color[i] = 1;
                    }
                }
                if (alpha * a[2] + beta * b[2] + gamma * c[2] < state.image_depth[(j * state.image_width) + i]) {
                    data_fragment temp_df;
                    temp_df.data = color;
                    data_output d;
                    state.fragment_shader(temp_df, d, state.uniform_data);
                    state.image_color[(j * state.image_width) + i] = make_pixel(255*d.output_color[0], 255*d.output_color[1], 255*d.output_color[2]);
                    state.image_depth[(j * state.image_width) + i] = alpha * a[2] + beta * b[2] + gamma * c[2];
                }
            }
        }
    }
}

