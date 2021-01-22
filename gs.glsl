#version 330
#extension GL_EXT_geometry_shader4 : enable

layout (triangle_strip, max_vertices = 3) out;

varying in mat4 mvp[];
varying in vec3 normal[];
varying out vec3 nor;
varying out vec3 col;

void main (void)
{
    vec3 norm = cross (vec3(gl_PositionIn[0] - gl_PositionIn[1]), vec3(gl_PositionIn[0] - gl_PositionIn[2]));
    // vec3 norm = normalize(normal[0] + normal[1] + normal[2]);
    // nor = normal[0];
    nor = norm;
    vec3 colour = gl_FrontColorIn[2].x > 0.5 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
    gl_Position = mvp[0] * gl_PositionIn[0];
    col = colour;
    EmitVertex();
    gl_Position = mvp[0] * gl_PositionIn[1];
    col = colour;
    EmitVertex();
    gl_Position = mvp[0] * gl_PositionIn[2];
    col = colour;
    EmitVertex();

    EndPrimitive ();
}
