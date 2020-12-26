#version 330
#extension GL_EXT_geometry_shader4 : enable

layout (triangle_strip, max_vertices = 3) out;

varying in mat4 mvp[];
varying in vec3 normal[];
varying out vec3 nor;

void main (void)
{
    vec3 norm = cross (vec3(gl_PositionIn[0] - gl_PositionIn[1]), vec3(gl_PositionIn[0] - gl_PositionIn[2]));
    // vec3 norm = normalize(normal[0] + normal[1] + normal[2]);
  // nor = normal[0];
    nor = norm;
  gl_Position = mvp[0] * gl_PositionIn[0];
  EmitVertex();
  gl_Position = mvp[0] * gl_PositionIn[1];
  EmitVertex();
  gl_Position = mvp[0] * gl_PositionIn[2];
  EmitVertex();

  EndPrimitive ();
}
