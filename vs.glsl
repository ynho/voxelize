varying vec3 normal;
varying mat4 mvp;

void main (void)
{
  // gl_FrontColor = vec4 (0.5, 0.0, 0.0, 0.0);
  // nor = normalize (gl_NormalMatrix * gl_Normal);
  normal = normalize (gl_Normal);
  mvp = gl_ModelViewProjectionMatrix;
  // gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
  gl_Position = gl_Vertex;
  gl_FrontColor = gl_Color;
}
