varying vec3 nor;

void main (void)
{
    vec3 col1 = vec3(0.8, 0.1, 0.2);
    vec3 col2 = vec3(0.2, 0.7, 0.3);
    vec3 col3 = vec3(0.1, 0.3, 0.9);
    gl_FragColor.xyz = col1 * clamp (dot (normalize (nor),
                                      normalize (vec3 (1.0, 1.0, 0.3))), 0.0, 1.0)
      +  col2 * clamp (dot (normalize (nor),
                            normalize (vec3 (-1.0, 1.0, -0.6))), 0.0, 1.0)
      +  col3 * clamp (dot (normalize (nor),
                            normalize (vec3 (0.0, -1.0, 0.6))), 0.0, 1.0);
}
