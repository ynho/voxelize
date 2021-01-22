#include <strings.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <SDL2/SDL.h>
/* #include <SDL2/SDL_opengl.h> */
#include <SCE/core/SCECore.h>
/* #include <GL/gl.h> */
#include <GL/glew.h>


struct mesh {
    int *vertices;
    SCEvertices *verticesf;
    SCEvertices *normals;
    SCE_TDCMaterial *materials;
    mdcindices_t *indices;
    int *v_flags;
    unsigned int n_vertices;
    unsigned int n_indices;
    int v_edgestart;
    int i_edgestart;
    int mode;
};

static void init_mesh (struct mesh *mesh) {
    mesh->vertices = NULL;
    mesh->verticesf = NULL;
    mesh->normals = NULL;
    mesh->materials = NULL;
    mesh->indices = NULL;
    mesh->v_flags = NULL;
    mesh->n_vertices = 0;
    mesh->n_indices = 0;
    mesh->v_edgestart = 0;
    mesh->i_edgestart = 0;
    mesh->mode = GL_QUADS;
}

static void duplicate_mesh (struct mesh *new, struct mesh *old) {
    size_t vsize = old->n_vertices * 3 * sizeof *new->verticesf;
    size_t isize = old->n_indices * sizeof *new->indices;
    new->verticesf = malloc (vsize);
    new->indices = malloc (isize);
    memcpy (new->verticesf, old->verticesf, vsize);
    if (old->materials) {
        new->materials = malloc (vsize);
        memcpy (new->materials, old->materials, vsize);
    }
    memcpy (new->indices, old->indices, isize);
    new->n_vertices = old->n_vertices;
    new->n_indices = old->n_indices;
    new->mode = old->mode;
}

static void proj (float m[16], float a, float r, float n, float f) {
    m[5] = 1.0f / tanf (a * 0.5f);
    m[0] = m[5] / r;
    m[10] = -f / (f - n * 2.0f);
    m[11] = -2.0f * n * (f / (f - n));
    m[14] = -1.0f;

    m[1] = m[2] = m[3] = m[4] = m[6] = m[7] =
        m[8] = m[9] = m[12] = m[13] = m[15] = 0.0f;
}

static void transpose (float m[16]) {
    float t;
    t = m[1]; m[1] = m[4]; m[4] = t;
    t = m[2]; m[2] = m[8]; m[8] = t;
    t = m[3]; m[3] = m[12]; m[12] = t;
    t = m[6]; m[6] = m[9]; m[9] = t;
    t = m[7]; m[7] = m[13]; m[13] = t;
    t = m[11]; m[11] = m[14]; m[14] = t;
}

#define SCREEN_W 2600
#define SCREEN_H 1900

#define VOX_W 160
#define VOX_H 160
#define VOX_D 80
#define NUM_VOX (VOX_W*VOX_H*VOX_D)

#define SIMP 1

/* #define RAD (0.0174532925) */
static void setup_view (int rx, int ry, int vox_w, int vox_h, int dist) {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef (0.0, 0.0, -dist);
    glRotatef (rx, 1.0, 0.0, 0.0);
    glRotatef (ry, 0.0, 0.0, 1.0);
    glTranslatef (-vox_w / 2.0, -vox_h / 2.0, 0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float matrix[16];
    proj (matrix, 70.0 * RAD, (float)SCREEN_W / SCREEN_H, 0.1, 1000.0);
    transpose (matrix);
    glLoadMatrixf (matrix);
}


static void draw (struct mesh *mesh, int use_float) {
    int coul[6] = {0x00550000,
                   0x00555000,
                   0x00005500,
                   0x00005550,
                   0x00000055,
                   0x00500055};
    /* glPolygonMode(GL_FRONT, GL_LINE); */
    glBegin (mesh->mode);
    for (uint32_t i = 0; i < mesh->n_indices; i++) {
        mdcindices_t k = mesh->indices[i];
        int c = 0x00222222;
        if (mesh->v_flags) {
            for (int j = 0; j < 6; j++) {
                if (mesh->v_flags[k] & (1 << j)) {
                    c += coul[j];
                }
            }
        }
        /* glColor3ub (c >> 16, c >> 8 & 0xFF, c & 0xFF); */
        /* float *nor = &mesh->normals[k * 3]; */
        /* glNormal3f (nor[0], nor[1], nor[2]); */
        if (mesh->materials)
            glColor3ub (mesh->materials[k], 0, 0);
        if (use_float) {
            float *vertex = &mesh->verticesf[k * 3];
            glVertex3f (vertex[0], vertex[1], vertex[2]);
        } else {
            int *vertex = &mesh->vertices[k * 3];
            glVertex3f (vertex[0], vertex[1], vertex[2]);
        }
    }
    glEnd();
}

#if 0
void print_binary(unsigned long number)
{
    if (number >> 1) {
        print_binary(number >> 1);
    }
    putc((number & 1) ? '1' : '0', stdout);
}

static void printb (unsigned long n) {
    int i;
    for (i = 11; i > -1; i--)
        printf ("\t%d", i);
    printf ("\n");
    for (i = 11; i > -1; i--) {
        printf ("\t");
        putc ((n >> (i*2 + 1) & 1) ? '1' : '0', stdout);
        putc ((n >> (i*2 + 0) & 1) ? '1' : '0', stdout);
    }
}
#endif

#define OFFSET(x, y, z, w, h) ((w) * ((z) * (h) + (y)) + (x))

#define N_OCTAVES 5

static const int n_octaves = N_OCTAVES;
#define OC_VOXELS0 256
#define OC_VOXELS1 79
#define OC_VOXELS2 33
#define OC_VOXELS3 12
#define OC_VOXELS4 3
static const double octave_freq[N_OCTAVES] = {
    1.0 / OC_VOXELS0,
    1.0 / OC_VOXELS1,
    1.0 / OC_VOXELS2,
    1.0 / OC_VOXELS3,
    1.0 / OC_VOXELS4
};
#define COEF0 38.0
#define COEF1 17.0
#define COEF2 12.0
#define COEF3 3.9
#define COEF4 0.17
static const double octave_coef[N_OCTAVES] = {
    COEF0, COEF1, COEF2, COEF3, COEF4
};

static const double ground_level = 10.0;

static double density (long x, long y, long z)
{
    double p[3], p2[3];
    int i;
    /* z -= 20; */
    double derp = (-z + ground_level) * 0.725;
    return derp;
    /* x += 200; */
    p[0] = x; p[1] = y; p[2] = z;
    SCE_Vector3_Copy (p2, p);
    SCE_Vector3_Operator1 (p2, *=, 0.003);
    float koi = SCE_Noise_Smooth3D (p2) * 72.0;
    SCE_Vector3_Operator1 (p, +=, koi);
    
    /* p[0] += SCE_Noise_Smooth3D (p2) * 8.0; */
    /* p2[0] += 10.0; */
    /* p2[1] += 7.0; */
    /* p[1] += SCE_Noise_Smooth3D (p2) * 8.0; */
    /* SCE_Vector3_Operator1 (p2, +=, 20.0); */
    /* p[2] += SCE_Noise_Smooth3D (p2) * 8.0; */

    int num_octaves = 5;

    for (i = 0; i < num_octaves; i++) {
        SCE_Vector3_Operator1v (p2, = octave_freq[i] *, p);
        derp += SCE_Noise_Smooth3D (p2) * octave_coef[i];
    }
    return derp;
}

static double densitym (long x, long y, long z)
{
    double p[3], p2[3];
    int i;
    double derp = 0.0;

    /* x += 200; */
    p[0] = x; p[1] = y; p[2] = z;
    SCE_Vector3_Copy (p2, p);
    SCE_Vector3_Operator1 (p2, *=, 0.003);
    float koi = SCE_Noise_Smooth3D (p2) * 72.0;
    SCE_Vector3_Operator1 (p, +=, koi);
    
    /* p[0] += SCE_Noise_Smooth3D (p2) * 8.0; */
    /* p2[0] += 10.0; */
    /* p2[1] += 7.0; */
    /* p[1] += SCE_Noise_Smooth3D (p2) * 8.0; */
    /* SCE_Vector3_Operator1 (p2, +=, 20.0); */
    /* p[2] += SCE_Noise_Smooth3D (p2) * 8.0; */

    int num_octaves = 5;

    for (i = 0; i < num_octaves; i++) {
        SCE_Vector3_Operator1v (p2, = octave_freq[i] *, p);
        derp += SCE_Noise_Smooth3D (p2) * octave_coef[i];
    }
    return derp;
}

static unsigned char voxel_density (long x, long y, long z) {
    float d = (double)SCE_Math_Clampf (density(x, y, z), -1.0, 1.0);
    return (d * 0.5 + 0.5) * 255;
}
static unsigned char voxelm_density (long x, long y, long z) {
    float d = (double)SCE_Math_Clampf (densitym(x, y, z), -1.0, 1.0);
    return (d * 0.5 + 0.5) * 255;
}

static void generate_voxels (uint8_t *voxels, SCE_TDCMaterial *materials) {
    for (int z = 0; z < VOX_D; z++) {
        for (int y = 0; y < VOX_H; y++) {
            for (int x = 0; x < VOX_W; x++) {
                size_t offset = OFFSET(x, y, z, VOX_W, VOX_H);
                voxels[offset] = voxel_density (x, y, z);
                materials[offset] = voxelm_density (x*4+100, y*3+100, z-30) > 127 ? 200 : 10;
            }
        }
    }
}


static void compute_gradient (SCE_TVector3 gradient, uint8_t *densities,
                              int i, int j, int k, SCEulong w, SCEulong h)
{
    gradient[0] = densities[OFFSET (i-1, j, k, w, h)] - densities[OFFSET (i+1, j, k, w, h)];
    gradient[1] = densities[OFFSET (i, j-1, k, w, h)] - densities[OFFSET (i, j+1, k, w, h)];
    gradient[2] = densities[OFFSET (i, j, k-1, w, h)] - densities[OFFSET (i, j, k+1, w, h)];
}

static float middle (uint8_t a, uint8_t b) {
#if 0
    return (127.0 - (float)a) / ((float)b - a);
#else
    float w = (127.0 - (float)a) / ((float)b - a);
    if (w < 0.0)
        return SCE_Math_Clampf (1.0 - w, 0.0, 1.0);
    else
        return SCE_Math_Clampf (w, 0.0, 1.0);
#endif
}

void generate_hermite (uint8_t *voxels, uint8_t *bitmap, SCE_SDCHermiteData *hermite) {
    for (int i = 0; i < NUM_VOX; i++) {
        bitmap[i / 8] |= (voxels[i] > 127 ? 1 : 0) << (i % 8);
        for (int j = 0; j < 3; j++)
            SCE_Vector4_Set (hermite[i].normalw[j], 0.0, 0.0, 1.0, 0.5);
    }

    float *gradients = malloc (3 * NUM_VOX * sizeof *gradients);
    for (int i = 0; i < NUM_VOX; i++)
        SCE_Vector3_Set (&gradients[i * 3], 0.0, 0.0, 0.0);

    for (int z = 1; z < VOX_D - 1; z++) {
        for (int y = 1; y < VOX_H - 1; y++) {
            for (int x = 1; x < VOX_W - 1; x++) {
                float w;
                size_t offset = OFFSET(x, y, z, VOX_W, VOX_H);
                w = middle (voxels[offset], voxels[OFFSET(x+1, y, z, VOX_W, VOX_H)]);
                hermite[offset].normalw[0][3] = w;
                w = middle (voxels[offset], voxels[OFFSET(x, y+1, z, VOX_W, VOX_H)]);
                hermite[offset].normalw[1][3] = w;
                w = middle (voxels[offset], voxels[OFFSET(x, y, z+1, VOX_W, VOX_H)]);
                hermite[offset].normalw[2][3] = w;
                compute_gradient (&gradients[OFFSET(x, y, z, VOX_W, VOX_H) * 3], voxels, x, y, z, VOX_W, VOX_H);
                SCE_Vector3_Operator1 (&gradients[OFFSET(x, y, z, VOX_W, VOX_H) * 3], *=, 0.5);
            }
        }
    }
    for (int z = 1; z < VOX_D - 2; z++) {
        for (int y = 1; y < VOX_H - 2; y++) {
            for (int x = 1; x < VOX_W - 2; x++) {
                float *o = &gradients[OFFSET(x, y, z, VOX_W, VOX_H) * 3];
                float *ox = &gradients[OFFSET(x+1, y, z, VOX_W, VOX_H) * 3];
                float *oy = &gradients[OFFSET(x, y+1, z, VOX_W, VOX_H) * 3];
                float *oz = &gradients[OFFSET(x, y, z+1, VOX_W, VOX_H) * 3];
                SCE_SDCHermiteData *h = &hermite[OFFSET(x, y, z, VOX_W, VOX_H)];
                SCE_Vector3_Operator2v (h->normalw[0], =, o, +, ox);

                if (SCE_Vector3_Length (h->normalw[0]) > 0.0000001)
                    SCE_Vector3_Normalize (h->normalw[0]);
                SCE_Vector3_Operator2v (h->normalw[1], =, o, +, oy);
                if (SCE_Vector3_Length (h->normalw[1]) > 0.0000001)
                    SCE_Vector3_Normalize (h->normalw[1]);
                SCE_Vector3_Operator2v (h->normalw[2], =, o, +, oz);
                if (SCE_Vector3_Length (h->normalw[2]) > 0.0000001)
                    SCE_Vector3_Normalize (h->normalw[2]);
            }
        }
    }
    free (gradients);
#if 0
    memset (bitmap, 0, NUM_VOX);
    /* bitmap[2] = 1 << 5 | 1 << 6; */
    /* bitmap[5] = 1 << 1 | 1 << 2; */
    bitmap[2] = 1 << 5;
    bitmap[3] = 1 << 2;
    bitmap[5] = 1 << 1;
    bitmap[4] = 1 << 6;
    /* SCE_Vector4_Set (hermite[OFFSET(1, 0, 1, 4, 4)].normalw[1], 0.0, -1.0, 0.0, 0.5 */
#endif
}

uint32_t readint (void) {
    uint32_t n;
    fread (&n, 1, sizeof n, stdin);
    return n;
}

static void load_hermite (uint8_t **bitmap, SCE_SDCHermiteData **hermite,
                          uint32_t *vox_w, uint32_t *vox_h, uint32_t *vox_d) {
    *vox_w = readint ();
    *vox_h = readint ();
    *vox_d = readint ();
    int n_vox = *vox_w * *vox_h * *vox_d;
    *bitmap = malloc (n_vox / 8 + 1);
    fread (*bitmap, 1, n_vox / 8 + 1, stdin);
    *hermite = malloc (n_vox * sizeof **hermite);
    fread (*hermite, sizeof **hermite, n_vox, stdin);
}


char* load_file (const char *fname) {
    char *data = NULL;
    FILE *fp = fopen(fname, "r");
    long size;
    fseek (fp, 0, SEEK_END);
    size = ftell (fp);
    fseek (fp, 0, SEEK_SET);
    data = malloc (size + 1);
    fread (data, 1, size, fp);
    fclose (fp);
    data[size] = 0;
    return data;
}

int mk_shader (const char *fname, int type) {
    int s = glCreateShader (type);
    char *src = load_file (fname);
    glShaderSource (s, 1, (const GLchar * const*)&src, NULL);
    glCompileShader (s);
    int st;
    glGetShaderiv (s, GL_COMPILE_STATUS, &st);
    if (st != GL_TRUE) {
        char buf[512] = {0};
        int size = 512;
        glGetShaderInfoLog (s, size, &size, buf);
        fprintf (stderr, "no compilo %s:\n%s\n", fname, buf);
        exit (42);
    }
    return s;
}

int load_shader (void) {
    int types[] = {GL_VERTEX_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER};
    char *src[] = {"vs.glsl", "gs.glsl", "ps.glsl"};
    /* int types[] = {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER}; */
    /* char *src[] = {"vs.glsl", "ps.glsl"}; */

    int p = glCreateProgram ();
    for (int i = 0; i < sizeof types / sizeof *types; i++) {
        int s = mk_shader (src[i], types[i]);
        glAttachShader (p, s);
    }
    glLinkProgram (p);
    int st;
    glGetProgramiv (p, GL_LINK_STATUS, &st);
    if (st != GL_TRUE) {
        char buf[512] = {0};
        int size = 512;
        glGetProgramInfoLog (p, size, &size, buf);
        fprintf (stderr, "no linko:\n%s\n", buf);
        exit (42);
    }
    return p;
}

static int get_config (uint8_t *bitmap, int x, int y, int z, int w, int h) {
    int config = 0;
    size_t c;
#define FETCH_BIT(dx, dy, dz, s)                                        \
    (c = OFFSET(x + dx, y + dy, z + dz, w, h), (bitmap[c / 8] >> c % 8 & 1) << s)
    config |= FETCH_BIT(0, 0, 0, 0);
    config |= FETCH_BIT(1, 0, 0, 1);
    config |= FETCH_BIT(0, 1, 0, 2);
    config |= FETCH_BIT(1, 1, 0, 3);
    config |= FETCH_BIT(0, 0, 1, 4);
    config |= FETCH_BIT(1, 0, 1, 5);
    config |= FETCH_BIT(0, 1, 1, 6);
    config |= FETCH_BIT(1, 1, 1, 7);
    return config;
}

static int zero(float a) { return SCE_Math_IsZero (a); }

static void draw_normal (SCE_SDCHermiteData *hermite, int x, int y, int z, size_t offset, int axis) {
    SCE_TVector3 p1, p2;
    SCE_Vector3_Set (p1, x, y, z);
    p1[axis] += hermite[offset].normalw[axis][3];
    SCE_Vector3_Operator2v (p2, =, p1, +, hermite[offset].normalw[axis]);

    glBegin (GL_LINES);
    glColor3f (1.0, 1.0, 1.0);
    glVertex3f (p1[0], p1[1], p1[2]);
    glColor3f (1.0, 0.0, 0.0);
    glVertex3f (p2[0], p2[1], p2[2]);
    glEnd ();
}

static void draw_hermite (uint8_t *bitmap, SCE_SDCHermiteData *hermite,
                          uint32_t vox_w, uint32_t vox_h, uint32_t vox_d) {

    for (int z = 0; z < vox_d - 1; z++) {
        for (int y = 0; y < vox_h - 1; y++) {
            for (int x = 0; x < vox_w - 1; x++) {
                int config = get_config (bitmap, x, y, z, vox_w, vox_h);
                if (config != 0 && config != 255) {
                    size_t offset = OFFSET(x, y, z, vox_w, vox_h);
                    if ((config & 1) != ((config >> 1) & 1))
                        draw_normal (hermite, x, y, z, offset, 0);
                    if ((config & 1) != ((config >> 2) & 1))
                        draw_normal (hermite, x, y, z, offset, 1);
                    if ((config & 1) != ((config >> 4) & 1))
                        draw_normal (hermite, x, y, z, offset, 2);
                }
            }
        }
    }
}


static void triangulate (struct mesh *mesh) {
    size_t num = 6 * mesh->n_indices / 4;
    uint32_t *indices = malloc (num * sizeof *indices);
    for (int i = 0, j = 0; i < mesh->n_indices; i += 4, j += 6) {
        indices[j] = mesh->indices[i];
        indices[j+1] = mesh->indices[i+1];
        indices[j+2] = mesh->indices[i+3];
        indices[j+3] = mesh->indices[i+1];
        indices[j+4] = mesh->indices[i+2];
        indices[j+5] = mesh->indices[i+3];
    }
    free(mesh->indices);
    mesh->indices=indices;
    mesh->n_indices = num;
    mesh->mode = GL_TRIANGLES;
}

int main (void) {
    SDL_Window *Window = NULL;
    SDL_GLContext glContext;
    const int ww = SCREEN_W, wh = SCREEN_H;
    uint32_t vox_w, vox_h, vox_d = VOX_D;

    SDL_Init (SDL_INIT_VIDEO);
    Window = SDL_CreateWindow ("LE MAO", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                               ww, wh, SDL_WINDOW_OPENGL);
    glContext = SDL_GL_CreateContext (Window);
    (void)glContext;
    SDL_ShowWindow (Window);
    glewInit ();

    glClearColor (0.0f, 0.0f, 0.0f, 0.0f);
    glClear (GL_COLOR_BUFFER_BIT);
    glViewport (0.0, 0.0, ww, wh);
    SDL_Event ev;
    int running = 1;

    srand(1547917722);
    SCE_Init_Core (stderr, 0);
    glEnable(GL_DEPTH_TEST);

    struct mesh mesh;

    init_mesh (&mesh);

    uint8_t *bitmap;
    {
        SCE_SDCHermiteData *hermite;

#if 1
        uint8_t *voxels;
        SCE_TDCMaterial *materials;
        vox_w = VOX_W, vox_h = VOX_H, vox_d = VOX_D;
        hermite = malloc (NUM_VOX * sizeof *hermite);
        bitmap = malloc (NUM_VOX);
        voxels = malloc (NUM_VOX);
        materials = malloc (NUM_VOX * sizeof *materials);
        memset (bitmap, 0, NUM_VOX);
        memset (voxels, 0, NUM_VOX);
        for (int i = 0; i < NUM_VOX; i++)
            materials[i] = SCE_DC_UNDEFINED_MATERIAL;
        generate_voxels (voxels, materials);
        generate_hermite (voxels, bitmap, hermite);
#else
        load_hermite (&bitmap, &hermite, &vox_w, &vox_h, &vox_d);
#endif

#if 0
        /* supposedly manifold but actually fucking isnt */
        SCE_SMDCGenerator gen;
        SCE_MDC_Init (&gen);
        SCE_MDC_Build (&gen, vox_w, vox_h, vox_d);
        mesh.n_vertices = SCE_MDC_ComputeNumVertices (&gen, bitmap);
        mesh.verticesf = malloc (mesh.n_vertices * sizeof(float) * 3);
        mesh.normals = malloc (mesh.n_vertices * sizeof(float) * 3);
        SCE_MDC_GenerateVertices (&gen, bitmap, hermite, mesh.verticesf);
        mesh.n_indices = SCE_MDC_ComputeNumIndices (&gen, bitmap);
        mesh.indices = malloc (mesh.n_indices * sizeof *mesh.indices);
        size_t test_i = SCE_MDC_GenerateIndices (&gen, bitmap, mesh.indices);
#else
        /* classic dual contouring (not manifold but proud not to be) */
        SCE_SDCGenerator gen;
        SCE_DC_Init (&gen);
        SCE_DC_Build (&gen, vox_w, vox_h, vox_d);
        mesh.n_vertices = SCE_DC_ComputeNumVertices (&gen, bitmap);
        mesh.verticesf = malloc (mesh.n_vertices * sizeof(float) * 3);
        /* mesh.normals = malloc (mesh.n_vertices * sizeof(float) * 3); */
        mesh.materials = malloc (mesh.n_vertices * sizeof *mesh.materials);
        SCE_DC_GenerateVertices (&gen, bitmap, hermite, mesh.verticesf);
        mesh.n_indices = SCE_DC_ComputeNumIndices (&gen, bitmap);
        mesh.indices = malloc (mesh.n_indices * sizeof *mesh.indices);
        size_t test_i = SCE_DC_GenerateIndicesAndMaterials (&gen, bitmap, materials, mesh.indices, mesh.materials);
        int lol = 0;
        for (int i = 0; i < mesh.n_vertices; i++) {
            if (mesh.materials[i] == SCE_DC_UNDEFINED_MATERIAL)
                lol++;
        }
        printf ("lol : %d\n");
#endif

        if (mesh.n_indices != test_i) {
            printf ("omg wtf %d %d\n", mesh.n_indices, test_i);
        }
        /* SCE_Geometry_ComputeQuadsNormals (mesh.verticesf, mesh.indices, mesh.n_vertices, */
        /*                                   mesh.n_indices, mesh.normals); */
        /* printf ("\n"); */
        printf ("n_indices = %d\n", mesh.n_indices);
        printf ("n_vertices = %d\n", mesh.n_vertices);
        /* for (int i = 0; i < mesh.n_indices; i++) { */
        /*     printf ("%d ", mesh.indices[i]); */
        /* } */
        /* printf ("\n"); */
        triangulate (&mesh);
    }

#if SIMP
    /* simplify */
    struct mesh simp;
    {
        SCE_SQEMMesh qm;

        init_mesh (&simp);
        duplicate_mesh (&simp, &mesh);
        /* triangulate (&simp); */

        SCE_QEMD_Init (&qm);
        SCE_QEMD_SetMaxVertices (&qm, simp.n_vertices);
        SCE_QEMD_SetMaxIndices (&qm, simp.n_indices);
        SCE_QEMD_Build (&qm);
        /* float *vert_r = malloc (simp.n_vertices * 3 * sizeof *vert_r); */
        /* uint32_t *ind_r = malloc (simp.n_indices * sizeof *ind_r); */
        SCE_QEMD_Set (&qm, simp.verticesf, NULL, simp.materials, NULL, simp.indices, simp.n_vertices, simp.n_indices);
        SCE_QEMD_Process (&qm, simp.n_vertices * 0.97);
        SCE_QEMD_Get (&qm, simp.verticesf, NULL, simp.materials, simp.indices, &simp.n_vertices, &simp.n_indices);
        printf ("\nn_indices = %d\n", simp.n_indices);
        printf ("n_vertices = %d\n", simp.n_vertices);
        /* the decimator conveniently takes care of the unused vertices generated by the MDC */
    }
#endif

    int p = load_shader ();

    glEnable (GL_CULL_FACE);
    /* glCullFace (GL_FRONT); */
    glClearColor (0.5, 0.5, 0.5, 0.0);

    int prev_x = 0, prev_y = 0, ry = 0, rx = 0, mouse_pressed = 0;
    float dist = 4.0;//vox_d;
    int wireframe = 0;
    while (running){
        while (SDL_PollEvent (&ev)) {
            switch (ev.type) {
            case SDL_QUIT:
                running = 0;
                break;
            case SDL_KEYDOWN:
                switch (ev.key.keysym.sym) {
				case SDLK_ESCAPE:
					running = 0;
					break;
                case SDLK_w:
                    wireframe = !wireframe;
                    if (wireframe)
                        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    else
                        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    break;
				default:
					break;
                }
                break;

            case SDL_MOUSEBUTTONDOWN:
                prev_x = ev.button.x;
                prev_y = ev.button.y;
                mouse_pressed = 1;
                break;

            case SDL_MOUSEBUTTONUP:
                mouse_pressed = 0;
                break;

            case SDL_MOUSEWHEEL:
                dist -= ev.wheel.y;
                break;

            case SDL_MOUSEMOTION:
                if (mouse_pressed) {
                    ry += ev.motion.x - prev_x;
                    rx += ev.motion.y - prev_y;
                    prev_x = ev.motion.x;
                    prev_y = ev.motion.y;
                }
                break;
            }
        }
        glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        setup_view (rx, ry, vox_w, vox_h, dist);

        glUseProgram (p);
#if SIMP
        glTranslatef (-100.0, 0.0, 0.0);
        draw (&mesh, 1);
        glTranslatef (200.0, 0.0, 0.0);
        draw (&simp, 1);
#else
        draw (&mesh, 1);
#endif
        /* glUseProgram (0); */
        /* draw_hermite (bitmap, hermite, vox_w, vox_h, vox_d); */

        SDL_GL_SwapWindow (Window);
    }

    return 0;
}
