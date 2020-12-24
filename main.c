#include <strings.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SCE/core/SCECore.h>


struct mesh {
    int *vertices;
    SCEvertices *verticesf;
    SCEindices *indices;
    int *v_flags;
    int n_vertices;
    int n_indices;
    int v_edgestart;
    int i_edgestart;
};

static void init_mesh (struct mesh *mesh) {
    mesh->vertices = NULL;
    mesh->indices = NULL;
    mesh->v_flags = NULL;
    mesh->n_vertices = 0;
    mesh->n_indices = 0;
    mesh->v_edgestart = 0;
    mesh->i_edgestart = 0;
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

#define SCREEN_W 1600
#define SCREEN_H 900

/* #define RAD (0.0174532925) */
static void setup_view (int rx, int ry, int dist) {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef (0.0, 0.0, -dist);
    glRotatef (rx, 1.0, 0.0, 0.0);
    glRotatef (ry, 0.0, 0.0, 1.0);
    /* glTranslatef (-10.0, -10.0, 0); */

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float matrix[16];
    proj (matrix, 70.0 * RAD, (float)SCREEN_W / SCREEN_H, 0.1, 1000.0);
    transpose (matrix);
    glLoadMatrixf (matrix);
}

#define GRID_W 4
#define GRID_H 4
#define GRID_D 4
/* maximum 4 vertices per cell */
#define n_vert (GRID_W * GRID_H * GRID_D * 4)
#define n_ind ((GRID_W - 1) * (GRID_H - 1) * 2 * 3)



static void draw (struct mesh *mesh, int use_float) {
    int coul[6] = {0x00550000,
                   0x00555000,
                   0x00005500,
                   0x00005550,
                   0x00000055,
                   0x00500055};
    /* glPolygonMode(GL_FRONT, GL_LINE); */
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_QUADS);
    for (int i = 0; i < mesh->n_indices; i++) {
        int k = mesh->indices[i];
        int c = 0x00222222;
        if (mesh->v_flags) {
            for (int j = 0; j < 6; j++) {
                if (mesh->v_flags[k] & (1 << j)) {
                    c += coul[j];
                }
            }
        }
        glColor3ub (c >> 16, c >> 8 & 0xFF, c & 0xFF);
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

void generate_hermite (uint8_t *bitmap, SCE_SMDCHermiteData *hermite) {
    memset (bitmap, 0, 4*4*4);
    /* bitmap[2] = 1 << 5 | 1 << 6; */
    /* bitmap[5] = 1 << 1 | 1 << 2; */
    bitmap[2] = 1 << 5;
    bitmap[3] = 1 << 2;
    bitmap[5] = 1 << 1;
    bitmap[4] = 1 << 6;
    /* SCE_Vector4_Set (hermite[OFFSET(1, 0, 1, 4, 4)].normalw[1], 0.0, -1.0, 0.0, 0.5 */
}

int main (void) {
    SDL_Window *Window = NULL;
    SDL_GLContext glContext;
    const int ww = SCREEN_W, wh = SCREEN_H;

    SDL_Init (SDL_INIT_VIDEO);
    Window = SDL_CreateWindow ("LE MAO", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                               ww, wh, SDL_WINDOW_OPENGL);
    glContext = SDL_GL_CreateContext (Window);
    SDL_ShowWindow (Window);

    glClearColor (0.0f, 0.0f, 0.0f, 0.0f);
    glClear (GL_COLOR_BUFFER_BIT);
    glViewport (0.0, 0.0, ww, wh);
    SDL_Event ev;
    int running = 1;

    srand(1547917722);

    glEnable(GL_DEPTH_TEST);

    int i;
    int grid_vertices[n_vert * 3];
    int grid_indices[n_ind];
    int grid_v_flags[n_vert];
    struct mesh mesh;

    init_mesh (&mesh);

    /* mesh.vertices = grid_vertices; */
    /* mesh.indices = grid_indices; */
    /* mesh.v_flags = grid_v_flags; */
    /* mesh.n_vertices = n_vert; */
    /* mesh.n_indices = n_ind; */

    {
        SCE_SMDCGenerator gen;
        uint8_t *bitmap = malloc (4 * 4 * 4);
        SCE_SMDCHermiteData *hermite = malloc (4 * 4 * 4 * sizeof *hermite);

        generate_hermite (bitmap, hermite);

        SCE_MDC_Init (&gen);
        SCE_MDC_Build (&gen, 4, 4, 4);

        mesh.n_vertices = SCE_MDC_ComputeNumVertices (&gen, bitmap);
        printf ("\n--\n");
        mesh.verticesf = malloc (mesh.n_vertices * sizeof(float) * 3);
        SCE_MDC_GenerateVertices (&gen, bitmap, hermite, mesh.verticesf);
        mesh.n_indices = SCE_MDC_ComputeNumIndices (&gen, bitmap);
        mesh.indices = malloc (mesh.n_indices * sizeof *mesh.indices);
        mesh.n_indices = SCE_MDC_GenerateIndices (&gen, bitmap, mesh.indices);
        printf ("\n");
        printf ("n_indices = %d\n", mesh.n_indices);
        printf ("n_vertices = %d\n", mesh.n_vertices);
        /* for (int i = 0; i < mesh.n_indices; i++) { */
        /*     printf ("%d ", mesh.indices[i]); */
        /* } */
        /* printf ("\n"); */
    }

    int prev_x = 0, prev_y = 0, ry = 0, rx = 0, mouse_pressed = 0;
    float dist = 10.0;
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

        setup_view (rx, ry, dist);

        draw (&mesh, 1);

        SDL_GL_SwapWindow (Window);
    }

    return 0;

}
