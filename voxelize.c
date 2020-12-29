#include <strings.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <SCE/core/SCECore.h>
#include <SCE/core/libwar.h>

void writeint (uint32_t n) {
    fwrite (&n, 1, sizeof n, stdout);
}

static void write_hermite (uint8_t *bitmap, SCE_SMDCHermiteData *hermite,
                           uint32_t vox_w, uint32_t vox_h, uint32_t vox_d) {
    writeint (vox_w);
    writeint (vox_h);
    writeint (vox_d);
    int n_vox = vox_w * vox_h * vox_d;
    fwrite (bitmap, 1, n_vox / 8 + 1, stdout);
    fwrite (hermite, sizeof *hermite, n_vox, stdout);
}

#if 0
static void manmade (uint8_t *bitmap, SCE_SMDCHermiteData *hermite) {
    bitmap[2] = (1 << 5) | (1 << 6);
    /* bitmap[3] = (1 << 1); */

    /* ########################## */
    hermite[17].normalw[1][3] = 0.5;
    SCE_Vector3_Set (hermite[17].normalw[1], 0.0, -1.0, 0.0);
    hermite[21].normalw[1][3] = 0.5;
    SCE_Vector3_Set (hermite[21].normalw[1], 0.0, 1.0, 0.0);

    hermite[20].normalw[0][3] = 0.5;
    SCE_Vector3_Set (hermite[20].normalw[0], -1.0, 0.0, 0.0);
    /* hermite[21].normalw[0][3] = 0.5; */
    /* SCE_Vector3_Set (hermite[21].normalw[0], 1.0, 0.0, 0.0); */
    /* SCE_Vector3_Normalize (hermite[21].normalw[0]); */

    hermite[5].normalw[2][3] = 0.5;
    SCE_Vector3_Set (hermite[5].normalw[2], 0.0, 0.0, -1.0);
    hermite[21].normalw[2][3] = 0.5;
    SCE_Vector3_Set (hermite[21].normalw[2], 0.0, 0.0, 1.0);
    /* ####################### */

    hermite[18].normalw[1][3] = 0.5;
    SCE_Vector3_Set (hermite[18].normalw[1], 0.0, -1.0, 0.0);
    hermite[22].normalw[1][3] = 0.3;
    SCE_Vector3_Set (hermite[22].normalw[1], 0.4, 1.0, 0.0);
    
    hermite[6].normalw[2][3] = 0.5;
    SCE_Vector3_Set (hermite[6].normalw[2], 0.0, 0.0, -1.0);
    hermite[22].normalw[2][3] = 0.5;
    SCE_Vector3_Set (hermite[22].normalw[2], 0.0, 0.0, 1.0);

    hermite[22].normalw[0][3] = 0.5;
    SCE_Vector3_Set (hermite[22].normalw[0], 1.0, 0.0, 0.0);
    /* ####################### */


    /* hermite[21].normalw[1][3] = 0.5; */
    /* hermite[20].normalw[0][3] = 0.5; */
}
#endif

static void mk_line (SCE_SLine3 *l, int x, int y, int z, int axis) {
    SCE_TVector3 o, n;
    SCE_Vector3_Set (n, 0.0, 0.0, 0.0);
    n[axis] = 1.0;
    SCE_Vector3_Set (o, x, y, z);
    SCE_Line3_SetOrigin (l, o);
    SCE_Line3_SetNormal (l, n);
}

static int intersects (WarMesh *mesh, SCE_SLine3 *l, SCE_TVector3 n, float *d, int axis) {
    SCE_TVector3 p, nor;
    int ret = 0;
    /* int index; */

    *d = 1.0;
    for (int i = 0; i < mesh->vcount; i += 3) {
        if (SCE_Plane_TriangleLineIntersection (&mesh->pos[i * 3], &mesh->pos[(i + 1) * 3],
                                                &mesh->pos[(i + 2) * 3], l, p)) {
            float dist = p[axis] - l->o[axis];
            if (dist < 1.0 && dist > 0.0) {
                SCE_SPlane tri;
                SCE_Plane_SetFromTriangle (&tri, &mesh->pos[i * 3], &mesh->pos[(i + 1) * 3], &mesh->pos[(i + 2) * 3]);
                SCE_Plane_GetNormalv (&tri, nor);
                if (nor[axis] < 0.0)
                    ret++;
                else
                    ret--;
                if (dist < *d) {
                    *d = dist;
                    /* index = i; */
                    SCE_Vector3_Copy (n, nor);
                }
            }
        }
    }
#if 0
    if (ret) {
        *d = min;
        SCE_Plane_SetFromTriangle (&tri, &mesh->pos[index * 3], &mesh->pos[(index + 1) * 3],
                                   &mesh->pos[(index + 2) * 3]);
        SCE_Plane_GetNormalv (&tri, n);
        //            fprintf (stderr, "inter: %.2f %.2f %.2f orig: %.2f %.2f %.2f fnor: %.2f %.2f %.2f lnor: %.2f %.2f %.2f  # %.2f\n", p[0], p[1], p[2], l->o[0], l->o[1], l->o[2], n[0], n[1], n[2],  l->n[0], l->n[1], l->n[2], *d);
    }
#endif
    return ret;
}

#define OFFSET(x, y, z, w, h) ((w) * ((z) * (h) + (y)) + (x))

static void voxelize_x (WarMesh *mesh, uint8_t *bitmap, SCE_SMDCHermiteData *hermite,
                        uint32_t vox_w, uint32_t vox_h, uint32_t vox_d) {
    int density;
    SCE_SLine3 l;
    float w;
    SCE_TVector3 n;
    int axis = 0;

    for (int z = 0; z < vox_d; z++) {
        for (int y = 0; y < vox_h; y++) {
            density = 0;
            for (int x = 0; x < vox_w; x++) {
                size_t offset = OFFSET(x, y, z, vox_w, vox_h);

                if (density > 0 && bitmap)
                    bitmap[offset / 8] |= 1 << (offset % 8);

                mk_line (&l, x, y, z, axis);
                int inter = intersects (mesh, &l, n, &w, axis);
                if (inter != 0) {
                    if ((density == 0 && inter > 0) || density == -inter) {
                        hermite[offset].normalw[axis][3] = w;
                        SCE_Vector3_Normalize (n);
                        SCE_Vector3_Copy (hermite[offset].normalw[axis], n);
                    }
                    density += inter;
                }
            }
        }
    }
}
static void voxelize_y (WarMesh *mesh, uint8_t *bitmap, SCE_SMDCHermiteData *hermite,
                        uint32_t vox_w, uint32_t vox_h, uint32_t vox_d) {
    int density;
    SCE_SLine3 l;
    float w;
    SCE_TVector3 n;
    int axis = 1;

    for (int z = 0; z < vox_d; z++) {
        for (int x = 0; x < vox_w; x++) {
            density = 0;
            for (int y = 0; y < vox_h; y++) {
                size_t offset = OFFSET(x, y, z, vox_w, vox_h);

                if (density > 0 && bitmap)
                    bitmap[offset / 8] |= 1 << (offset % 8);

                mk_line (&l, x, y, z, axis);
                int inter = intersects (mesh, &l, n, &w, axis);
                if (inter != 0) {
                    if ((density == 0 && inter > 0) || density == -inter) {
                        hermite[offset].normalw[axis][3] = w;
                        SCE_Vector3_Normalize (n);
                        SCE_Vector3_Copy (hermite[offset].normalw[axis], n);
                    }
                    density += inter;
                }
            }
        }
    }
}
static void voxelize_z (WarMesh *mesh, uint8_t *bitmap, SCE_SMDCHermiteData *hermite,
                        uint32_t vox_w, uint32_t vox_h, uint32_t vox_d) {
    int density;
    SCE_SLine3 l;
    float w;
    SCE_TVector3 n;
    int axis = 2;

    for (int y = 0; y < vox_h; y++) {
        for (int x = 0; x < vox_w; x++) {
            density = 0;
            for (int z = 0; z < vox_d; z++) {
                size_t offset = OFFSET(x, y, z, vox_w, vox_h);

                if (density > 0 && bitmap)
                    bitmap[offset / 8] |= 1 << (offset % 8);

                mk_line (&l, x, y, z, axis);
                int inter = intersects (mesh, &l, n, &w, axis);
                if (inter != 0) {
                    if ((density == 0 && inter > 0) || density == -inter) {
                        hermite[offset].normalw[axis][3] = w;
                        SCE_Vector3_Normalize (n);
                        SCE_Vector3_Copy (hermite[offset].normalw[axis], n);
                    }
                    density += inter;
                }
            }
        }
    }
}

#define DIM 30

int main (void) {
    uint32_t vox_w = DIM, vox_h = DIM, vox_d = DIM;
    uint32_t n_vox = vox_w * vox_h * vox_d;

    uint8_t *bitmap = malloc (n_vox / 8 + 1);
    SCE_SMDCHermiteData *hermite = malloc (n_vox * sizeof *hermite);

    memset (bitmap, 0, n_vox / 8 + 1);

    for (int i = 0; i < n_vox; i++) {
        for (int j = 0; j < 3; j++) {
            hermite[i].normalw[j][3] = 0.0;
            SCE_Vector3_Set (hermite[i].normalw[j], 0.0, 0.0, 0.0);
        }
    }

    WarMesh *mesh = war_read (stdin, 0, 0);

    voxelize_x (mesh, bitmap, hermite, vox_w, vox_h, vox_d);
    voxelize_y (mesh, NULL, hermite, vox_w, vox_h, vox_d);
    voxelize_z (mesh, NULL, hermite, vox_w, vox_h, vox_d);
    write_hermite (bitmap, hermite, vox_w, vox_h, vox_d);

    return 0;
}
