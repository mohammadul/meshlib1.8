#ifndef HEADER_BC297A1C5F3E1461
#define HEADER_BC297A1C5F3E1461

/**
 * @file meshlib.h
 * @author Sk. Mohammadul Haque
 * @version 1.8.0.0
 * @copyright
 * Copyright (c) 2013-2021 Sk. Mohammadul Haque.
 * @brief This header file contains declarations of all functions of meshlib.
 */

#ifndef __MESHLIB__
#define __MESHLIB__

#define _CRT_SECURE_NO_DEPRECATE

#ifdef __cplusplus
#define __MESH__CPP__
extern "C"
{
#endif


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#define MESHLIBAPI  extern

#if defined (GCC) || defined (__GNUC__)
typedef FILE *FILEPOINTER; /**< File pointer */
#else
typedef struct _iobuf *FILEPOINTER; /**< File pointer */
#endif


#define MESH_INTDATA_TYPE 1 /**< Integer datatype selector */
#define MESH_FLOATDATA_TYPE 1 /**< Float datatype selector */

#if MESH_INTDATA_TYPE==0
#define INTDATA int32_t /* do not change this, careful see meshload fscanf and other functions */ /**< Integer datatype */
#define MESH_INTDATA_MIN (INT32_MIN) /**< Integer datatype minimum */
#define MESH_INTDATA_MAX (INT32_MAX) /**< Integer datatype maximum */
#else
#define INTDATA int64_t /* do not change this, careful see meshload fscanf and other functions */ /**< Integer datatype */
#define MESH_INTDATA_MIN (INT64_MIN) /**< Integer datatype minimum */
#define MESH_INTDATA_MAX (INT64_MAX) /**< Integer datatype maximum */
#endif

#if MESH_FLOATDATA_TYPE==0
#define FLOATDATA float /* do not change this, careful see meshload fscanf and other functions */ /**< Float datatype */
#define MESH_FLOATDATA_MIN (FLT_MIN) /**< Float datatype minimum */
#define MESH_FLOATDATA_MAX (FLT_MAX) /**< Float datatype maximum */
#else
#define FLOATDATA double /* do not change this, careful see meshload fscanf and other functions */ /**< Float datatype */
#define MESH_FLOATDATA_MIN (DBL_MIN) /**< Float datatype minimum */
#define MESH_FLOATDATA_MAX (DBL_MAX) /**< Float datatype maximum */
#endif

#define MESH_ORIGIN_TYPE_BUILD 00 /**< Mesh origin type - create new */

#define MESH_ORIGIN_TYPE_OFF 11 /**< Mesh origin type - OFF file */
#define MESH_ORIGIN_TYPE_NOFF 12 /**< Mesh origin type - NOFF file */
#define MESH_ORIGIN_TYPE_COFF 13 /**< Mesh origin type - COFF file */
#define MESH_ORIGIN_TYPE_NCOFF 14 /**< Mesh origin type - NCOFF file */

#define MESH_ORIGIN_TYPE_XYZ 20 /**< Mesh origin type - XYZ file */

#define MESH_ORIGIN_TYPE_PLY_ASCII 30 /**< Mesh origin type - PLY ascii file */
#define MESH_ORIGIN_TYPE_PLY_BINARY_LITTLE_ENDIAN 31 /**< Mesh origin type - PLY binary LE file */
#define MESH_ORIGIN_TYPE_PLY_BINARY_BIG_ENDIAN 32 /**< Mesh origin type - PLY binary BE file */

#define MESH_ORIGIN_TYPE_BINV1 40 /**< Mesh origin type - CloudCompare BINv1file */
#define MESH_ORIGIN_TYPE_BUNDLE_OUT 50 /**< Mesh origin type - Bundle OUT file */
#define MESH_ORIGIN_TYPE_NVM 60 /**< Mesh origin type - NVM file */
#define MESH_ORIGIN_TYPE_BINCOLMAP 70 /**< Mesh origin type - COLMAP BIN file */

#define MESH_ERR_MALLOC 0 /**< Mesh error type - allocation */
#define MESH_ERR_SIZE_MISMATCH 1 /**< Mesh error type - size mismatch */
#define MESH_ERR_FNOTOPEN 2 /**< Mesh error type - file open */
#define MESH_ERR_INCOMPATIBLE 3 /**< Mesh error type - incompatible data */
#define MESH_ERR_UNKNOWN 4 /**< Mesh error type - unknown */

#define MESH_PI (3.14159265359) /**< \f$ \pi \f$ */
#define MESH_TWOPI (6.28318530718) /**< \f$ 2\pi \f$ */

/* 0xABCD, A3&D are vertices, B&C are faces,  */
#define MESH_PROPS_VERTICES (0x01) /**< Mesh vertices */
/** \cond HIDDEN_SYMBOLS */
#define __MESH_PROPS_VNORMALS (0x02)
#define __MESH_PROPS_VCOLORS (0x04)
#define __MESH_PROPS_VFACES (0x08)
#define __MESH_PROPS_VSCALARS (0x4000)
/** \endcond */
#define MESH_PROPS_VNORMALS (MESH_PROPS_VERTICES | __MESH_PROPS_VNORMALS) /**< Mesh vertices and vertex normals */
#define MESH_PROPS_VCOLORS (MESH_PROPS_VERTICES | __MESH_PROPS_VCOLORS) /**< Mesh vertices and vertex colors */
#define MESH_PROPS_VFACES (MESH_PROPS_VERTICES | __MESH_PROPS_VFACES) /**< Mesh vertices and vertex face adjacency */
#define MESH_PROPS_VSCALARS (MESH_PROPS_VERTICES | __MESH_PROPS_VSCALARS) /**< Mesh vertices and vertex scalars */
#define MESH_PROPS_V_ALL_PROPS (0xC00F) /**< Mesh all vertex properties */
/** \cond HIDDEN_SYMBOLS */
#define __MESH_PROPS_FACES (0x10)
#define __MESH_PROPS_FNORMALS (0x20)
#define __MESH_PROPS_FCOLORS (0x40)
#define __MESH_PROPS_FAREAS (0x80)
#define __MESH_PROPS_FFACES (0x100)
#define __MESH_PROPS_FSCALARS (0x200)
#define __MESH_PROPS_F_ALL_PROPS (0xFF0)

#define __MESH_PROPS_EDGES (0x1000)

/** \endcond */
#define MESH_PROPS_FACES (MESH_PROPS_VERTICES | __MESH_PROPS_FACES) /**< Mesh faces */
#define MESH_PROPS_FNORMALS (MESH_PROPS_FACES | __MESH_PROPS_FNORMALS) /**< Mesh faces and face normals */
#define MESH_PROPS_FCOLORS (MESH_PROPS_FACES | __MESH_PROPS_FCOLORS) /**< Mesh faces and face colors */
#define MESH_PROPS_FAREAS (MESH_PROPS_FACES | __MESH_PROPS_FAREAS) /**< Mesh faces and face areas */
#define MESH_PROPS_FFACES (MESH_PROPS_FACES | __MESH_PROPS_FFACES) /**< Mesh faces and face face adjacency */
#define MESH_PROPS_FSCALARS (MESH_PROPS_FACES | __MESH_PROPS_FSCALARS) /**< Mesh vertices and face scalars */
#define MESH_PROPS_F_ALL_PROPS (MESH_PROPS_FACES | __MESH_PROPS_F_ALL_PROPS) /**< Mesh all face properties */

#define MESH_PROPS_EDGES (MESH_PROPS_VERTICES | __MESH_PROPS_FACES | __MESH_PROPS_EDGES) /**< Mesh edges */
#define MESH_PROPS_ALL_PROPS (0xFFFF) /**< Mesh all properties */

#define MESH_EPS20 (1e-20) /**< Large epsilon */
#define MESH_EPS12 (1e-12) /**< Medium epsilon */
#define MESH_EPS8 (1e-8) /**< Small epsilon */
#if MESH_FLOATDATA_TYPE==0
#define MESH_EPSM (FLT_EPSILON) /**< Machine epsilon */
#else
#define MESH_EPSM (DBL_EPSILON) /**< Machine epsilon */
#endif

#define MESH_MIN(a,b) (((a)<(b))?(a):(b)) /**< Minimum of two variables */
#define MESH_MAX(a,b) (((a)>(b))?(a):(b)) /**< Maximum of two variables */

/* deprecated macros below */
/** \cond HIDDEN_SYMBOLS */
/* 0xABCD, A3&D are vertices, B&C are faces,  */
#define MESH_CLONE_VERTICES MESH_PROPS_VERTICES

#define __MESH_CLONE_VNORMALS __MESH_PROPS_VNORMALS
#define __MESH_CLONE_VCOLORS __MESH_PROPS_VCOLORS
#define __MESH_CLONE_VFACES __MESH_PROPS_VFACES
#define __MESH_CLONE_VSCALARS __MESH_PROPS_VSCALARS

#define MESH_CLONE_VNORMALS MESH_PROPS_VNORMALS
#define MESH_CLONE_VCOLORS MESH_PROPS_VCOLORS
#define MESH_CLONE_VFACES MESH_PROPS_VFACES
#define MESH_CLONE_VSCALARS MESH_PROPS_VSCALARS
#define MESH_CLONE_V_ALL_PROPS MESH_PROPS_V_ALL_PROPS

#define __MESH_CLONE_FACES __MESH_PROPS_FACES
#define __MESH_CLONE_FNORMALS __MESH_PROPS_FNORMALS
#define __MESH_CLONE_FCOLORS __MESH_PROPS_FCOLORS
#define __MESH_CLONE_FAREAS __MESH_PROPS_FAREAS
#define __MESH_CLONE_FFACES __MESH_PROPS_FFACES
#define __MESH_CLONE_FSCALARS __MESH_PROPS_FSCALARS
#define __MESH_CLONE_F_ALL_PROPS __MESH_PROPS_F_ALL_PROPS

#define __MESH_CLONE_EDGES __MESH_PROPS_EDGES

#define MESH_CLONE_FACES MESH_PROPS_FACES
#define MESH_CLONE_FNORMALS MESH_PROPS_FNORMALS
#define MESH_CLONE_FCOLORS MESH_PROPS_FCOLORS
#define MESH_CLONE_FAREAS MESH_PROPS_FAREAS
#define MESH_CLONE_FFACES MESH_PROPS_FFACES
#define MESH_CLONE_FSCALARS MESH_PROPS_FSCALARS
#define MESH_CLONE_F_ALL_PROPS MESH_PROPS_F_ALL_PROPS

#define MESH_CLONE_EDGES MESH_PROPS_EDGES
#define MESH_CLONE_ALL_PROPS MESH_PROPS_ALL_PROPS

/** \endcond */

#define MESH_ALIGN_GLOBAL_POSITION (0x01) /**< Mesh alignment global position */
#define MESH_ALIGN_GLOBAL_ORIENTATION (0x02) /**< Mesh alignment global orientation */
#define MESH_ALIGN_GLOBAL_SCALE (0x04) /**< Mesh alignment global scale */
#define MESH_ALIGN_GLOBAL_ALL (0x07) /**< Mesh alignment all global properties */
#define MESH_ALIGN_GLOBAL_DO_TRANSFORM (0x08) /**< Mesh alignment global do transform */

typedef INTDATA INTDATA2[2]; /**< 2- element INTDATA */
typedef INTDATA INTDATA3[3]; /**< 3- element INTDATA */

typedef struct mesh_vector2
{
    FLOATDATA x; /**< x co-ordinate */
    FLOATDATA y; /**< y co-ordinate */
} mesh_vector2; /**< Generic 2-d vector */
typedef mesh_vector2* MESH_VECTOR2; /**< Generic 2-d vector pointer */

typedef struct mesh_vector3
{
    FLOATDATA x; /**< x co-ordinate */
    FLOATDATA y; /**< y co-ordinate */
    FLOATDATA z; /**< z co-ordinate */
} mesh_vector3; /**< Generic 3-d vector */
typedef mesh_vector3* MESH_VECTOR3; /**< Generic 3-d vector pointer */

typedef mesh_vector3 mesh_vertex; /**< Vertex */
typedef mesh_vertex* MESH_VERTEX; /**< Vertex pointer */

typedef mesh_vector3 mesh_normal; /**< Normal */
typedef mesh_normal* MESH_NORMAL; /**< Normal pointer */

typedef struct mesh_color
{
    FLOATDATA r; /**< Red channel */
    FLOATDATA g; /**< Blue channel */
    FLOATDATA b; /**< Green channel */
    FLOATDATA a; /**< Alpha channel */
} mesh_color;
typedef mesh_color* MESH_COLOR; /**< Color */

typedef struct mesh_struct
{
    INTDATA num_items; /**< Number of items */
    INTDATA *items; /**< Pointer to INTDATA items */
} mesh_struct; /**< INTDATA Structure */
typedef mesh_struct* MESH_STRUCT; /**< INTDATA Structure pointer */

typedef struct mesh_struct2
{
    INTDATA num_items; /**< Number of items */
    INTDATA2 *items; /**< Pointer to INTDATA2 items */
} mesh_struct2; /**< INTDATA2 Structure */
typedef mesh_struct2* MESH_STRUCT2; /**< INTDATA2 Structure pointer */

typedef struct mesh_struct3
{
    INTDATA num_items; /**< Number of items */
    INTDATA3 *items; /**< Pointer to INTDATA3 items */
} mesh_struct3; /**< INTDATA3 Structure */
typedef mesh_struct3* MESH_STRUCT3; /**< INTDATA3 Structure pointer */

typedef FLOATDATA mesh_scalar; /**< Scalar */
typedef mesh_scalar* MESH_SCALAR; /**< Pointer to Scalar */

typedef struct mesh_face
{
    INTDATA num_vertices; /**< Number of vertices */
    INTDATA* vertices; /**< Pointer to vertex indices */
} mesh_face; /**< Face */
typedef mesh_face* MESH_FACE; /**< Pointer to face */

typedef struct mesh_edge
{
    INTDATA vertices[2]; /**< Edge vertices */
    INTDATA faces[2]; /**< Edge faces */
} mesh_edge; /**< Edge */
typedef struct mesh_edge* MESH_EDGE; /**< Pointer to edge */

typedef struct mesh_adjface
{
    INTDATA num_faces; /**< Number of adjacent faces */
    INTDATA *faces; /**< Pointer to adjacent face indices */
} mesh_adjface; /**< Adjacent face structure */

typedef struct mesh_adjface mesh_vface; /**< Vertex adjacent faces */
typedef mesh_vface* MESH_VFACE; /**< Pointer to vertex adjacent faces */

typedef struct mesh_adjface mesh_fface; /**< Face adjacent faces */
typedef mesh_fface* MESH_FFACE; /**< Pointer to face adjacent faces */


typedef struct mesh_rotation
{
    FLOATDATA data[9]; /**< 3X3 Row-major matrix data */
} mesh_rotation; /**< Rotation */
typedef mesh_rotation* MESH_ROTATION; /**< Pointer to rotation */

typedef struct mesh_affine
{
    FLOATDATA data[12]; /**< 3X4 Row-major matrix data */
} mesh_affine; /**< Affine transformation */
typedef mesh_affine* MESH_AFFINE; /**< Pointer to affine transformation */

typedef struct mesh_affine mesh_rigid; /**< Rigid transformation */
typedef mesh_rigid* MESH_RIGID; /**< Pointer to rigid transformation */

typedef struct mesh
{
    uint8_t origin_type; /**< Origin type */
    uint8_t is_loaded; /**< Is loaded? */
    uint8_t is_vertices; /**< Has vertices? */
    uint8_t is_faces; /**< Has faces? */
    uint8_t is_edges; /**< Has edges? */

    uint8_t is_vnormals; /**< Has vertex normals? */
    uint8_t is_fnormals; /**< Has face normals? */
    uint8_t is_vcolors; /**< Has vertex colors? */
    uint8_t is_fcolors; /**< Has face colors? */
    uint8_t is_vfaces; /**< Has vertex adjacent faces? */
    uint8_t is_ffaces; /**< Has face adjacent faces? */
    uint8_t is_fareas; /**< Has face areas? */
    uint8_t is_vscalars; /**< Has vertex scalar field scalars? */
    uint8_t is_fscalars; /**< Has face scalar field scalars? */

    INTDATA num_vertices; /**< Number of vertices */
    INTDATA num_faces; /**< Number of faces */
    INTDATA num_edges; /**< Number of edges */

    MESH_VERTEX vertices; /**< Pointer to vertices */
    MESH_FACE faces; /**< Pointer to faces */
    MESH_EDGE edges; /**< Pointer to edges */
    MESH_NORMAL vnormals; /**< Pointer to vertex normals */
    MESH_NORMAL fnormals; /**< Pointer to face normals */
    MESH_COLOR vcolors; /**< Pointer to vertex colors */
    MESH_COLOR fcolors; /**< Pointer to face colors */

    MESH_VFACE vfaces; /**< Pointer to vertex adjacent faces */
    MESH_FFACE ffaces; /**< Pointer to face adjacent faces */
    MESH_SCALAR fareas; /**< Pointer to face areas */
    MESH_SCALAR vscalars; /**< Pointer to vertex scalar field scalars */
    MESH_SCALAR fscalars; /**< Pointer to face scalar field scalars */

    uint8_t is_trimesh; /**< Is trimesh? */
    uint8_t dummy;
} mesh; /**< Mesh */
typedef mesh* MESH; /**< Pointer to mesh */

MESHLIBAPI void mesh_error(int type);

MESHLIBAPI MESH mesh_create_mesh_new();
MESHLIBAPI void mesh_free_mesh(MESH m);
MESHLIBAPI MESH mesh_create_mesh_new_grid(MESH_VECTOR3 sz, MESH_VECTOR3 pos, INTDATA m, INTDATA n);
MESHLIBAPI MESH mesh_create_mesh_new_cuboid(MESH_VECTOR3 sz, MESH_VECTOR3 pos);
MESHLIBAPI MESH mesh_create_mesh_new_ellipsoid(MESH_VECTOR3 sz, MESH_VECTOR3 pos);
MESHLIBAPI MESH mesh_create_mesh_new_uniform_ellipsoid(MESH_VECTOR3 sz, MESH_VECTOR3 pos, int n);
MESHLIBAPI MESH mesh_create_mesh_new_cylinder(MESH_VECTOR3 sz, MESH_VECTOR3 pos);
MESHLIBAPI MESH mesh_create_mesh_new_cone(MESH_VECTOR3 sz, MESH_VECTOR3 pos);
MESHLIBAPI MESH mesh_create_mesh_new_rectangle_flat(MESH_VECTOR3 sz, MESH_VECTOR3 pos);
MESHLIBAPI MESH mesh_create_mesh_new_ellipse_flat(MESH_VECTOR3 sz, MESH_VECTOR3 pos);

MESHLIBAPI MESH mesh_clone_mesh(MESH m, uint16_t flags);
MESHLIBAPI MESH mesh_combine_mesh(MESH m1, MESH m2);
MESHLIBAPI int mesh_alloc_mesh_props(MESH m, INTDATA nv, INTDATA nf, INTDATA ne, uint16_t flags);
MESHLIBAPI int mesh_free_mesh_props(MESH m, uint16_t flags);
MESHLIBAPI int mesh_alloc_face_vertices(MESH_FACE mf, INTDATA nv);
MESHLIBAPI int mesh_free_face_vertices(MESH_FACE mf);

MESHLIBAPI MESH mesh_load_file(const char* fname);

MESHLIBAPI MESH mesh_load_off(const char* fname);
/** \cond HIDDEN_SYMBOLS */
MESHLIBAPI MESH __mesh_parse_off_header(MESH m, FILEPOINTER fp);
MESHLIBAPI MESH __mesh_parse_off_vertices(MESH m, FILEPOINTER fp);
MESHLIBAPI MESH __mesh_parse_off_faces(MESH m, FILEPOINTER fp);
/** \endcond */

MESHLIBAPI MESH mesh_load_xyz(const char* fname);
/** \cond HIDDEN_SYMBOLS */
MESHLIBAPI MESH __mesh_parse_xyz_data(MESH m, FILEPOINTER fp);
/** \endcond */

MESHLIBAPI MESH mesh_load_ply(const char* fname);
/** \cond HIDDEN_SYMBOLS */
MESHLIBAPI MESH __mesh_parse_ply_header(MESH m, FILEPOINTER fp, int* dtypes);
MESHLIBAPI MESH __mesh_parse_ply_body(MESH m, FILEPOINTER fp, int* dtypes);
MESHLIBAPI MESH __mesh_parse_ply_vertices(MESH m, FILEPOINTER fp, int* dtypes);
MESHLIBAPI MESH __mesh_parse_ply_faces(MESH m, FILEPOINTER fp, int* dtypes);
/** \endcond */

MESHLIBAPI MESH mesh_load_bin(const char* fname);
MESHLIBAPI MESH mesh_load_out(const char* fname);
MESHLIBAPI MESH mesh_load_nvm(const char* fname);
MESHLIBAPI MESH mesh_load_colmap(const char* fname);

MESHLIBAPI int mesh_save_file(MESH m, const char* fname);

MESHLIBAPI int mesh_save_off(MESH m, const char* fname);
MESHLIBAPI int mesh_save_xyz(MESH m, const char* fname);
MESHLIBAPI int mesh_save_ply(MESH m, const char* fname);
MESHLIBAPI int mesh_save_bin(MESH m, const char* fname);
MESHLIBAPI int mesh_save_obj(MESH m, const char* fname);

MESHLIBAPI int mesh_calc_vertex_normals(MESH m);
MESHLIBAPI int mesh_calc_face_normals(MESH m);
MESHLIBAPI int mesh_calc_edges(MESH m);
MESHLIBAPI int mesh_calc_vertex_adjacency(MESH m);
MESHLIBAPI int mesh_calc_face_adjacency(MESH m);
MESHLIBAPI int mesh_upsample(MESH m, int iters);
MESHLIBAPI int mesh_upsample_loop(MESH m, int iters);
MESHLIBAPI int mesh_upsample_tarea_adaptive(MESH m, int miters, FLOATDATA e);
MESHLIBAPI void mesh_cross_vector3(MESH_VECTOR3 x, MESH_VECTOR3 y, MESH_VECTOR3 z);
MESHLIBAPI void mesh_cross_normal(MESH_NORMAL x, MESH_NORMAL y, MESH_NORMAL z);
MESHLIBAPI FLOATDATA mesh_calc_triangle_area(MESH_VERTEX a, MESH_VERTEX b, MESH_VERTEX c);
MESHLIBAPI void mesh_calc_face_normal(MESH_VERTEX v1, MESH_VERTEX v2, MESH_VERTEX v3, MESH_NORMAL n);
MESHLIBAPI void mesh_calc_aabb(MESH m, MESH_VECTOR3 minv, MESH_VECTOR3 maxv, MESH_VECTOR3 center);
MESHLIBAPI int mesh_calc_signed_area(MESH m, MESH_VECTOR3 area);
MESHLIBAPI FLOATDATA mesh_calc_area(MESH m);
MESHLIBAPI FLOATDATA mesh_calc_volume(MESH m);

MESHLIBAPI INTDATA mesh_find(MESH_STRUCT s, INTDATA q);
MESHLIBAPI INTDATA mesh_find2(MESH_STRUCT2 s, INTDATA q);
MESHLIBAPI INTDATA mesh_find3(MESH_STRUCT3 s, INTDATA q);

MESHLIBAPI int mesh_remove_boundary_vertices(MESH m, int iters);
MESHLIBAPI int mesh_remove_boundary_faces(MESH m, int iters);
MESHLIBAPI int mesh_remove_triangles_with_small_area(MESH m, FLOATDATA area);
MESHLIBAPI int mesh_remove_unreferenced_vertices(MESH m);
MESHLIBAPI int mesh_remove_zero_area_faces(MESH m);
MESHLIBAPI int mesh_remove_close_vertices(MESH m, FLOATDATA r);
MESHLIBAPI int mesh_remove_ear_faces(MESH m, int niters);
MESHLIBAPI int mesh_remove_non_manifold_vertices(MESH m);

MESHLIBAPI int mesh_isnumeric(FILEPOINTER fp);
MESHLIBAPI int mesh_go_next_word(FILEPOINTER fp);
MESHLIBAPI int mesh_read_word(FILEPOINTER fp, char *c_word, int sz);
MESHLIBAPI int mesh_read_word_skip_comment(FILEPOINTER fp, char *c_word, int sz);
MESHLIBAPI int mesh_read_word_only(FILEPOINTER fp, char *c_word, int sz);
MESHLIBAPI int mesh_read_word_only_skip_comment(FILEPOINTER fp, char *c_word, int sz);
MESHLIBAPI int mesh_count_words_in_line(FILEPOINTER fp, int *count);
MESHLIBAPI int mesh_skip_line(FILEPOINTER fp);

MESHLIBAPI int mesh_filter_bilateral(MESH m, FLOATDATA sigma_c, FLOATDATA sigma_s, int niters);
MESHLIBAPI int mesh_filter_laplacian(MESH m, FLOATDATA r);
MESHLIBAPI int mesh_filter_laplacian_restricted(MESH m, FLOATDATA r, FLOATDATA ang);
MESHLIBAPI int mesh_filter_laplacian_depth(MESH m, FLOATDATA r, MESH_VECTOR3 vp);
MESHLIBAPI int mesh_filter_taubin(MESH m, FLOATDATA lambd, FLOATDATA mu, int niters);
MESHLIBAPI int mesh_filter_vertex_color_bilateral(MESH m, FLOATDATA sigma_k, FLOATDATA sigma_c, FLOATDATA sigma_s, int niters);
MESHLIBAPI int mesh_filter_vertex_color_laplacian(MESH m, FLOATDATA r);
MESHLIBAPI int mesh_filter_vertex_color_min(MESH m, INTDATA niters);
MESHLIBAPI int mesh_filter_vertex_color_max(MESH m, INTDATA niters);

MESHLIBAPI MESH_ROTATION mesh_rotation_create();
MESHLIBAPI void mesh_rotation_free(MESH_ROTATION r);
MESHLIBAPI MESH_ROTATION mesh_rotation_set_matrix(FLOATDATA *mat, MESH_ROTATION r);
MESHLIBAPI MESH_ROTATION mesh_rotation_set_angleaxis(FLOATDATA ang, MESH_NORMAL axis, MESH_ROTATION r);

MESHLIBAPI int mesh_translate(MESH m, FLOATDATA x, FLOATDATA y, FLOATDATA z);
MESHLIBAPI int mesh_translate_vector(MESH m, MESH_VECTOR3 v);
MESHLIBAPI int mesh_scale(MESH m, FLOATDATA sx, FLOATDATA sy, FLOATDATA sz);

MESHLIBAPI MESH_VERTEX mesh_vertex_rotate(MESH_VERTEX v, MESH_ROTATION r);
MESHLIBAPI int mesh_rotate(MESH m, MESH_ROTATION r);

MESHLIBAPI MESH_AFFINE mesh_affine_create();
MESHLIBAPI void mesh_affine_free(MESH_AFFINE tx);
MESHLIBAPI MESH_AFFINE mesh_affine_set_matrix(FLOATDATA *mat, MESH_AFFINE r);
MESHLIBAPI int mesh_transform(MESH m, MESH_AFFINE tx);
MESHLIBAPI MESH_AFFINE mesh_affine_set_rotation_translation(MESH_ROTATION r, MESH_VECTOR3 t, MESH_AFFINE tx);

MESHLIBAPI void mesh_set_seed(int seed);
MESHLIBAPI FLOATDATA __mesh_rand(void);
MESHLIBAPI FLOATDATA __mesh_randn(void);
MESHLIBAPI FLOATDATA __mesh_randexp(void);
MESHLIBAPI FLOATDATA __mesh_randfun(FLOATDATA(*fun)(FLOATDATA), FLOATDATA xmin, FLOATDATA xmax);

MESHLIBAPI int mesh_add_noise_uniform(MESH m, FLOATDATA sigma);
MESHLIBAPI int mesh_add_noise_gaussian(MESH m, FLOATDATA sigma);
MESHLIBAPI int mesh_add_noise_exp(MESH m, FLOATDATA sigma);
MESHLIBAPI int mesh_add_noise_func(MESH m, FLOATDATA sigma, FLOATDATA(*fun)(FLOATDATA), FLOATDATA xmin, FLOATDATA xmax);
MESHLIBAPI int mesh_add_noise_uniform_normal(MESH m, FLOATDATA sigma);
MESHLIBAPI int mesh_add_noise_gaussian_normal(MESH m, FLOATDATA sigma);
MESHLIBAPI int mesh_add_noise_exp_normal(MESH m, FLOATDATA sigma);
MESHLIBAPI int mesh_add_noise_func_normal(MESH m, FLOATDATA sigma, FLOATDATA(*fun)(FLOATDATA), FLOATDATA xmin, FLOATDATA xmax);
MESHLIBAPI int mesh_add_noise_uniform_tangent(MESH m, FLOATDATA sigma);
MESHLIBAPI int mesh_add_noise_gaussian_tangent(MESH m, FLOATDATA sigma);
MESHLIBAPI int mesh_add_noise_exp_tangent(MESH m, FLOATDATA sigma);
MESHLIBAPI int mesh_add_noise_func_tangent(MESH m, FLOATDATA sigma, FLOATDATA(*fun)(FLOATDATA), FLOATDATA xmin, FLOATDATA xmax);

MESHLIBAPI void mesh_draw_mesh(MESH m);
MESHLIBAPI void mesh_draw_mesh_smooth(MESH m);
MESHLIBAPI void mesh_draw_point_cloud(MESH m);

MESHLIBAPI int mesh_align_global(MESH m1, MESH m2, int flags, MESH_AFFINE tx);

/* deprecated function names */
#define mesh_bilateral_filter mesh_filter_bilateral
#define mesh_laplacian_filter mesh_filter_laplacian
#define mesh_restricted_laplacian_filter mesh_filter_laplacian_restricted
#define mesh_depth_laplacian_filter mesh_filter_laplacian_depth
#define mesh_taubin_filter mesh_filter_taubin
#define mesh_bilateral_vertex_color_filter mesh_filter_vertex_color_bilateral
#define mesh_laplacian_vertex_color_filter mesh_filter_vertex_color_laplacian
#define mesh_min_vertex_color_filter mesh_filter_vertex_color_min
#define mesh_max_vertex_color_filter mesh_filter_vertex_color_max

#define mesh_write_file mesh_save_file
#define mesh_write_off mesh_save_off
#define mesh_write_xyz mesh_save_xyz
#define mesh_write_ply mesh_save_ply
#define mesh_write_bin mesh_save_bin
#define mesh_write_obj mesh_save_obj

#define mesh_read_file mesh_load_file
#define mesh_read_off mesh_load_off
#define mesh_read_xyz mesh_load_xyz
#define mesh_read_ply mesh_load_ply
#define mesh_read_bin mesh_load_bin
#define mesh_read_out mesh_load_out
#define mesh_read_nvm mesh_load_nvm
#define mesh_read_colmap mesh_load_colmap

#ifdef __cplusplus
}
#endif

#endif
#endif // header guard

