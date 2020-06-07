/**
 * @file meshlib.h
 * @author Sk. Mohammadul Haque
 * @version 1.4.2.0
 * @copyright
 * Copyright (c) 2013, 2014, 2015, 2016 Sk. Mohammadul Haque.
 * @brief This header file contains declarations of all functions of meshlib.
 */

#ifndef __MESHLIB__
#define __MESHLIB__

#define _CRT_SECURE_NO_DEPRECATE

/*! \mainpage Meshlib
 *
 * \section intro_sec Introduction
 * Meshlib is a simple mesh library written in C.
 *
 * \section build_sec Build
 * To build the whole project, Code::blocks is required.
 *
 * \section content_sec Contents
 * Load/Write PLY, OFF, ASC files.
 *
 * Basic Vertex Manipulations.
 *
 * Basic Vertex Transformations.
 *
 * Basic Face Manipulations.
 *
 * Bilateral Filtering.
 *
 * Laplacian Filtering.
 *
 * Mesh Cleaning Algorithms.
 *
 *
 */



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
#define MESHLIBAPI  extern

#if defined (GCC) || defined (__GNUC__)
typedef FILE *FILEPOINTER; /**< File pointer */
#else
typedef struct _iobuf *FILEPOINTER; /**< File pointer */
#endif


#define MESH_INTDATA_TYPE 0 /**< Integer datatype selector */
#define MESH_FLOATDATA_TYPE 1 /**< Float datatype selector */

#if MESH_INTDATA_TYPE==0
#define INTDATA int32_t /* do not change this, careful see meshload fscanf and other functions */ /**< Integer datatype */
#else
#define INTDATA int64_t /* do not change this, careful see meshload fscanf and other functions */ /**< Integer datatype */
#endif

#if MESH_FLOATDATA_TYPE==0
#define FLOATDATA float /* do not change this, careful see meshload fscanf and other functions */ /**< Float datatype */
#else
#define FLOATDATA double /* do not change this, careful see meshload fscanf and other functions */ /**< Float datatype */
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

#define MESH_ERR_MALLOC 0 /**< Mesh error type - allocation */
#define MESH_ERR_SIZE_MISMATCH 1 /**< Mesh error type - size mismatch */
#define MESH_ERR_FNOTOPEN 2 /**< Mesh error type - file open */
#define MESH_ERR_INCOMPATIBLE 3 /**< Mesh error type - incompatible data */
#define MESH_ERR_UNKNOWN 4 /**< Mesh error type - unknown */

#define MESH_PI (3.14159265359) /**< \f$ \pi \f$ */
#define MESH_TWOPI (6.28318530718) /**< \f$ 2\pi \f$ */


#define MESH_CLONE_VERTICES (0x01) /**< Clone mesh vertices */
/** \cond HIDDEN_SYMBOLS */
#define __MESH_CLONE_VNORMALS (0x02)
#define __MESH_CLONE_VCOLORS (0x04)
#define __MESH_CLONE_VFACES (0x08)
/** \endcond */
#define MESH_CLONE_VNORMALS (MESH_CLONE_VERTICES | __MESH_CLONE_VNORMALS) /**< Clone mesh vertices and vertex normals */
#define MESH_CLONE_VCOLORS (MESH_CLONE_VERTICES | __MESH_CLONE_VCOLORS) /**< Clone mesh vertices and vertex colors */
#define MESH_CLONE_VFACES (MESH_CLONE_VERTICES | __MESH_CLONE_VFACES) /**< Clone mesh vertices and vertex face adjacency */
#define MESH_CLONE_V_ALL_PROPS (0x0F) /**< Clone mesh all vertex properties */
/** \cond HIDDEN_SYMBOLS */
#define __MESH_CLONE_FACES (0x10)
#define __MESH_CLONE_FNORMALS (0x20)
#define __MESH_CLONE_FCOLORS (0x40)
#define __MESH_CLONE_FAREAS (0x80)
#define __MESH_CLONE_FFACES (0x100)
#define __MESH_CLONE_F_ALL_PROPS (0xFF0)

#define __MESH_CLONE_EDGES (0x1000)

/** \endcond */
#define MESH_CLONE_FACES (MESH_CLONE_VERTICES | __MESH_CLONE_FACES) /**< Clone mesh faces */
#define MESH_CLONE_FNORMALS (MESH_CLONE_FACES | __MESH_CLONE_FNORMALS) /**< Clone mesh faces and face normals */
#define MESH_CLONE_FCOLORS (MESH_CLONE_FACES | __MESH_CLONE_FCOLORS) /**< Clone mesh faces and face colors */
#define MESH_CLONE_FAREAS (MESH_CLONE_FACES | __MESH_CLONE_FAREAS) /**< Clone mesh faces and face areas */
#define MESH_CLONE_FFACES (MESH_CLONE_FACES | __MESH_CLONE_FFACES) /**< Clone mesh faces and face face adjacency */
#define MESH_CLONE_F_ALL_PROPS (MESH_CLONE_FACES | __MESH_CLONE_F_ALL_PROPS) /**< Clone mesh all face properties */

#define MESH_CLONE_EDGES (MESH_CLONE_VERTICES | __MESH_CLONE_FACES | __MESH_CLONE_EDGES) /**< Clone mesh edges */
#define MESH_CLONE_ALL_PROPS (0xFFFF) /**< Clone mesh all properties */


typedef INTDATA INTDATA2[2]; /**< 2- element INTDATA */
typedef INTDATA INTDATA3[3]; /**< 3- element INTDATA */

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
    FLOATDATA data[9]; /**< Matrix data */
} mesh_rotation; /**< Rotation */
typedef mesh_rotation* MESH_ROTATION; /**< Pointer to rotation */

typedef struct mesh_transform
{
    FLOATDATA *data; /**< Matrix data */
} mesh_transform; /**< Transformation */
typedef mesh_transform* MESH_TRANSFORM; /**< Pointer to transformation */

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
    FLOATDATA* fareas; /**< Pointer to face areas */

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
MESHLIBAPI MESH mesh_create_mesh_new_cylinder(MESH_VECTOR3 sz, MESH_VECTOR3 pos);
MESHLIBAPI MESH mesh_create_mesh_new_cone(MESH_VECTOR3 sz, MESH_VECTOR3 pos);

MESHLIBAPI MESH mesh_clone_mesh(MESH m, uint16_t flags);
MESHLIBAPI MESH mesh_combine_mesh(MESH m1, MESH m2);

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

MESHLIBAPI int mesh_write_file(MESH m, const char* fname);

MESHLIBAPI int mesh_write_off(MESH m, const char* fname);
MESHLIBAPI int mesh_write_xyz(MESH m, const char* fname);
MESHLIBAPI int mesh_write_ply(MESH m, const char* fname);

MESHLIBAPI int mesh_calc_vertex_normals(MESH m);
MESHLIBAPI int mesh_calc_face_normals(MESH m);
MESHLIBAPI int mesh_calc_edges(MESH m);
MESHLIBAPI int mesh_calc_vertex_adjacency(MESH m);
MESHLIBAPI int mesh_calc_face_adjacency(MESH m);
MESHLIBAPI int mesh_upsample(MESH m, int iters);
MESHLIBAPI void mesh_cross_vector3(MESH_VECTOR3 x, MESH_VECTOR3 y, MESH_VECTOR3 z);
MESHLIBAPI void mesh_cross_normal(MESH_NORMAL x, MESH_NORMAL y, MESH_NORMAL z);
MESHLIBAPI FLOATDATA mesh_calc_triangle_area(MESH_VERTEX a, MESH_VERTEX b, MESH_VERTEX c);
MESHLIBAPI void mesh_calc_face_normal(MESH_VERTEX v1, MESH_VERTEX v2, MESH_VERTEX v3, MESH_NORMAL n);

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
MESHLIBAPI int mesh_read_word_only(FILEPOINTER fp, char *c_word, int sz);
MESHLIBAPI int mesh_count_words_in_line(FILEPOINTER fp, int *count);
MESHLIBAPI int mesh_skip_line(FILEPOINTER fp);

MESHLIBAPI int mesh_bilateral_filter(MESH m, FLOATDATA sigma_c, FLOATDATA sigma_s, int niters);
MESHLIBAPI int mesh_laplacian_filter(MESH m, FLOATDATA r);
MESHLIBAPI int mesh_restricted_laplacian_filter(MESH m, FLOATDATA r, FLOATDATA ang);

MESHLIBAPI MESH_ROTATION mesh_rotation_create();
MESHLIBAPI void mesh_rotation_free(MESH_ROTATION r);
MESHLIBAPI MESH_ROTATION mesh_rotation_set_matrix(FLOATDATA *mat, MESH_ROTATION r);
MESHLIBAPI MESH_ROTATION mesh_rotation_set_angleaxis(FLOATDATA ang, MESH_NORMAL axis, MESH_ROTATION r);

MESHLIBAPI int mesh_translate(MESH m, FLOATDATA x, FLOATDATA y, FLOATDATA z);
MESHLIBAPI int mesh_translate_vector(MESH m, MESH_VERTEX v);
MESHLIBAPI int mesh_scale(MESH m, FLOATDATA sx, FLOATDATA sy, FLOATDATA sz);

MESHLIBAPI MESH_VERTEX mesh_vertex_rotate(MESH_VERTEX v, MESH_ROTATION r);
MESHLIBAPI int mesh_rotate(MESH m, MESH_ROTATION r);

MESHLIBAPI void mesh_draw_mesh(MESH m);
MESHLIBAPI void mesh_draw_mesh_smooth(MESH m);
MESHLIBAPI void mesh_draw_point_cloud(MESH m);

#ifdef __cplusplus
}
#endif

#endif
