/**
 * @file meshtransform.c
 * @author Sk. Mohammadul Haque
 * @version 1.8.0.0
 * @copyright
 * Copyright (c) 2013-2021 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to different mesh transformations.
 */

#include "../include/meshlib.h"
#include <math.h>


/** \brief Creates a new rotation
 *
 * \return Output rotation
 *
 */

MESH_ROTATION mesh_rotation_create()
{
    MESH_ROTATION r;
    if((r = (MESH_ROTATION)malloc(sizeof(mesh_rotation)))==NULL) mesh_error(MESH_ERR_MALLOC);
    return r;
}

/** \brief Frees a given rotation
 *
 * \param r Input rotation
 * \return NULL
 *
 */

void mesh_rotation_free(MESH_ROTATION r)
{
    free(r);
}


/** \brief Creates a new affine transformation
 *
 * \return Output affine transformation
 *
 */

MESH_AFFINE mesh_affine_create()
{
    MESH_AFFINE tx;
    if((tx = (MESH_AFFINE)malloc(sizeof(mesh_affine)))==NULL) mesh_error(MESH_ERR_MALLOC);
    return tx;
}

/** \brief Frees a given affine transformation
 *
 * \param tx Input affine transformation
 * \return NULL
 *
 */

void mesh_affine_free(MESH_AFFINE tx)
{
    free(tx);
}

/** \brief Sets rotation from a matrix
 *
 * \param[in] mat Input matrix
 * \param[out] r Input rotation
 * \return Output rotation
 *
 */

MESH_ROTATION mesh_rotation_set_matrix(FLOATDATA *mat, MESH_ROTATION r)
{
    int k;
    if(r==NULL) r = mesh_rotation_create();
    for(k=0; k<9; ++k) r->data[k] = mat[k];
    return r;
}

/** \brief Sets affine transformation from a matrix
 *
 * \param[in] mat Input matrix
 * \param[out] r Input affine transformation
 * \return Output affine transformation
 *
 */

MESH_AFFINE mesh_affine_set_matrix(FLOATDATA *mat, MESH_AFFINE r)
{
    int k;
    if(r==NULL) r = mesh_affine_create();
    for(k=0; k<12; ++k) r->data[k] = mat[k];
    return r;
}

/** \brief Sets rotation from angle axis
 *
 * \param[in] ang Input angle of rotation
 * \param[out] axis Input axis of rotation
 * \param[out] r Input rotation
 * \return Output rotation
 *
 */

MESH_ROTATION mesh_rotation_set_angleaxis(FLOATDATA ang, MESH_NORMAL axis, MESH_ROTATION r)
{
    FLOATDATA c, s, tmp0, tmp1, tmp2;
    if(r==NULL) r = mesh_rotation_create();

    tmp0 = sqrt(axis->x*axis->x+axis->y*axis->y+axis->z*axis->z);
    if(tmp0>0)
    {
        axis->x /= tmp0;
        axis->y /= tmp0;
        axis->z /= tmp0;
    }
    c = cos(ang);
    s = sin(ang);
    tmp0 = 1-c;
    r->data[0] = c+axis->x*axis->x*tmp0;
    r->data[4] = c+axis->y*axis->y*tmp0;
    r->data[8] = c+axis->z*axis->z*tmp0;

    tmp1 = axis->x*axis->y*tmp0;
    tmp2 = axis->z*s;
    r->data[1] = tmp1+tmp2;
    r->data[3] = tmp1-tmp2;
    tmp1 = axis->x*axis->z*tmp0;
    tmp2 = axis->y*s;
    r->data[2] = tmp1-tmp2;
    r->data[6] = tmp1+tmp2;
    tmp1 = axis->y*axis->z*tmp0;
    tmp2 = axis->x*s;
    r->data[5] = tmp1+tmp2;
    r->data[7] = tmp1-tmp2;
    return r;
}

/** \brief Translates a mesh by x, y and z amounts
 *
 * \param[in] m Input mesh
 * \param[in] x X component
 * \param[in] y Y component
 * \param[in] z Z component
 * \return Error code
 *
 */

int mesh_translate(MESH m, FLOATDATA x, FLOATDATA y, FLOATDATA z)
{
    INTDATA i;
    if(m==NULL) return 1;
    for(i=0; i<m->num_vertices; ++i)
    {
        m->vertices[i].x += x;
        m->vertices[i].y += y;
        m->vertices[i].z += z;
    }
    return 0;
}

/** \brief Translates a mesh by a given 3-d vector
 *
 * \param[in] m Input mesh
 * \param[in] v Input vector
 * \return Error code
 *
 */

int mesh_translate_vector(MESH m, MESH_VECTOR3 v)
{
    INTDATA i;
    if(m==NULL) return 1;
    for(i=0; i<m->num_vertices; ++i)
    {
        m->vertices[i].x += v->x;
        m->vertices[i].y += v->y;
        m->vertices[i].z += v->z;
    }
    return 0;
}

/** \brief Scales a mesh by x, y and z amounts
 *
 * \param[in] m Input mesh
 * \param[in] sx X component
 * \param[in] sy Y component
 * \param[in] sz Z component
 * \return Error code
 *
 */

int mesh_scale(MESH m, FLOATDATA sx, FLOATDATA sy, FLOATDATA sz)
{
    INTDATA i;
    if(m==NULL) return 1;
    for(i=0; i<m->num_vertices; ++i)
    {
        m->vertices[i].x *= sx;
        m->vertices[i].y *= sy;
        m->vertices[i].z *= sz;
    }
    return 0;
}

/** \brief Rotates a vertex by a given rotation
 *
 * \param[in] v Input vertex
 * \param[in] r Input rotation
 * \return Output vertex
 *
 */

MESH_VERTEX mesh_vertex_rotate(MESH_VERTEX v, MESH_ROTATION r)
{
    FLOATDATA x, y, z;
    x = r->data[0]*v->x+r->data[3]*v->y+r->data[6]*v->z;
    y = r->data[1]*v->x+r->data[4]*v->y+r->data[7]*v->z;
    z = r->data[2]*v->x+r->data[5]*v->y+r->data[8]*v->z;

    v->x = x;
    v->y = y;
    v->z = z;
    return v;
}

/** \brief Rotates a mesh by a given rotation
 *
 * \param[in] m Input vertex
 * \param[in] r Input rotation
 * \return Error code
 *
 */

int mesh_rotate(MESH m, MESH_ROTATION r)
{
    INTDATA i;
    if(m==NULL) return 1;
    for(i=0; i<m->num_vertices; ++i)
    {
        FLOATDATA x, y, z;
        x = r->data[0]*m->vertices[i].x+r->data[3]*m->vertices[i].y+r->data[6]*m->vertices[i].z;
        y = r->data[1]*m->vertices[i].x+r->data[4]*m->vertices[i].y+r->data[7]*m->vertices[i].z;
        z = r->data[2]*m->vertices[i].x+r->data[5]*m->vertices[i].y+r->data[8]*m->vertices[i].z;

        m->vertices[i].x = x;
        m->vertices[i].y = y;
        m->vertices[i].z = z;
    }
    if(m->is_vnormals) mesh_calc_vertex_normals(m);
    if(m->is_fnormals) mesh_calc_face_normals(m);

    return 0;
}

/** \brief Transforms a mesh by a given affine transformation
 *
 * \param[in] m Input vertex
 * \param[in] tx Input affine transformation
 * \return Error code
 *
 */

int mesh_transform(MESH m, MESH_AFFINE tx)
{
    INTDATA i;
    if(m==NULL) return 1;
    for(i=0; i<m->num_vertices; ++i)
    {
        FLOATDATA x, y, z;
        x = tx->data[0]*m->vertices[i].x+tx->data[3]*m->vertices[i].y+tx->data[6]*m->vertices[i].z;
        y = tx->data[1]*m->vertices[i].x+tx->data[4]*m->vertices[i].y+tx->data[7]*m->vertices[i].z;
        z = tx->data[2]*m->vertices[i].x+tx->data[5]*m->vertices[i].y+tx->data[8]*m->vertices[i].z;

        m->vertices[i].x = x+tx->data[9];
        m->vertices[i].y = y+tx->data[10];
        m->vertices[i].z = z+tx->data[11];
    }
    if(m->is_vnormals) mesh_calc_vertex_normals(m);
    if(m->is_fnormals) mesh_calc_face_normals(m);

    return 0;
}

/** \brief Sets an affine transformation with rotation and translation
 *
 * \param[in] r Input rotation
 * \param[out] t Input translation
 * \param[out] tx Input affine transformation
 * \return Output affine transformation
 *
 */

MESH_AFFINE mesh_transform_set_rotation_translation(MESH_ROTATION r, MESH_VECTOR3 t, MESH_AFFINE tx)
{
    if(tx==NULL) tx = mesh_affine_create();
    tx->data[0] = r->data[0];
    tx->data[1] = r->data[1];
    tx->data[2] = r->data[2];
    tx->data[3] = r->data[3];
    tx->data[4] = r->data[4];
    tx->data[5] = r->data[5];
    tx->data[6] = r->data[6];
    tx->data[7] = r->data[7];
    tx->data[8] = r->data[8];
    tx->data[9] = t->x;
    tx->data[10] = t->y;
    tx->data[11] = t->z;
    return tx;
}

/** \brief Sets an affine transformation with rotation and translation
 *
 * \param[in] m1 Input mesh
 * \param[out] m2 Input mesh to align
 * \param[in] flags (MESH_ALIGN_GLOBAL_POSITION/MESH_ALIGN_GLOBAL_ORIENTATION/MESH_ALIGN_GLOBAL_SCALE/MESH_ALIGN_GLOBAL_ALL/MESH_ALIGN_GLOBAL_DO_TRANSFORM)
 * \param[out] tx Output affine transformation, if not null
 * \return Error code
 *
 */

int mesh_align_global(MESH m1, MESH m2, int flags, MESH_AFFINE tx)
{

    return -1;
}

