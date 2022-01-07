/**
 * @file meshcreate.c
 * @author Sk. Mohammadul Haque
 * @version 1.8.0.0
 * @copyright
 * Copyright (c) 2013-2021 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to mesh creation and freeing.
 */

#include "../include/meshlib.h"
#include <string.h>

/** \brief Creates a new mesh
 *
 * \return Output mesh
 *
 */

MESH mesh_create_mesh_new()
{
    MESH m = NULL;
    if((m = (MESH)malloc(sizeof(mesh)))==NULL)
        mesh_error(MESH_ERR_MALLOC);
    m->is_vertices = 0;
    m->is_faces = 0;
    m->is_edges = 0;
    m->is_vcolors = 0;
    m->is_fcolors = 0;
    m->is_vnormals = 0;
    m->is_fnormals = 0;
    m->is_vfaces = 0;
    m->is_ffaces = 0;

    m->is_fareas = 0;
    m->is_vscalars = 0;
    m->is_fscalars = 0;

    m->num_vertices = 0;
    m->num_faces = 0;
    m->num_edges = 0;

    m->vertices = NULL;
    m->faces = NULL;
    m->edges = NULL;

    m->vcolors = NULL;
    m->fcolors = NULL;

    m->vnormals = NULL;
    m->fnormals = NULL;

    m->vfaces = NULL;
    m->ffaces = NULL;
    m->fareas = NULL;

    m->vscalars = NULL;
    m->fscalars = NULL;

    m->is_trimesh = 0;

    return m;
}

/** \brief Frees a mesh
 *
 * \param[in] m Input mesh
 * \return NULL
 *
 */

void mesh_free_mesh(MESH m)
{
    INTDATA i;
    if(m->is_vertices)
        free(m->vertices);
    if(m->is_faces)
    {
        for(i=0; i<m->num_faces; ++i)
            free(m->faces[i].vertices);
        free(m->faces);
    }
    if(m->is_edges)
        free(m->edges);
    if(m->is_vcolors)
        free(m->vcolors);
    if(m->is_fcolors)
        free(m->fcolors);
    if(m->is_vnormals)
        free(m->vnormals);
    if(m->is_fnormals)
        free(m->fnormals);
    if(m->is_vfaces)
    {
        for(i=0; i<m->num_vertices; ++i)
        {
            if(m->vfaces[i].faces!=NULL)
                free(m->vfaces[i].faces);
        }
        free(m->vfaces);
    }
    if(m->is_ffaces)
    {
        for(i=0; i<m->num_faces; ++i)
        {
            if(m->ffaces[i].faces!=NULL)
                free(m->ffaces[i].faces);
        }
        free(m->ffaces);
    }
    if(m->is_fareas)
        free(m->fareas);

    if(m->is_vscalars)
        free(m->vscalars);

    if(m->is_fscalars)
        free(m->fscalars);

    free(m);
}

/** \brief Creates a grid mesh
 *
 * \param[in] sz Size vector
 * \param[in] pos Position vector
 * \param[in] m Number of x-samples
 * \param[in] n Number of y-samples
 * \return Output mesh
 *
 */

MESH mesh_create_mesh_new_grid(MESH_VECTOR3 sz, MESH_VECTOR3 pos, INTDATA m, INTDATA n)
{
    MESH m1 = mesh_create_mesh_new();
    INTDATA mplus1 = m+1, nplus1 = n+1, nv = mplus1*nplus1, nf = 2*m*n;
    INTDATA nfby2 = m*n, i;
    FLOATDATA szx = -sz->x/2.0+pos->x, szy = -sz->y/2.0+pos->y, szz = pos->z;
    FLOATDATA dx = sz->x/m, dy = sz->y/n;
    MESH_VERTEX mv = NULL;
    MESH_FACE mf = NULL;

    if((m1->vertices = (MESH_VERTEX) malloc(nv*sizeof(mesh_vertex))) ==NULL)
        mesh_error(MESH_ERR_MALLOC);
    if((m1->faces = (MESH_FACE) malloc(nf*sizeof(mesh_face))) ==NULL)
        mesh_error(MESH_ERR_MALLOC);
    mv = m1->vertices;
    mf = m1->faces;
    m1->is_trimesh = 1;
    m1->num_vertices = nv;
    m1->is_vertices = 1;
    #pragma omp parallel for firstprivate(szx, szy, szz, dx, dy)
    for(i=0; i<nv; ++i)
    {
        MESH_VERTEX cv = mv+i;
        INTDATA y = i/mplus1;
        INTDATA x = i-y*mplus1;
        cv->x = szx+dx*x;
        cv->y = szy+dy*y;
        cv->z = szz;
    }

    m1->is_faces = 1;
    m1->num_faces = nf;
    #pragma omp parallel for firstprivate(m, n, mplus1, nplus1)
    for(i=0; i<nfby2; ++i)
    {
        INTDATA y = i/m, x = (i-m*y);
        MESH_FACE cf1 = mf+2*i, cf2 = mf+2*i+1;
        INTDATA xt = x%2, yt = y%2;
        INTDATA* cf1v = NULL, *cf2v = NULL;
        if((cf1->vertices = (INTDATA*) malloc(3*sizeof(INTDATA)))==NULL)
            mesh_error(MESH_ERR_MALLOC);
        if((cf2->vertices = (INTDATA*) malloc(3*sizeof(INTDATA)))==NULL)
            mesh_error(MESH_ERR_MALLOC);
        cf1v = cf1->vertices;
        cf2v = cf2->vertices;
        if(xt==yt)
        {
            /* first triangle */
            cf1->num_vertices = 3;
            cf1v[0] = x+mplus1*y;
            cf1v[2] = cf1v[0]+mplus1;
            cf1v[1] = cf1v[2]+1;

            /* second triangle */
            cf2->num_vertices = 3;
            cf2v[0] = cf1v[0];
            cf2v[1] = cf1v[0]+1;
            cf2v[2] = cf1v[1];
        }
        else
        {
            /* first triangle */
            cf1->num_vertices = 3;
            cf1v[0] = x+mplus1*y;
            cf1v[1] = cf1v[0]+1;
            cf1v[2] = cf1v[0]+mplus1;

            /* second triangle */
            cf2->num_vertices = 3;
            cf2v[0] = cf1v[1];
            cf2v[1] = cf1v[2]+1;
            cf2v[2] = cf1v[2];
        }
    }
    m1->is_loaded = 1;
    return m1;
}

/** \brief Creates a cuboid mesh
 *
 * \param[in] sz Size vector
 * \param[in] pos Position vector
 * \return Output mesh
 *
 */

MESH mesh_create_mesh_new_cuboid(MESH_VECTOR3 sz, MESH_VECTOR3 pos)
{
    MESH m = NULL;
    INTDATA i;
    m = mesh_create_mesh_new();

    m->num_vertices = 8;
    if((m->vertices = (MESH_VERTEX)malloc(8*sizeof(mesh_vertex)))==NULL)
        mesh_error(MESH_ERR_MALLOC);
    m->vertices[0].x = -sz->x*0.5;
    m->vertices[0].y = -sz->y*0.5;
    m->vertices[0].z = -sz->z*0.5;
    m->vertices[1].x = -sz->x*0.5;
    m->vertices[1].y = -sz->y*0.5;
    m->vertices[1].z = sz->z*0.5;
    m->vertices[2].x = -sz->x*0.5;
    m->vertices[2].y = sz->y*0.5;
    m->vertices[2].z = -sz->z*0.5;
    m->vertices[3].x = -sz->x*0.5;
    m->vertices[3].y = sz->y*0.5;
    m->vertices[3].z = sz->z*0.5;
    m->vertices[4].x = sz->x*0.5;
    m->vertices[4].y = -sz->y*0.5;
    m->vertices[4].z = -sz->z*0.5;
    m->vertices[5].x = sz->x*0.5;
    m->vertices[5].y = -sz->y*0.5;
    m->vertices[5].z = sz->z*0.5;
    m->vertices[6].x = sz->x*0.5;
    m->vertices[6].y = sz->y*0.5;
    m->vertices[6].z = -sz->z*0.5;
    m->vertices[7].x = sz->x*0.5;
    m->vertices[7].y = sz->y*0.5;
    m->vertices[7].z = sz->z*0.5;

    m->is_vertices = 1;
    m->is_trimesh = 1;
    m->is_loaded = 1;
    if((m->faces = (MESH_FACE)malloc(12*sizeof(mesh_face)))==NULL)
        mesh_error(MESH_ERR_MALLOC);
    m->num_faces = 12;
    for(i=0; i<m->num_faces; ++i)
    {
        m->faces[i].num_vertices = 3;
        if((m->faces[i].vertices = (INTDATA*)malloc(3*sizeof(INTDATA)))==NULL)
            mesh_error(MESH_ERR_MALLOC);
    }
    m->faces[0].vertices[0] = 0;
    m->faces[0].vertices[1] = 1;
    m->faces[0].vertices[2] = 2;
    m->faces[1].vertices[0] = 1;
    m->faces[1].vertices[1] = 3;
    m->faces[1].vertices[2] = 2;
    m->faces[2].vertices[0] = 0;
    m->faces[2].vertices[1] = 4;
    m->faces[2].vertices[2] = 1;
    m->faces[3].vertices[0] = 1;
    m->faces[3].vertices[1] = 4;
    m->faces[3].vertices[2] = 5;
    m->faces[4].vertices[0] = 0;
    m->faces[4].vertices[1] = 2;
    m->faces[4].vertices[2] = 4;
    m->faces[5].vertices[0] = 2;
    m->faces[5].vertices[1] = 6;
    m->faces[5].vertices[2] = 4;
    m->faces[6].vertices[0] = 1;
    m->faces[6].vertices[1] = 5;
    m->faces[6].vertices[2] = 3;
    m->faces[7].vertices[0] = 3;
    m->faces[7].vertices[1] = 5;
    m->faces[7].vertices[2] = 7;
    m->faces[8].vertices[0] = 2;
    m->faces[8].vertices[1] = 3;
    m->faces[8].vertices[2] = 6;
    m->faces[9].vertices[0] = 3;
    m->faces[9].vertices[1] = 7;
    m->faces[9].vertices[2] = 6;
    m->faces[10].vertices[0] = 4;
    m->faces[10].vertices[1] = 6;
    m->faces[10].vertices[2] = 5;
    m->faces[11].vertices[0] = 5;
    m->faces[11].vertices[1] = 6;
    m->faces[11].vertices[2] = 7;
    m->is_faces = 1;
    if(pos!=NULL)
        mesh_translate_vector(m, pos);
    return m;
}

/** \brief Creates an ellipsoid mesh
 *
 * \param[in] sz Size vector
 * \param[in] pos Position vector
 * \return Output mesh
 *
 */

MESH mesh_create_mesh_new_ellipsoid(MESH_VECTOR3 sz, MESH_VECTOR3 pos)
{
    MESH m = NULL;
    INTDATA i, j, k, vc, vnx, tmp, nhz = 36, nvr = 16;
    FLOATDATA a1, a2, sin1, sin2, cos1, cos2;
    m = mesh_create_mesh_new();
    m->num_vertices = nhz*nvr+2;
    if((m->vertices = (MESH_VERTEX)malloc((m->num_vertices)*sizeof(mesh_vertex)))==NULL)
        mesh_error(MESH_ERR_MALLOC);
    m->vertices[0].x = 0.0;
    m->vertices[0].y = sz->y*0.5;
    m->vertices[0].z = 0.0;
    for(j=0; j<nvr; ++j)
    {
        a1 = MESH_PI*(FLOATDATA)(j+1)/(nvr+1);
        sin1 = sin(a1);
        cos1 = cos(a1);
        for(i=0; i<nhz; ++i)
        {
            a2 = 2*MESH_PI*((FLOATDATA)i)/nhz;
            sin2 = sin(a2);
            cos2 = cos(a2);
            tmp = i+j*(nhz)+1;
            m->vertices[tmp].x = sin1*cos2*sz->x*0.5;
            m->vertices[tmp].y = cos1*sz->y*0.5;
            m->vertices[tmp].z = sin1*sin2*sz->z*0.5;
        }
    }
    m->vertices[m->num_vertices-1].x = 0.0;
    m->vertices[m->num_vertices-1].y = -sz->y*0.5;
    m->vertices[m->num_vertices-1].z = 0.0;

    m->is_vertices = 1;
    m->is_trimesh = 1;
    m->is_loaded = 1;
    m->num_faces = nhz*nvr*2;
    if((m->faces = (MESH_FACE)malloc(m->num_faces*sizeof(mesh_face)))==NULL)
        mesh_error(MESH_ERR_MALLOC);
    for(i=0; i<m->num_faces; ++i)
    {
        m->faces[i].num_vertices = 3;
        if((m->faces[i].vertices = (INTDATA*)malloc(3*sizeof(INTDATA)))==NULL)
            mesh_error(MESH_ERR_MALLOC);
    }
    k = 0;
    for(i=0; i<nhz-1; ++i)
    {
        m->faces[k].vertices[0] = 0;
        m->faces[k].vertices[1] = i+2;
        m->faces[k].vertices[2] = i+1;
        ++k;
    }
    m->faces[k].vertices[0] = 0;
    m->faces[k].vertices[1] = 1;
    m->faces[k].vertices[2] = nhz;
    ++k;
    for(j=0; j<nvr-1; ++j)
    {
        for(i=0; i<nhz-1; ++i)
        {
            vc = i+j*nhz+1;
            vnx = vc+nhz;
            m->faces[k].vertices[0] = vc;
            m->faces[k].vertices[1] = vc+1;
            m->faces[k].vertices[2] = vnx+1;
            ++k;
            m->faces[k].vertices[0] = vc;
            m->faces[k].vertices[1] = vnx+1;
            m->faces[k].vertices[2] = vnx;
            ++k;
        }
        vc = nhz+j*nhz;
        vnx = vc+nhz;
        m->faces[k].vertices[0] = vc;
        m->faces[k].vertices[1] = vc+1-nhz;
        m->faces[k].vertices[2] = vnx+1-nhz;
        ++k;
        m->faces[k].vertices[0] = vc;
        m->faces[k].vertices[1] = vnx+1-nhz;
        m->faces[k].vertices[2] = vnx;
        ++k;
    }
    for(i=0; i<nhz-1; ++i)
    {
        m->faces[k].vertices[0] = m->num_vertices-i-2;
        m->faces[k].vertices[1] = m->num_vertices-1;
        m->faces[k].vertices[2] = m->num_vertices-i-3;
        ++k;
    }
    m->faces[k].vertices[0] = m->num_vertices-nhz-1;
    m->faces[k].vertices[1] = m->num_vertices-1;
    m->faces[k].vertices[2] = m->num_vertices-2;

    m->is_faces = 1;
    if(pos!=NULL)
        mesh_translate_vector(m, pos);
    return m;
}

/** \brief Creates a uniformly sampled ellipsoid mesh
 *
 * \param[in] sz Size vector
 * \param[in] pos Position vector
 * \param[in] n Number of upsample ~3
 * \return Output mesh
 *
 */

MESH mesh_create_mesh_new_uniform_ellipsoid(MESH_VECTOR3 sz, MESH_VECTOR3 pos, int n)
{
    MESH m = NULL;
    INTDATA i;
    const int nv = 12, nf = 20;
    FLOATDATA q = (sqrt(5.0)+1.0)/2.0;
    const mesh_vertex verts[] = {{0, q, 1.0}, {0, q, -1.0}, {0, -q, 1.0}, {0, -q, -1.0},
        {q, 1.0, 0.0}, {q, -1.0, 0.0}, {-q, 1.0, 0.0}, {-q, -1.0, 0.0},
        {1.0, 0.0, q}, {-1.0, 0.0, q}, {1.0, 0.0,-q}, {-1.0, 0.0,-q}
    };
    const int faces[][3] =
    {
        {1,0,4},{0,1,6},{2,3,5},{3,2,7},
        {4,5,10},{5,4,8},{6,7,9},{7,6,11},
        {8,9,2},{9,8,0},{10,11,1},{11,10,3},
        {0,8,4},{0,6,9},{1,4,10},{1,11,6},
        {2,5,8},{2,9,7},{3,10,5},{3,7,11}
    };
    m = mesh_create_mesh_new();
    m->num_vertices = nv;
    if((m->vertices = (MESH_VERTEX)malloc((m->num_vertices)*sizeof(mesh_vertex)))==NULL)
        mesh_error(MESH_ERR_MALLOC);

    memcpy(m->vertices, verts, sizeof(mesh_vertex)*nv);
    m->is_vertices = 1;
    m->is_trimesh = 1;
    m->is_loaded = 1;
    m->num_faces = nf;
    if((m->faces = (MESH_FACE)malloc(m->num_faces*sizeof(mesh_face)))==NULL)
        mesh_error(MESH_ERR_MALLOC);
    for(i=0; i<m->num_faces; ++i)
    {
        const int *cface = faces[i];
        MESH_FACE cf = m->faces+i;
        if((cf->vertices = (INTDATA*)malloc(3*sizeof(INTDATA)))==NULL)
            mesh_error(MESH_ERR_MALLOC);
        cf->num_vertices = 3;

        cf->vertices[0] = cface[0];
        cf->vertices[1] = cface[1];
        cf->vertices[2] = cface[2];
    }
    m->is_faces = 1;
    mesh_upsample(m, n);
    for(i=0; i<m->num_vertices; ++i)
    {
        MESH_VERTEX cv = m->vertices+i;
        FLOATDATA l = sqrt(cv->x*cv->x+cv->y*cv->y+cv->z*cv->z);
        cv->x /= l;
        cv->y /= l;
        cv->z /= l;
    }
    mesh_scale(m, sz->x, sz->y, sz->z);
    if(pos!=NULL)
        mesh_translate_vector(m, pos);
    return m;
}

/** \brief Creates a cylinder mesh
 *
 * \param[in] sz Size vector
 * \param[in] pos Position vector
 * \return Output mesh
 *
 */

MESH mesh_create_mesh_new_cylinder(MESH_VECTOR3 sz, MESH_VECTOR3 pos)
{
    MESH m = NULL;
    INTDATA i, k, vc, nhz = 36;
    FLOATDATA a1, sin1, cos1;
    m = mesh_create_mesh_new();
    m->num_vertices = nhz*2+2;
    if((m->vertices = (MESH_VERTEX)malloc((m->num_vertices)*sizeof(mesh_vertex)))==NULL)
        mesh_error(MESH_ERR_MALLOC);
    m->vertices[0].x = 0.0;
    m->vertices[0].y = sz->y*0.5;
    m->vertices[0].z = 0.0;
    k = 1;
    for(i=0; i<nhz; ++i)
    {
        a1 = 2*MESH_PI*((FLOATDATA)i)/nhz;
        sin1 = sin(a1);
        cos1 = cos(a1);
        m->vertices[k].x = sin1*sz->x*0.5;
        m->vertices[k].y = sz->y*0.5;
        m->vertices[k].z = cos1*sz->z*0.5;
        ++k;
        m->vertices[k].x = sin1*sz->x*0.5;
        m->vertices[k].y = -sz->y*0.5;
        m->vertices[k].z = cos1*sz->z*0.5;
        ++k;
    }
    m->vertices[m->num_vertices-1].x = 0.0;
    m->vertices[m->num_vertices-1].y = -sz->y*0.5;
    m->vertices[m->num_vertices-1].z = 0.0;

    m->is_vertices = 1;
    m->is_trimesh = 1;
    m->is_loaded = 1;
    m->num_faces = nhz*4;
    if((m->faces = (MESH_FACE)malloc(m->num_faces*sizeof(mesh_face)))==NULL)
        mesh_error(MESH_ERR_MALLOC);
    for(i=0; i<m->num_faces; ++i)
    {
        m->faces[i].num_vertices = 3;
        if((m->faces[i].vertices = (INTDATA*)malloc(3*sizeof(INTDATA)))==NULL)
            mesh_error(MESH_ERR_MALLOC);
    }
    k = 0;
    for(i=0; i<nhz-1; ++i)
    {
        m->faces[k].vertices[0] = 0;
        m->faces[k].vertices[1] = 2*i+1;
        m->faces[k].vertices[2] = 2*i+3;
        ++k;
    }
    m->faces[k].vertices[0] = 0;
    m->faces[k].vertices[1] = 2*nhz-1;
    m->faces[k].vertices[2] = 1;
    ++k;
    for(i=0; i<nhz-1; ++i)
    {
        vc = 2*i+1;
        m->faces[k].vertices[0] = vc;
        m->faces[k].vertices[1] = vc+1;
        m->faces[k].vertices[2] = vc+2;
        ++k;
        m->faces[k].vertices[0] = vc+2;
        m->faces[k].vertices[1] = vc+1;
        m->faces[k].vertices[2] = vc+3;
        ++k;
    }
    vc = 2*nhz-1;
    m->faces[k].vertices[0] = vc;
    m->faces[k].vertices[1] = vc+1;
    m->faces[k].vertices[2] = vc+2-2*nhz;
    ++k;
    m->faces[k].vertices[0] = vc+2-2*nhz;
    m->faces[k].vertices[1] = vc+1;
    m->faces[k].vertices[2] = vc+3-2*nhz;
    ++k;

    for(i=0; i<nhz-1; ++i)
    {
        m->faces[k].vertices[0] = m->num_vertices-2*i-4;
        m->faces[k].vertices[1] = m->num_vertices-1;
        m->faces[k].vertices[2] = m->num_vertices-2*i-2;
        ++k;
    }
    m->faces[k].vertices[0] = m->num_vertices-2*nhz;
    m->faces[k].vertices[1] = m->num_vertices-2;
    m->faces[k].vertices[2] = m->num_vertices-1;

    m->is_faces = 1;
    if(pos!=NULL)
        mesh_translate_vector(m, pos);
    return m;
}

/** \brief Creates a cone mesh
 *
 * \param[in] sz Size vector
 * \param[in] pos Position vector
 * \return Output mesh
 *
 */

MESH mesh_create_mesh_new_cone(MESH_VECTOR3 sz, MESH_VECTOR3 pos)
{
    MESH m = NULL;
    INTDATA i, k, nhz = 36;
    FLOATDATA a1, sin1, cos1;
    m = mesh_create_mesh_new();
    m->num_vertices = nhz+2;
    if((m->vertices = (MESH_VERTEX)malloc((m->num_vertices)*sizeof(mesh_vertex)))==NULL)
        mesh_error(MESH_ERR_MALLOC);
    m->vertices[0].x = 0.0;
    m->vertices[0].y = sz->y*0.5;
    m->vertices[0].z = 0.0;
    k = 1;
    for(i=0; i<nhz; ++i)
    {
        a1 = 2*MESH_PI*((FLOATDATA)i)/nhz;
        sin1 = sin(a1);
        cos1 = cos(a1);
        m->vertices[k].x = sin1*sz->x*0.5;
        m->vertices[k].y = -sz->y*0.5;
        m->vertices[k].z = cos1*sz->z*0.5;
        ++k;
    }
    m->vertices[m->num_vertices-1].x = 0.0;
    m->vertices[m->num_vertices-1].y = -sz->y*0.5;
    m->vertices[m->num_vertices-1].z = 0.0;

    m->is_vertices = 1;
    m->is_trimesh = 1;
    m->is_loaded = 1;
    m->num_faces = nhz*2;
    if((m->faces = (MESH_FACE)malloc(m->num_faces*sizeof(mesh_face)))==NULL)
        mesh_error(MESH_ERR_MALLOC);
    for(i=0; i<m->num_faces; ++i)
    {
        m->faces[i].num_vertices = 3;
        if((m->faces[i].vertices = (INTDATA*)malloc(3*sizeof(INTDATA)))==NULL)
            mesh_error(MESH_ERR_MALLOC);
    }
    k = 0;
    for(i=1; i<nhz; ++i)
    {
        m->faces[k].vertices[0] = 0;
        m->faces[k].vertices[1] = i;
        m->faces[k].vertices[2] = i+1;
        ++k;
    }
    m->faces[k].vertices[0] = 0;
    m->faces[k].vertices[1] = nhz;
    m->faces[k].vertices[2] = 1;
    ++k;
    for(i=1; i<nhz; ++i)
    {
        m->faces[k].vertices[0] = m->num_vertices-i-2;
        m->faces[k].vertices[1] = m->num_vertices-1;
        m->faces[k].vertices[2] = m->num_vertices-i-1;
        ++k;
    }
    m->faces[k].vertices[0] = m->num_vertices-nhz-1;
    m->faces[k].vertices[1] = m->num_vertices-2;
    m->faces[k].vertices[2] = m->num_vertices-1;

    m->is_faces = 1;
    if(pos!=NULL)
        mesh_translate_vector(m, pos);
    return m;
}

/** \brief Creates a flat rectangle mesh
 *
 * \param[in] sz Size vector
 * \param[in] pos Position vector
 * \return Output mesh
 *
 */

MESH mesh_create_mesh_new_rectangle_flat(MESH_VECTOR3 sz, MESH_VECTOR3 pos)
{
    MESH m = mesh_create_mesh_new();
    INTDATA nv = 4, nf = 2;
    mesh_alloc_mesh_props(m, nv, nf, 0, MESH_PROPS_VERTICES|MESH_PROPS_FACES);
    MESH_VERTEX mv = m->vertices;
    mv[0].x = 0.5*sz->x+pos->x;
    mv[0].y = -0.5*sz->y+pos->y;
    mv[0].z = pos->z;

    mv[1].x = 0.5*sz->x+pos->x;
    mv[1].y = 0.5*sz->y+pos->y;
    mv[1].z = pos->z;

    mv[2].x = -0.5*sz->x+pos->x;
    mv[2].y = 0.5*sz->y+pos->y;
    mv[2].z = pos->z;

    mv[3].x = -0.5*sz->x+pos->x;
    mv[3].y = -0.5*sz->y+pos->y;
    mv[3].z = pos->z;


    MESH_FACE mf = m->faces;
    if((mf[0].vertices = (INTDATA*)malloc(3*sizeof(INTDATA)))==NULL)
    {
        mesh_error(MESH_ERR_MALLOC);
    }
    mf[0].num_vertices = 3;
    mf[0].vertices[0] = 0;
    mf[0].vertices[1] = 1;
    mf[0].vertices[2] = 2;

    if((mf[1].vertices = (INTDATA*)malloc(3*sizeof(INTDATA)))==NULL)
    {
        mesh_error(MESH_ERR_MALLOC);
    }
    mf[1].num_vertices = 3;
    mf[1].vertices[0] = 2;
    mf[1].vertices[1] = 3;
    mf[1].vertices[2] = 0;

    return m;
}

/** \brief Creates a flat ellipse mesh
 *
 * \param[in] sz Size vector
 * \param[in] pos Position vector
 * \return Output mesh
 *
 */

MESH mesh_create_mesh_new_ellipse_flat(MESH_VECTOR3 sz, MESH_VECTOR3 pos)
{
    MESH m = mesh_create_mesh_new();
    INTDATA nv = 37, nf = 36, i;
    mesh_alloc_mesh_props(m, nv, nf, 0, MESH_PROPS_VERTICES|MESH_PROPS_FACES);
    MESH_VERTEX mv = m->vertices;
    mv[0].x = pos->x;
    mv[0].y = pos->y;
    mv[0].z = pos->z;
    for(i=1; i<nv; ++i)
    {
        FLOATDATA th = MESH_TWOPI*(i-1)/((FLOATDATA)(nv-1.));
        mv[i].x = sz->x*cos(th)+pos->x;
        mv[i].y = sz->y*sin(th)+pos->y;
        mv[i].z = pos->z;
    }

    MESH_FACE mf = m->faces;
    for(i=0; i<nf; ++i)
    {
        INTDATA d1 = i+1;
        INTDATA d2 = i+2;
        if(i==(nf-1))
        {
            d2 = 1;
        }
        if((mf[i].vertices = (INTDATA*)malloc(3*sizeof(INTDATA)))==NULL)
        {
            mesh_error(MESH_ERR_MALLOC);
        }
        mf[i].num_vertices = 3;
        mf[i].vertices[0] = 0;
        mf[i].vertices[1] = d1;
        mf[i].vertices[2] = d2;
    }
    return m;
}


