/**
 * @file meshops.c
 * @author Sk. Mohammadul Haque
 * @version 1.4.2.0
 * @copyright
 * Copyright (c) 2013, 2014, 2015, 2016 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to mesh combinatorial operations.
 */

#include <string.h>
#include <omp.h>
#include "../include/meshlib.h"


/** \brief Clones a given mesh into another mesh
 *
 * \param[in] m Input mesh to clone
 * \param[in] flags Flags to copy which properties (MESH_CLONE_VERTICES/MESH_CLONE_VNORMALS/MESH_CLONE_VCOLORS/MESH_CLONE_VFACES/MESH_CLONE_V_ALL_PROPS/MESH_CLONE_FACES/MESH_CLONE_FNORMALS/MESH_CLONE_FCOLORS/MESH_CLONE_FAREAS/MESH_CLONE_F_ALL_PROPS/MESH_CLONE_ALL_PROPS)
 * \return Output cloned mesh
 *
 */

MESH mesh_clone_mesh(MESH m, uint16_t flags)
{
    MESH m2 = NULL;
    INTDATA i;
    m2 = mesh_create_mesh_new();
    m2->origin_type = m->origin_type;
    m2->is_loaded = m->is_loaded;
    if((flags&MESH_CLONE_VERTICES)&&m->is_vertices)
    {
        if((m2->vertices = (MESH_VERTEX)malloc(sizeof(mesh_vertex)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy(m2->vertices, m->vertices, sizeof(mesh_vertex)*(m->num_vertices));
        m2->is_vertices = m->is_vertices;
        m2->num_vertices = m->num_vertices;
    }

    if((flags&__MESH_CLONE_VNORMALS)&&m->is_vnormals)
    {
        if((m2->vnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy(m2->vnormals, m->vnormals, sizeof(mesh_normal)*(m->num_vertices));
        m2->is_vnormals = m->is_vnormals;
    }

    if((flags&__MESH_CLONE_VCOLORS)&&m->is_vcolors)
    {
        if((m2->vcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy(m2->vcolors, m->vcolors, sizeof(mesh_color)*(m->num_vertices));
        m2->is_vcolors = m->is_vcolors;
    }

    if((flags&__MESH_CLONE_VFACES)&&m->is_vfaces)
    {
        if((m2->vfaces = (MESH_VFACE)malloc(sizeof(mesh_vface)*m->num_vertices))==NULL) mesh_error(MESH_ERR_MALLOC);
        #pragma omp parallel for shared(m, m2)
        for(i=0; i<m->num_vertices; ++i)
        {
            if((m2->vfaces[i].faces = (INTDATA *)malloc(sizeof(INTDATA)*(m->vfaces[i].num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
            memcpy(m2->vfaces[i].faces, m->vfaces[i].faces, sizeof(INTDATA)*(m->vfaces[i].num_faces));
            m2->vfaces[i].num_faces = m->vfaces[i].num_faces;
        }
    }

    if((flags&__MESH_CLONE_FACES)&&m->is_faces)
    {
        if((m2->faces = (MESH_FACE)malloc(sizeof(mesh_face)*m->num_faces))==NULL) mesh_error(MESH_ERR_MALLOC);
        #pragma omp parallel for shared(m, m2)
        for(i=0; i<m->num_faces; ++i)
        {
            if((m2->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*(m->faces[i].num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
            memcpy(m2->faces[i].vertices, m->faces[i].vertices, sizeof(INTDATA)*(m->faces[i].num_vertices));
            m2->faces[i].num_vertices = m->faces[i].num_vertices;
        }
        m2->is_faces = m->is_faces;
        m2->num_faces = m->num_faces;
        m2->is_trimesh = m->is_trimesh;
    }

    if((flags&__MESH_CLONE_FNORMALS)&&m->is_fnormals)
    {
        if((m2->fnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy(m2->fnormals, m->fnormals, sizeof(mesh_normal)*(m->num_faces));
        m2->is_fnormals = m->is_fnormals;
    }

    if((flags&__MESH_CLONE_FCOLORS)&&m->is_fcolors)
    {
        if((m2->fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy(m2->fcolors, m->fcolors, sizeof(mesh_color)*(m->num_faces));
        m2->is_fcolors = m->is_fcolors;
    }

    if((flags&__MESH_CLONE_FFACES)&&m->is_ffaces)
    {
        if((m2->ffaces = (MESH_FFACE)malloc(sizeof(mesh_fface)*m->num_faces))==NULL) mesh_error(MESH_ERR_MALLOC);
        #pragma omp parallel for shared(m, m2)
        for(i=0; i<m->num_faces; ++i)
        {
            if((m2->ffaces[i].faces = (INTDATA *)malloc(sizeof(INTDATA)*(m->ffaces[i].num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
            memcpy(m2->ffaces[i].faces, m->ffaces[i].faces, sizeof(INTDATA)*(m->ffaces[i].num_faces));
            m2->ffaces[i].num_faces = m->ffaces[i].num_faces;
        }
    }

    if((flags&__MESH_CLONE_FAREAS)&&m->is_fareas)
    {
        if((m2->fareas = (FLOATDATA*)malloc(sizeof(FLOATDATA)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy(m2->fareas, m->fareas, sizeof(FLOATDATA)*(m->num_faces));
        m2->is_fareas = m->is_fareas;
    }

    if((flags&__MESH_CLONE_EDGES)&&m->is_edges)
    {
        if((m2->edges = (MESH_EDGE)malloc(sizeof(mesh_edge)*(m->num_edges)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy(m2->edges, m->edges, sizeof(mesh_edge)*(m->num_edges));
        m2->is_edges = m->is_edges;
        m2->num_edges = m->num_edges;
    }

    return m2;
}

/** \brief Combines a given mesh with another given mesh
 *
 * \param[in] m1 Input mesh to combine with
 * \param[in] m2 Input mesh to combine
 * \return Output combined mesh
 *
 */

MESH mesh_combine_mesh(MESH m1, MESH m2)
{
    INTDATA i, m1_nv, m1_nf, m1_ne;
    if(m1==NULL) return mesh_clone_mesh(m2, MESH_CLONE_ALL_PROPS);
    if(m1->is_vertices>m2->is_vertices) return m1;
    if(m1->is_vertices<m2->is_vertices)
    {
        mesh_free_mesh(m1);
        return mesh_clone_mesh(m2, MESH_CLONE_ALL_PROPS);
    }
    m1->is_loaded |= m2->is_loaded;
    m1_nv = m1->num_vertices;
    m1_nf = m1->num_faces;
    m1_ne = m1->num_edges;

    if((m1->is_vnormals>m2->is_vnormals)||(m1->is_vcolors>m2->is_vcolors)|| (m1->is_vfaces>m2->is_vfaces)|| (m1->is_faces>m2->is_faces)||(m1->is_fnormals>m2->is_fnormals)||(m1->is_fcolors>m2->is_fcolors)|| (m1->is_ffaces>m2->is_ffaces)|| (m1->is_fareas>m2->is_fareas)|| (m1->is_edges>m2->is_edges)) mesh_error(MESH_ERR_INCOMPATIBLE);

    if(m1->is_vertices)
    {
        if((m1->vertices = (MESH_VERTEX)realloc(m1->vertices, sizeof(mesh_vertex)*(m1_nv+m2->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy((m1->vertices)+m1_nv, m2->vertices, sizeof(mesh_vertex)*(m2->num_vertices));
        m1->num_vertices += m2->num_vertices;
    }

    if(m1->is_vnormals)
    {
        if((m1->vnormals = (MESH_NORMAL)realloc(m1->vnormals, sizeof(mesh_normal)*(m1_nv+m2->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy((m1->vnormals)+m1_nv, m2->vnormals, sizeof(mesh_normal)*(m2->num_vertices));
    }

    if(m1->is_vcolors)
    {
        if((m1->vcolors = (MESH_COLOR)realloc(m1->vcolors, sizeof(mesh_color)*(m1_nv+m2->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy((m1->vcolors)+m1_nv, m2->vcolors, sizeof(mesh_color)*(m2->num_vertices));
    }

    if(m1->is_vfaces)
    {
        if((m1->vfaces = (MESH_VFACE)realloc(m1->vfaces, sizeof(mesh_vface)*(m1_nv+m2->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        #pragma omp parallel for shared(m1, m2, m1_nv, m1_nf)
        for(i=0; i<m2->num_vertices; ++i)
        {
            INTDATA j;
            if((m1->vfaces[m1_nv+i].faces = (INTDATA *)malloc(sizeof(INTDATA)*(m2->vfaces[i].num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
            m1->vfaces[m1_nv+i].num_faces = m2->vfaces[i].num_faces;
            for(j=0; j<m2->vfaces[i].num_faces; ++j)
            {
                m1->vfaces[m1_nv+i].faces[j] = m2->vfaces[i].faces[j]+m1_nf;
            }
        }
    }

    if(m1->is_faces)
    {
        if((m1->faces = (MESH_FACE)realloc(m1->faces, sizeof(mesh_face)*(m1_nf+m2->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        #pragma omp parallel for shared(m1, m2, m1_nv, m1_nf)
        for(i=0; i<m2->num_faces; ++i)
        {
            INTDATA j;
            if((m1->faces[m1_nf+i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*(m2->faces[i].num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
            m1->faces[m1_nf+i].num_vertices = m2->faces[i].num_vertices;
            for(j=0; j<m2->faces[i].num_vertices; ++j)
            {
                m1->faces[m1_nf+i].vertices[j] = m2->faces[i].vertices[j]+m1_nv;
            }
        }
        m1->num_faces += m2->num_faces;
        m1->is_trimesh = m1->is_trimesh && m2->is_trimesh;
    }

    if(m1->is_fnormals)
    {
        if((m1->fnormals = (MESH_NORMAL)realloc(m1->fnormals, sizeof(mesh_normal)*(m1_nf+m2->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy((m1->fnormals)+m1_nf, m2->fnormals, sizeof(mesh_normal)*(m2->num_faces));
    }

    if(m1->is_fcolors)
    {
        if((m1->fcolors = (MESH_COLOR)realloc(m1->fcolors, sizeof(mesh_color)*(m1_nf+m2->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy((m1->fcolors)+m1_nf, m2->fcolors, sizeof(mesh_color)*(m2->num_faces));
    }

    if(m1->is_ffaces)
    {
        if((m1->ffaces = (MESH_FFACE)realloc(m1->ffaces, sizeof(mesh_fface)*(m1_nf+m2->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        #pragma omp parallel for shared(m1, m2, m1_nv, m1_nf)
        for(i=0; i<m2->num_faces; ++i)
        {
            INTDATA j;
            if((m1->ffaces[m1_nf+i].faces = (INTDATA *)malloc(sizeof(INTDATA)*(m2->ffaces[i].num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
            m1->ffaces[m1_nf+i].num_faces = m2->ffaces[i].num_faces;
            for(j=0; j<m2->ffaces[i].num_faces; ++j)
            {
                m1->ffaces[m1_nf+i].faces[j] = m2->ffaces[i].faces[j]+m1_nf;
            }
        }
    }

    if(m2->is_fareas)
    {
        if((m1->fareas = (FLOATDATA*)realloc(m1->fareas, sizeof(FLOATDATA)*(m1_nf+m2->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memcpy((m1->fareas)+m1_nf, m2->fareas, sizeof(FLOATDATA)*(m2->num_faces));
    }

     if(m1->is_edges)
    {
        if((m1->edges = (MESH_EDGE)realloc(m1->edges, sizeof(mesh_edge)*(m1_ne+m2->num_edges)))==NULL) mesh_error(MESH_ERR_MALLOC);
        #pragma omp parallel for shared(m1, m2, m1_nv, m1_nf, m1_ne)
        for(i=0; i<m2->num_edges; ++i)
        {
            m1->edges[m1_ne+i].vertices[0] = m2->edges[i].vertices[0]+m1_nv;
            m1->edges[m1_ne+i].vertices[1] = m2->edges[i].vertices[1]+m1_nv;
            m1->edges[m1_ne+i].faces[0] = m2->edges[i].faces[0]+m1_nf;
            m1->edges[m1_ne+i].faces[1] = m2->edges[i].faces[1]+m1_nf;
        }
    }
    return m1;
}

