/**
 * @file meshclean.c
 * @author Sk. Mohammadul Haque
 * @version 1.4.2.0
 * @copyright
 * Copyright (c) 2013, 2014, 2015, 2016 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to different mesh cleaning algorithms.
 */

#include <string.h>
#include "../include/meshlib.h"

/** \cond HIDDEN_SYMBOLS */
static __inline INTDATA __mesh_find_parent(INTDATA i, INTDATA* parents);

static __inline FLOATDATA __mesh_calc_triangle_area(MESH_VERTEX a, MESH_VERTEX b, MESH_VERTEX c)
{
    mesh_vertex p, q, r;
    FLOATDATA area;
    p.x = a->x - b->x;
    p.y = a->y - b->y;
    p.z = a->z - b->z;
    q.x = a->x - c->x;
    q.y = a->y - c->y;
    q.z = a->z - c->z;
    mesh_cross_vector3(&p, &q, &r);
    area = 0.5*sqrt(r.x*r.x+r.y*r.y+r.z*r.z);
    return area;
}

static __inline FLOATDATA __mesh_calc_vertex_distance_squared(MESH_VERTEX a, MESH_VERTEX b)
{
    FLOATDATA dx, dy, dz;
    dx = a->x-b->x;
    dy = a->y-b->y;
    dz = a->z-b->z;
    return (dx*dx+dy*dy+dz*dz);
}

INTDATA __mesh_find_parent(INTDATA i, INTDATA* parents)
{
    if(parents[i]!=i) parents[i] = __mesh_find_parent(parents[i], parents);
    return parents[i];
}

static __inline void __mesh_union(INTDATA i, INTDATA j, INTDATA* parents, INTDATA* ranks)
{
    INTDATA ir, jr;
    ir = __mesh_find_parent(i, parents);
    jr = __mesh_find_parent(j, parents);
    if(ir==jr) return;
    if(ranks[i]<ranks[j]) parents[i] = j;
    else if(ranks[i]>ranks[j]) parents[j] = i;
    else
    {
        parents[j] = i;
        ++ranks[i];
    }
}

static __inline INTDATA __mesh_find2(MESH_STRUCT2 s, INTDATA q)
{
    INTDATA k;
    for(k=0; k<s->num_items; ++k)
    {
        if(s->items[k][0]==q) return k;
    }
    return -1;
}

#define __mesh_rm_vertices (0)
#define __mesh_rm_faces (1)

static int __mesh_remove_boundary_elements(MESH m, int iters, int type)
{
    MESH_FACE new_faces = NULL;
    MESH_COLOR new_fcolors = NULL;
    MESH_NORMAL new_fnormals = NULL;
    signed char *fflags = NULL;
    INTDATA i, j, k, s, num_deleted;
    uint8_t is_ffaces;

    if(m==NULL) return 1;
    if(m->is_faces==0) return 2;
    if(m->is_trimesh==0) return 3;

    for(s=0; s<iters; ++s)
    {
        num_deleted = 0;
        if(m->is_edges!=1) mesh_calc_edges(m);
        if((fflags = (signed char *)malloc(sizeof(char)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memset(fflags, 0, sizeof(char)*(m->num_faces));


        if(type==__mesh_rm_vertices)
        {
            for(i=0; i<m->num_edges; ++i)
            {
                if(m->edges[i].faces[0]<0||m->edges[i].faces[1]<0)
                {
                    for(k=0; k<m->vfaces[m->edges[i].vertices[0]].num_faces; ++k)
                    {
                        fflags[m->vfaces[m->edges[i].vertices[0]].faces[k]]= 1;
                    }
                    for(k=0; k<m->vfaces[m->edges[i].vertices[1]].num_faces; ++k)
                    {
                        fflags[m->vfaces[m->edges[i].vertices[1]].faces[k]]= 1;
                    }
                }
            }
        }
        else
        {
            for(i=0; i<m->num_edges; ++i)
            {
                if(m->edges[i].faces[0]<0||m->edges[i].faces[1]<0)
                {
                    MESH_FACE curr_face;
                    for(k=0; k<m->vfaces[m->edges[i].vertices[0]].num_faces; ++k)
                    {
                        curr_face = &(m->faces[m->vfaces[m->edges[i].vertices[0]].faces[k]]);
                        if(curr_face->vertices[0]==m->edges[i].vertices[1]||curr_face->vertices[1]==m->edges[i].vertices[1]||curr_face->vertices[2]==m->edges[i].vertices[1])
                        {
                            fflags[m->vfaces[m->edges[i].vertices[0]].faces[k]]= 1;
                            break;
                        }
                    }
                }
            }
        }

        for(i=0; i<m->num_faces; ++i) if(fflags[i]==1) ++num_deleted;
        if(num_deleted>0)
        {
            free(m->edges);
            m->num_edges = 0;
            m->is_edges = 0;
            if(m->is_ffaces)
            {
                is_ffaces = 1;
                #pragma omp parallel for shared(m)
                for(i=0; i<m->num_faces; ++i)
                {
                    if(m->ffaces[i].faces!=NULL) free(m->ffaces[i].faces);
                }
                free(m->ffaces);
                m->ffaces = NULL;
            }
            m->is_ffaces = 0;

            if((new_faces = (MESH_FACE)malloc(sizeof(mesh_face)*(m->num_faces-num_deleted)))==NULL) mesh_error(MESH_ERR_MALLOC);
            j = 0;
            for(i=0; i<m->num_faces; ++i)
            {
                if(fflags[i]!=1)
                {
                    new_faces[j].num_vertices = 3;
                    if((new_faces[j].vertices = (INTDATA *)malloc(sizeof(INTDATA)*3))==NULL) mesh_error(MESH_ERR_MALLOC);
                    new_faces[j].vertices[0] = m->faces[i].vertices[0];
                    new_faces[j].vertices[1] = m->faces[i].vertices[1];
                    new_faces[j].vertices[2] = m->faces[i].vertices[2];
				    ++j;
                }
                free(m->faces[i].vertices);
            }
            if(m->is_fcolors)
            {
                if((new_fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_faces-num_deleted)))==NULL) mesh_error(MESH_ERR_MALLOC);
                j = 0;
                for(i=0; i<m->num_faces; ++i)
                {
                    if(fflags[i]!=1)
                    {
                        new_fcolors[j].r = m->fcolors[i].r;
                        new_fcolors[j].g = m->fcolors[i].g;
                        new_fcolors[j].b = m->fcolors[i].b;
                        new_fcolors[j].a = m->fcolors[i].a;
                        ++j;
                    }
                }
                free(m->fcolors);
                m->fcolors = new_fcolors;
            }

            if(m->is_fnormals)
            {
                if((new_fnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_faces-num_deleted)))==NULL) mesh_error(MESH_ERR_MALLOC);
                j = 0;
                for(i=0; i<m->num_faces; ++i)
                {
                    if(fflags[i]!=1)
                    {
                        new_fnormals[j].x = m->fnormals[i].x;
                        new_fnormals[j].y = m->fnormals[i].y;
                        new_fnormals[j].z = m->fnormals[i].z;
                    }
                }
                free(m->fnormals);
                m->fnormals = new_fnormals;
            }

            if(m->is_vfaces)
            {
                for(i=0; i<m->num_vertices; ++i)
                {
                    if(m->vfaces[i].faces!=NULL) free(m->vfaces[i].faces);
                }
                free(m->vfaces);
                m->vfaces = NULL;
            }
            m->is_vfaces = 0;

            m->num_faces -= num_deleted;
            free(m->faces);
            m->faces = new_faces;
            mesh_remove_unreferenced_vertices(m);
            if(is_ffaces==1) mesh_calc_face_adjacency(m);

            mesh_calc_edges(m);
        }
        free(fflags);
    }
    return 0;
}

/** \endcond */

/** \brief Removes boundary vertices and connecting elements
 *
 * \param[in] m Input mesh
 * \param[in] iters Number of iterations
 * \return Error code
 *
 */

int mesh_remove_boundary_vertices(MESH m, int iters)
{
    return __mesh_remove_boundary_elements(m, iters, __mesh_rm_vertices);
}

/** \brief Removes boundary faces and connecting elements
 *
 * \param[in] m Input mesh
 * \param[in] iters Number of iterations
 * \return Error code
 *
 */

int mesh_remove_boundary_faces(MESH m, int iters)
{
    return __mesh_remove_boundary_elements(m, iters, __mesh_rm_faces);
}

/** \brief Removes triangles with area smaller than a given value
 *
 * \param[in] m Input mesh
 * \param[in] area Given area
 * \return Error code
 *
 */

int mesh_remove_triangles_with_small_area(MESH m, FLOATDATA area)
{
    signed char *fflags = NULL;
    MESH_FACE new_faces = NULL;
    MESH_COLOR new_fcolors = NULL;
    MESH_NORMAL new_fnormals = NULL;
    INTDATA i, j, k, num_deleted = 0;
    uint8_t is_ffaces = 0;
    if(area==0) area = 1e-18;
    if(m==NULL) return 1;
    if(m->is_faces==0) return 2;
    if(m->is_trimesh==0) return 3;

    if((fflags = (signed char *)malloc(sizeof(signed char)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
    memset(fflags, 0, sizeof(signed char)*(m->num_faces));
    for(i=0; i<m->num_faces; ++i)
    {
        if(mesh_calc_triangle_area(&(m->vertices[m->faces[i].vertices[0]]),&(m->vertices[m->faces[i].vertices[1]]),&(m->vertices[m->faces[i].vertices[2]]))<area)
        {
            fflags[i]= 1;
            ++num_deleted;
        }
    }
    if(num_deleted>0)
    {
        if(m->is_ffaces)
        {
            is_ffaces = 1;
            #pragma omp parallel for shared(m)
            for(i=0; i<m->num_faces; ++i)
            {
                if(m->ffaces[i].faces!=NULL) free(m->ffaces[i].faces);
            }
            free(m->ffaces);
            m->ffaces = NULL;
        }
        if(m->is_edges)
        {
            k = 0;
            for(i=0; i<m->num_edges; ++i)
            {
                if(fflags[m->edges[i].faces[0]]==0 &&fflags[m->edges[i].faces[1]]==0)
                {
                    m->edges[k].vertices[0] = m->edges[i].vertices[0];
                    m->edges[k].vertices[1] = m->edges[i].vertices[1];
                    m->edges[k].faces[0] = m->edges[i].faces[0];
                    m->edges[k].faces[1] = m->edges[i].faces[1];
                    ++k;
                }
            }
            if((m->edges = (MESH_EDGE)realloc(m->edges,sizeof(mesh_edge)*k))==NULL) mesh_error(MESH_ERR_MALLOC);
        }
        if(m->is_vfaces)
        {
            for(i=0; i<m->num_vertices; ++i)
            {
                k = 0;
                for(j=0; j<m->vfaces[i].num_faces; ++j)
                {
                    if(fflags[m->vfaces[i].faces[j]]!=1)
                    {
                        m->vfaces[i].faces[k] = m->vfaces[i].faces[j];
                        ++k;
                    }
                }
                m->vfaces[i].num_faces = k;
                if((m->vfaces[i].faces = (INTDATA*) realloc(m->vfaces[i].faces, sizeof(INTDATA)*k))==NULL && k!=0) mesh_error(MESH_ERR_MALLOC);
            }
        }
        if((new_faces = (MESH_FACE)malloc(sizeof(mesh_face)*(m->num_faces-num_deleted)))==NULL) mesh_error(MESH_ERR_MALLOC);
        j = 0;
        for(i=0; i<m->num_faces; ++i)
        {
            if(fflags[i]!=1)
            {
                new_faces[j].num_vertices = 3;
                if((new_faces[j].vertices = (INTDATA *)malloc(sizeof(INTDATA)*3))==NULL) mesh_error(MESH_ERR_MALLOC);
                new_faces[j].vertices[0] = m->faces[i].vertices[0];
                new_faces[j].vertices[1] = m->faces[i].vertices[1];
                new_faces[j].vertices[2] = m->faces[i].vertices[2];
			    ++j;
            }
            free(m->faces[i].vertices);
        }
        if(m->is_fcolors)
        {
            if((new_fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_faces-num_deleted)))==NULL) mesh_error(MESH_ERR_MALLOC);
            j = 0;
            for(i=0; i<m->num_faces; ++i)
            {
                if(fflags[i]!=1)
                {
                    new_fcolors[j].r = m->fcolors[i].r;
                    new_fcolors[j].g = m->fcolors[i].g;
                    new_fcolors[j].b = m->fcolors[i].b;
                    new_fcolors[j].a = m->fcolors[i].a;
                    ++j;
                }
            }
            free(m->fcolors);
            m->fcolors = new_fcolors;
        }

        if(m->is_fnormals)
        {
            if((new_fnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_faces-num_deleted)))==NULL) mesh_error(MESH_ERR_MALLOC);
            for(i=0; i<m->num_faces; ++i)
            {
                if(fflags[i]!=1)
                {
                    new_fnormals[j].x = m->fnormals[i].x;
                    new_fnormals[j].y = m->fnormals[i].y;
                    new_fnormals[j].z = m->fnormals[i].z;
                }
            }
            free(m->fnormals);
            m->fnormals = new_fnormals;
        }

        m->num_faces -= num_deleted;
        free(m->faces);
        m->faces = new_faces;
        free(fflags);
        mesh_calc_vertex_adjacency(m);
        if(is_ffaces==1) mesh_calc_face_adjacency(m);
    }
    return 0;
}

/** \brief Removes triangles with zero area
 *
 * \param[in] m Input mesh
 * \return Error code
 *
 */

int mesh_remove_zero_area_faces(MESH m)
{
    return mesh_remove_triangles_with_small_area(m, 0.0);
}

/** \brief Removes unreferenced vertices
 *
 * \param[in] m Input mesh
 * \return Error code
 *
 */

int mesh_remove_unreferenced_vertices(MESH m)
{
    signed char *vflags = NULL;
    INTDATA *vindx = NULL;
    MESH_VERTEX new_vertices = NULL;
    MESH_COLOR new_vcolors = NULL;
    MESH_NORMAL new_vnormals = NULL;
    INTDATA num_valid_flags = 0, i, j;
    uint8_t is_vfaces, is_edges;
    if(m==NULL) return 1;

    if((vflags = (signed char *)malloc(sizeof(signed char)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    if((vindx = (INTDATA *)malloc(sizeof(INTDATA)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    memset(vflags, 0, sizeof(signed char)*(m->num_vertices));
    memset(vindx, 0, sizeof(INTDATA)*(m->num_vertices));

    for(i=0; i<m->num_faces; ++i)
    {
        for(j=0; j<m->faces[i].num_vertices; ++j)
        {
            vflags[m->faces[i].vertices[j]] = 1;
        }
    }
    vindx[0] = vflags[0]-1;
    for(i=1; i<m->num_vertices; ++i)
    {
        vindx[i] = vindx[i-1]+vflags[i];
    }
    num_valid_flags = vindx[m->num_vertices-1]+1;
    if(num_valid_flags<m->num_vertices)
    {
        if(m->is_vfaces)
        {
            is_vfaces = 1;
            #pragma omp parallel for shared(m)
            for(i=0; i<m->num_vertices; ++i)
            {
                if(m->vfaces[i].faces!=NULL) free(m->vfaces[i].faces);
            }
            free(m->vfaces);
            m->vfaces = NULL;
            m->is_vfaces = 0;
        }
        if(m->is_edges)
        {
            is_edges = 1;
            free(m->edges);
            m->edges = NULL;
            m->num_edges = 0;
        }

        for(i=0; i<m->num_faces; ++i)
        {
            for(j=0; j<m->faces[i].num_vertices; ++j)
            {
                m->faces[i].vertices[j] = vindx[m->faces[i].vertices[j]];
            }
        }

        if((new_vertices = (MESH_VERTEX)malloc(sizeof(mesh_vertex)*(num_valid_flags)))==NULL) mesh_error(MESH_ERR_MALLOC);
        for(i=0; i<m->num_vertices; ++i)
        {
            if(vflags[i]==1)
            {
                new_vertices[vindx[i]].x = m->vertices[i].x;
                new_vertices[vindx[i]].y = m->vertices[i].y;
                new_vertices[vindx[i]].z = m->vertices[i].z;
            }
        }
        if(m->is_vcolors)
        {
            if((new_vcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(num_valid_flags)))==NULL) mesh_error(MESH_ERR_MALLOC);
            for(i=0; i<m->num_vertices; ++i)
            {
                if(vflags[i]==1)
                {
                    new_vcolors[vindx[i]].r = m->vcolors[i].r;
                    new_vcolors[vindx[i]].g = m->vcolors[i].g;
                    new_vcolors[vindx[i]].b = m->vcolors[i].b;
                    new_vcolors[vindx[i]].a = m->vcolors[i].a;
                }
            }
            free(m->vcolors);
            m->vcolors = new_vcolors;
        }

        if(m->is_vnormals)
        {
            if((new_vnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(num_valid_flags)))==NULL) mesh_error(MESH_ERR_MALLOC);
            for(i=0; i<m->num_vertices; ++i)
            {
                if(vflags[i]==1)
                {
                    new_vnormals[vindx[i]].x = m->vnormals[i].x;
                    new_vnormals[vindx[i]].y = m->vnormals[i].y;
                    new_vnormals[vindx[i]].z = m->vnormals[i].z;
                }
            }
            free(m->vnormals);
            m->vnormals = new_vnormals;
        }

        m->num_vertices = num_valid_flags;
        free(m->vertices);
        m->vertices = new_vertices;
        if(is_vfaces==1) mesh_calc_vertex_adjacency(m);
        if(is_edges==1) mesh_calc_edges(m);
    }
    free(vflags);
    free(vindx);

    return 0;
}

/** \brief Removes ear faces and connecting vertices
 *
 * \param[in] m Input mesh
 * \param[in] niters Number of iterations
 * \return Error code
 *
 */

int mesh_remove_ear_faces(MESH m, int niters)
{
    signed char *fflags = NULL;
    MESH_FACE new_faces = NULL;
    MESH_COLOR new_fcolors = NULL;
    MESH_NORMAL new_fnormals = NULL;
    INTDATA i, j, num_deleted;
    uint8_t is_ffaces, is_edges;
    int iters;
    if(m==NULL) return 1;
    if(m->is_faces==0) return 2;
    if(m->is_ffaces)
    {
        is_ffaces = 1;
        #pragma omp parallel for shared(m)
        for(i=0; i<m->num_faces; ++i)
        {
            if(m->ffaces[i].faces!=NULL) free(m->ffaces[i].faces);
        }
        free(m->ffaces);
        m->ffaces = NULL;
    }
    if(m->is_edges)
    {
        is_edges = 1;
        free(m->edges);
        m->edges = NULL;
        m->num_edges = 0;
    }
    if(!m->is_vfaces) mesh_calc_vertex_adjacency(m);

    for(iters=0; iters<niters; ++iters)
    {
        if((fflags = (signed char *)malloc(sizeof(signed char)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        memset(fflags, 0, sizeof(signed char)*(m->num_faces));
        num_deleted = 0;
        for(i=0; i<m->num_vertices; ++i)
        {
            if(m->vfaces[i].num_faces==1 && fflags[m->vfaces[i].faces[0]]==0)
            {
                fflags[m->vfaces[i].faces[0]]= 1;
                ++num_deleted;
            }
        }
        if(num_deleted>0)
        {
            if((new_faces = (MESH_FACE)malloc(sizeof(mesh_face)*(m->num_faces-num_deleted)))==NULL) mesh_error(MESH_ERR_MALLOC);
            j = 0;
            for(i=0; i<m->num_faces; ++i)
            {
                if(fflags[i]!=1)
                {
                    new_faces[j].num_vertices = 3;
                    if((new_faces[j].vertices = (INTDATA *)malloc(sizeof(INTDATA)*3))==NULL) mesh_error(MESH_ERR_MALLOC);
                    new_faces[j].vertices[0] = m->faces[i].vertices[0];
                    new_faces[j].vertices[1] = m->faces[i].vertices[1];
                    new_faces[j].vertices[2] = m->faces[i].vertices[2];
					++j;
                }
               	free(m->faces[i].vertices);
            }

            if(m->is_fcolors)
            {
                if((new_fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_faces-num_deleted)))==NULL) mesh_error(MESH_ERR_MALLOC);
                j = 0;
                for(i=0; i<m->num_faces; ++i)
                {
                    if(fflags[i]!=1)
                    {
                        new_fcolors[j].r = m->fcolors[i].r;
                        new_fcolors[j].g = m->fcolors[i].g;
                        new_fcolors[j].b = m->fcolors[i].b;
                        new_fcolors[j].a = m->fcolors[i].a;
                        ++j;
                    }
                }
                free(m->fcolors);
                m->fcolors = new_fcolors;
            }

            if(m->is_fnormals)
            {
                if((new_fnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_faces-num_deleted)))==NULL) mesh_error(MESH_ERR_MALLOC);
                for(i=0; i<m->num_faces; ++i)
                {
                    if(fflags[i]!=1)
                    {
                        new_fnormals[j].x = m->fnormals[i].x;
                        new_fnormals[j].y = m->fnormals[i].y;
                        new_fnormals[j].z = m->fnormals[i].z;
                    }
                }
                free(m->fnormals);
                m->fnormals = new_fnormals;
            }

            if(m->is_vfaces)
            {
                for(i=0; i<m->num_vertices; ++i)
                {
                    if(m->vfaces[i].faces!=NULL) free(m->vfaces[i].faces);
                }
                free(m->vfaces);
                m->vfaces = NULL;
            }
            m->is_vfaces = 0;

            m->num_faces -= num_deleted;
            free(m->faces);
            m->faces = new_faces;
            free(fflags);
            mesh_calc_vertex_adjacency(m);
        }
        else break;
    }
    mesh_remove_unreferenced_vertices(m);
    if(is_ffaces==1) mesh_calc_face_adjacency(m);
    if(is_edges==1) mesh_calc_edges(m);
    return 0;

}

/** \brief Removes close vertices
 *
 * \param[in] m Input mesh
 * \param[in] r Maximum distance between two vertices
 * \return Error code
 *
 */

int mesh_remove_close_vertices(MESH m, FLOATDATA r)
{
    INTDATA i, j, k;
    INTDATA *vparents = NULL;
    INTDATA *vranks = NULL;
    if(!m->is_vfaces) mesh_calc_vertex_adjacency(m);
    if((vparents = (INTDATA *)malloc(sizeof(INTDATA)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    if((vranks = (INTDATA *)malloc(sizeof(INTDATA)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    r *= r;
    for(i=0; i<m->num_vertices; ++i)
    {
        vparents[i] = i;
        vranks[i] = 0;
    }
    for(i=0; i<m->num_vertices; ++i)
    {
        for(j=0; j<m->vfaces[i].num_faces; ++j)
        {
            for(k=0; k<m->faces[m->vfaces[i].faces[j]].num_vertices; ++k)
            {
                if(i==m->faces[m->vfaces[i].faces[j]].vertices[k]) continue;
                if(__mesh_calc_vertex_distance_squared(&(m->vertices[i]), &(m->vertices[m->faces[m->vfaces[i].faces[j]].vertices[k]]))<r)
                {
                    __mesh_union(i, m->faces[m->vfaces[i].faces[j]].vertices[k], vparents, vranks);
                }
            }
        }
    }

    for(i=0; i<m->num_faces; ++i)
    {
        for(j=0; j<m->faces[i].num_vertices; ++j)
        {
            m->faces[i].vertices[j] = __mesh_find_parent(m->faces[i].vertices[j], vparents);
        }
    }
    free(vparents);
    free(vranks);
    mesh_remove_zero_area_faces(m);
    mesh_remove_unreferenced_vertices(m);
    return 0;
}

/** \brief Removes non-manifold vertices
 *
 * \param[in] m Input mesh
 * \return Error code
 *
 */

int mesh_remove_non_manifold_vertices(MESH m)
{
    signed char *vflags = NULL, *fflags = NULL;
    INTDATA *vindx = NULL;
    FLOATDATA* new_fareas = NULL;
    MESH_VERTEX new_vertices = NULL;
    MESH_FACE new_faces = NULL;
    MESH_COLOR new_vcolors = NULL, new_fcolors = NULL;
    MESH_NORMAL new_vnormals = NULL, new_fnormals = NULL;
    MESH_EDGE me;
    INTDATA mne, mnv;
    INTDATA num_valid_flags = 0, i, j, k;
    uint8_t is_vfaces = 0, is_edges = 0, is_ffaces = 0;
    if(m==NULL) return 1;
    if(m->is_edges==0) mesh_calc_edges(m);
    else is_edges = 1;

    if((vflags = (signed char *)malloc(sizeof(signed char)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    if((vindx = (INTDATA *)malloc(sizeof(INTDATA)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    memset(vflags, 0, sizeof(signed char)*(m->num_vertices));
    memset(vindx, 0, sizeof(INTDATA)*(m->num_vertices));

    me = m->edges;
    mne = m->num_edges;
    mnv = m->num_vertices;
    for(i=0; i<mne; ++i)
    {
        MESH_EDGE ce = me+i;
        if(ce->faces[0]<0||ce->faces[1]<0)
        {
            ++vflags[ce->vertices[0]];
            ++vflags[ce->vertices[1]];
        }
    }
    free(m->edges);
    m->edges = NULL;
    m->num_edges = 0;
    m->is_edges = 0;

    for(i=0; i<mnv; ++i)
    {
        if(vflags[i]>2) vflags[i] = 0;
        else vflags[i] = 1;
    }

    vindx[0] = vflags[0]-1;
    for(i=1; i<m->num_vertices; ++i)
    {
        vindx[i] = vindx[i-1]+vflags[i];
    }
    num_valid_flags = vindx[m->num_vertices-1]+1;
    if(num_valid_flags<m->num_vertices)
    {
        if(m->is_vfaces)
        {
            is_vfaces = 1;
            #pragma omp parallel for shared(m)
            for(i=0; i<m->num_vertices; ++i)
            {
                if(m->vfaces[i].faces!=NULL) free(m->vfaces[i].faces);
            }
            free(m->vfaces);
            m->vfaces = NULL;
            m->is_vfaces = 0;
        }

        if((new_vertices = (MESH_VERTEX)malloc(sizeof(mesh_vertex)*(num_valid_flags)))==NULL) mesh_error(MESH_ERR_MALLOC);
        for(i=0; i<m->num_vertices; ++i)
        {
            if(vflags[i]==1)
            {
                new_vertices[vindx[i]].x = m->vertices[i].x;
                new_vertices[vindx[i]].y = m->vertices[i].y;
                new_vertices[vindx[i]].z = m->vertices[i].z;
            }
        }

        if(m->is_vcolors)
        {
            if((new_vcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(num_valid_flags)))==NULL) mesh_error(MESH_ERR_MALLOC);
            for(i=0; i<m->num_vertices; ++i)
            {
                if(vflags[i]==1)
                {
                    new_vcolors[vindx[i]].r = m->vcolors[i].r;
                    new_vcolors[vindx[i]].g = m->vcolors[i].g;
                    new_vcolors[vindx[i]].b = m->vcolors[i].b;
                    new_vcolors[vindx[i]].a = m->vcolors[i].a;
                }
            }
            free(m->vcolors);
            m->vcolors = new_vcolors;
        }

        if(m->is_vnormals)
        {
            if((new_vnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(num_valid_flags)))==NULL) mesh_error(MESH_ERR_MALLOC);
            for(i=0; i<m->num_vertices; ++i)
            {
                if(vflags[i]==1)
                {
                    new_vnormals[vindx[i]].x = m->vnormals[i].x;
                    new_vnormals[vindx[i]].y = m->vnormals[i].y;
                    new_vnormals[vindx[i]].z = m->vnormals[i].z;
                }
            }
            free(m->vnormals);
            m->vnormals = new_vnormals;
        }

        m->num_vertices = num_valid_flags;
        free(m->vertices);
        m->vertices = new_vertices;

        if(m->is_faces)
        {
            if((fflags = (signed char *)malloc(sizeof(signed char)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
            memset(fflags, 1, sizeof(signed char)*(m->num_faces));
            num_valid_flags = m->num_faces;
            for(i=0; i<m->num_faces; ++i)
            {
                INTDATA* cv = m->faces[i].vertices;
                INTDATA cvn = m->faces[i].num_vertices;
                for(j=0; j<cvn; ++j)
                {
                    if(vflags[cv[j]]==0)
                    {
                        fflags[i] = 0;
                        --num_valid_flags;
                        break;
                    }
                }
            }

            if((new_faces = (MESH_FACE) malloc(sizeof(mesh_face)*(num_valid_flags)))==NULL) mesh_error(MESH_ERR_MALLOC);
            k = 0;
            for(i=0; i<m->num_faces; ++i)
            {
                MESH_FACE cfo = m->faces+i;
                if(fflags[i]==1)
                {
                    MESH_FACE cf = new_faces+k;

                    if((cf->vertices = (INTDATA*) malloc(sizeof(INTDATA)*(cfo->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
                    cf->num_vertices = cfo->num_vertices;
                    for(j=0; j<cf->num_vertices; ++j)
                    {
                        new_faces[k].vertices[j] = vindx[m->faces[i].vertices[j]];
                    }
					++k;
                }
                free(cfo->vertices);
            }
            free(m->faces);
            m->faces = new_faces;

            if(m->is_fcolors)
            {
                if((new_fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(num_valid_flags)))==NULL) mesh_error(MESH_ERR_MALLOC);
                k = 0;
                for(i=0; i<m->num_faces; ++i)
                {
                    if(fflags[i]==1)
                    {
                        new_fcolors[k].r = m->fcolors[i].r;
                        new_fcolors[k].g = m->fcolors[i].g;
                        new_fcolors[k].b = m->fcolors[i].b;
                        new_fcolors[k].a = m->fcolors[i].a;
                        ++k;
                    }
                }
                free(m->fcolors);
                m->fcolors = new_fcolors;
            }

            if(m->is_fnormals)
            {
                if((new_fnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(num_valid_flags)))==NULL) mesh_error(MESH_ERR_MALLOC);
                k = 0;
                for(i=0; i<m->num_faces; ++i)
                {
                    if(fflags[i]==1)
                    {
                        new_fnormals[k].x = m->fnormals[i].x;
                        new_fnormals[k].y = m->fnormals[i].y;
                        new_fnormals[k].z = m->fnormals[i].z;
                        ++k;
                    }
                }
                free(m->fnormals);
                m->fnormals = new_fnormals;
            }

            if(m->is_fareas)
            {
                if((new_fareas = (FLOATDATA*)malloc(sizeof(FLOATDATA)*(num_valid_flags)))==NULL) mesh_error(MESH_ERR_MALLOC);
                k = 0;
                for(i=0; i<m->num_faces; ++i)
                {
                    if(fflags[i]==1)
                    {
                        new_fareas[k] = m->fareas[i];
                        ++k;
                    }
                }
                free(m->fareas);
                m->fareas = new_fareas;
            }

            if(m->is_ffaces)
            {
                is_ffaces = 1;
                #pragma omp parallel for shared(m)
                for(i=0; i<m->num_faces; ++i)
                {
                    if(m->ffaces[i].faces!=NULL) free(m->ffaces[i].faces);
                }
                free(m->ffaces);
                m->ffaces = NULL;
                m->is_ffaces = 0;
            }

            m->num_faces = num_valid_flags;
            free(fflags);
        }

        if(is_vfaces==1) mesh_calc_vertex_adjacency(m);
        if(is_ffaces==1) mesh_calc_face_adjacency(m);
        if(is_edges==1) mesh_calc_edges(m);
    }
    free(vflags);
    free(vindx);

    return 0;
}


