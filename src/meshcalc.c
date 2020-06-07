/**
 * @file meshcalc.c
 * @author Sk. Mohammadul Haque
 * @version 1.4.2.0
 * @copyright
 * Copyright (c) 2013, 2014, 2015, 2016 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to different mesh computations.
 */

#include <string.h>
#include <omp.h>
#include "../include/meshlib.h"

/** \cond HIDDEN_SYMBOLS */

static __inline void __mesh_cross(MESH_NORMAL x, MESH_NORMAL y, MESH_NORMAL z)
{
    /* FLOATDATA n; */
    z->x = x->y*y->z - y->y*x->z;
    z->y = x->z*y->x - y->z*x->x;
    z->z = x->x*y->y - y->x*x->y;
    /*    n = sqrt(z->x*z->x+z->y*z->y+z->z*z->z);
        if(n>0)
        {
            z->x /=n;
            z->y /=n;
            z->z /=n;
        }*/
}

static __inline INTDATA __mesh_find(MESH_STRUCT s, INTDATA q)
{
    INTDATA k;
    for(k=0; k<s->num_items; ++k)
    {
        if(s->items[k]==q) return k;
    }
    return -1;
}

static __inline INTDATA __mesh_find3(MESH_STRUCT3 s, INTDATA q)
{
    INTDATA k;
    for(k=0; k<s->num_items; ++k)
    {
        if(s->items[k][0]==q) return k;
    }
    return -1;
}

static __inline FLOATDATA __mesh_calc_vertex_distance_squared(MESH_VERTEX a, MESH_VERTEX b)
{
    FLOATDATA dx, dy, dz;
    dx = a->x-b->x;
    dy = a->y-b->y;
    dz = a->z-b->z;
    return (dx*dx+dy*dy+dz*dz);
}

#ifndef __mesh_isnan
#define __mesh_isnan(x) ((x)!=(x))
#endif


/** \endcond */

/** \brief Computes the cross product of two 3-d vectors
 *
 * \param[in] x First vector
 * \param[in] y Second vector
 * \param[out] z Output cross product \f$ \mathbf{x}\times \mathbf{y} \f$
 * \return NULL
 *
 */

void mesh_cross_vector3(MESH_VECTOR3 x, MESH_VECTOR3 y, MESH_VECTOR3 z)
{
    z->x = x->y*y->z - y->y*x->z;
    z->y = x->z*y->x - y->z*x->x;
    z->z = x->x*y->y - y->x*x->y;
}

/** \brief Computes the normalized cross product of two normals
 *
 * \param[in] x First normal
 * \param[in] y Second normal
 * \param[out] z Output cross product \f$ \frac{\mathbf{x}\times \mathbf{y}}{\|\mathbf{x}\times \mathbf{y}\|_2} \f$
 * \return NULL
 *
 */

void mesh_cross_normal(MESH_NORMAL x, MESH_NORMAL y, MESH_NORMAL z)
{
    FLOATDATA n;
    z->x = x->y*y->z - y->y*x->z;
    z->y = x->z*y->x - y->z*x->x;
    z->z = x->x*y->y - y->x*x->y;
    n = sqrt(z->x*z->x+z->y*z->y+z->z*z->z);
    if(n>0.0)
    {
        z->x /=n;
        z->y /=n;
        z->z /=n;
    }
}

/** \brief Computes the face normal given 3 vertices
 *
 * \param[in] v1 First vertex
 * \param[in] v2 Second vertex
 * \param[in] v3 Third vertex
 * \param[out] n Output face normal \f$ \mathbf{n}_f \f$
 * \return NULL
 *
 */

void mesh_calc_face_normal(MESH_VERTEX v1, MESH_VERTEX v2, MESH_VERTEX v3, MESH_NORMAL n)
{
    FLOATDATA qx, qy, qz, px, py, pz, t;
    px = v2->x-v1->x;
    py = v2->y-v1->y;
    pz = v2->z-v1->z;
    qx = v3->x-v1->x;
    qy = v3->y-v1->y;
    qz = v3->z-v1->z;
    n->x = py*qz - pz*qy;
    n->y = pz*qx - px*qz;
    n->z = px*qy - py*qx;
    t = sqrt(n->x*n->x+n->y*n->y+n->z*n->z);
    if(t>0.0)
    {
        n->x /= t;
        n->y /= t;
        n->z /= t;
    }
}

/** \brief Computes vertex normals of a given mesh
 *
 * \param[in] m Input mesh
 * \return Error code
 *
 */

int mesh_calc_vertex_normals(MESH m)
{
    INTDATA i;

    if(m==NULL) return 1;
    if(m->is_faces==0) return 2;
    if(m->vfaces==0) mesh_calc_vertex_adjacency(m);
    if(!m->is_vnormals)
    {
        if((m->vnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    }
    memset(m->vnormals, 0, sizeof(mesh_normal)*(m->num_vertices));
    m->is_vnormals = 1;

    #pragma omp parallel for shared(m)
    for(i=0; i<m->num_faces; ++i)
    {
        INTDATA in0, in1, in2;
        FLOATDATA l0, l1, l2, l01, l12, l20;
        mesh_normal e0, e1, e2, e3;
        in0 = m->faces[i].vertices[0];
        in1 = m->faces[i].vertices[1];
        in2 = m->faces[i].vertices[2];
        e2.x = m->vertices[in0].x - m->vertices[in1].x;
        e2.y = m->vertices[in0].y - m->vertices[in1].y;
        e2.z = m->vertices[in0].z - m->vertices[in1].z;

        e0.x = m->vertices[in1].x - m->vertices[in2].x;
        e0.y = m->vertices[in1].y - m->vertices[in2].y;
        e0.z = m->vertices[in1].z - m->vertices[in2].z;

        e1.x = m->vertices[in2].x - m->vertices[in0].x;
        e1.y = m->vertices[in2].y - m->vertices[in0].y;
        e1.z = m->vertices[in2].z - m->vertices[in0].z;

        __mesh_cross(&e2, &e0, &e3);

        l0 = e0.x*e0.x+e0.y*e0.y+e0.z*e0.z;
        l1 = e1.x*e1.x+e1.y*e1.y+e1.z*e1.z;
        l2 = e2.x*e2.x+e2.y*e2.y+e2.z*e2.z;
        l01 = 1.0/(l0*l1);
        l12 = 1.0/(l1*l2);
        l20 = 1.0/(l2*l0);

        #pragma omp critical
        {
            m->vnormals[in0].x +=e3.x*l12;
            m->vnormals[in0].y +=e3.y*l12;
            m->vnormals[in0].z +=e3.z*l12;

            m->vnormals[in1].x +=e3.x*l20;
            m->vnormals[in1].y +=e3.y*l20;
            m->vnormals[in1].z +=e3.z*l20;

            m->vnormals[in2].x +=e3.x*l01;
            m->vnormals[in2].y +=e3.y*l01;
            m->vnormals[in2].z +=e3.z*l01;
        }
    }
    #pragma omp parallel for
    for(i=0; i<m->num_vertices; ++i)
    {
        FLOATDATA t = sqrt(m->vnormals[i].x*m->vnormals[i].x+m->vnormals[i].y*m->vnormals[i].y+m->vnormals[i].z*m->vnormals[i].z);
		if(__mesh_isnan(t))
		{
			m->vnormals[i].x = 0.0;
            m->vnormals[i].y = 0.0;
            m->vnormals[i].z = 0.0;
		}
        if(t>0.0)
        {
            m->vnormals[i].x /= t;
            m->vnormals[i].y /= t;
            m->vnormals[i].z /= t;
        }
    }
    return 0;
}

/** \brief Computes face normals of a given mesh
 *
 * \param[in] m Input mesh
 * \return Error code
 *
 */

int mesh_calc_face_normals(MESH m)
{
    INTDATA i;
    MESH_VERTEX v1, v2, v3;
    FLOATDATA qx, qy, qz, px, py, pz, t;

    if(m==NULL) return 1;
    if(m->is_faces==0) return 2;
    if(!m->is_fnormals)
    {
        if((m->fnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
    }
    m->is_fnormals = 1;

    #pragma omp parallel for shared(m) private(v1, v2, v3, qx, qy, qz, px, py, pz, t)
    for(i=0; i<m->num_faces; ++i)
    {
        v1 = &(m->vertices[m->faces[i].vertices[0]]);
        v2 = &(m->vertices[m->faces[i].vertices[1]]);
        v3 = &(m->vertices[m->faces[i].vertices[2]]);
        px = v2->x-v1->x;
        py = v2->y-v1->y;
        pz = v2->z-v1->z;
        qx = v3->x-v1->x;
        qy = v3->y-v1->y;
        qz = v3->z-v1->z;
        m->fnormals[i].x = py*qz - pz*qy;
        m->fnormals[i].y = pz*qx - px*qz;
        m->fnormals[i].z = px*qy - py*qx;
        t = sqrt((m->fnormals[i].x)*(m->fnormals[i].x)+(m->fnormals[i].y)*(m->fnormals[i].y)+(m->fnormals[i].z)*(m->fnormals[i].z));
        if(t>0.0)
        {
            m->fnormals[i].x /= t;
            m->fnormals[i].y /= t;
            m->fnormals[i].z /= t;
        }
    }
    return 0;
}

/** \brief Computes edges of a given mesh
 *
 * \param[in] m Input mesh
 * \return Error code
 *
 */

int mesh_calc_edges(MESH m)
{
    MESH_STRUCT3 e_table = NULL;
    int fg_tmp;
    INTDATA i, j, num_edges = 0;
    INTDATA i_01, i_12, i_02;
    INTDATA v0, v1, v2, v_tmp;
    if(m->is_edges) free(m->edges);
    m->num_edges = 0;
    m->edges = NULL;
    if(m==NULL) return 1;
    if(m->is_faces==0) return 2;

    if(m->is_vfaces!=1) mesh_calc_vertex_adjacency(m);

    e_table = (MESH_STRUCT3)malloc(m->num_vertices*sizeof(mesh_struct3));
    if(e_table==NULL) mesh_error(MESH_ERR_MALLOC);
    for(i=0; i<m->num_vertices; ++i)
    {
        e_table[i].num_items = 0;
        e_table[i].items = NULL;
    }
    for(i=0; i<m->num_faces; ++i)
    {
        fg_tmp = 0;

        v0 = m->faces[i].vertices[0];
        v1 = m->faces[i].vertices[1];
        v2 = m->faces[i].vertices[2];

        if(v0>v1)
        {
            v_tmp = v0;
            v0 = v1;
            v1 = v_tmp;
            fg_tmp = 1-fg_tmp;
        }
        if(v0>v2)
        {
            v_tmp = v0;
            v0 = v2;
            v2 = v_tmp;
            fg_tmp = 1-fg_tmp;
        }
        if(v1>v2)
        {
            v_tmp = v1;    /* v0<v1<v2 */
            v1 = v2;
            v2 = v_tmp;
            fg_tmp = 1-fg_tmp;
        }

        i_01 = __mesh_find3(&e_table[v0], v1);

        if(i_01<0)
        {
            if((e_table[v0].items = (INTDATA3*)realloc(e_table[v0].items, sizeof(INTDATA3)*(e_table[v0].num_items+1)))==NULL) mesh_error(MESH_ERR_MALLOC);
            e_table[v0].num_items += 1;
            e_table[v0].items[e_table[v0].num_items-1][0] = v1;

            e_table[v0].items[e_table[v0].num_items-1][1+fg_tmp] = i;
            e_table[v0].items[e_table[v0].num_items-1][2-fg_tmp] = -1;
            ++num_edges;
        }
        else
        {
            if(e_table[v0].items[i_01][2-fg_tmp]>=0) e_table[v0].items[i_01][1+fg_tmp] = i;
            else e_table[v0].items[i_01][2-fg_tmp] = i;
        }

        i_12 = __mesh_find3(&e_table[v1], v2);

        if(i_12<0)
        {
            if((e_table[v1].items = (INTDATA3*)realloc(e_table[v1].items, sizeof(INTDATA3)*(e_table[v1].num_items+1)))==NULL) mesh_error(MESH_ERR_MALLOC);
            e_table[v1].num_items += 1;
            e_table[v1].items[e_table[v1].num_items-1][0] = v2;

            e_table[v1].items[e_table[v1].num_items-1][1+fg_tmp] = i;
            e_table[v1].items[e_table[v1].num_items-1][2-fg_tmp] = -1;
            ++num_edges;
        }
        else
        {
            if(e_table[v1].items[i_12][2-fg_tmp]>=0) e_table[v1].items[i_12][1+fg_tmp] = i;
            else e_table[v1].items[i_12][2-fg_tmp] = i;
        }

        i_02 = __mesh_find3(&e_table[v0], v2);

        if(i_02<0)
        {
            if((e_table[v0].items = (INTDATA3*)realloc(e_table[v0].items, sizeof(INTDATA3)*(e_table[v0].num_items+1)))==NULL) mesh_error(MESH_ERR_MALLOC);
            e_table[v0].num_items += 1;
            e_table[v0].items[e_table[v0].num_items-1][0] = v2;

            e_table[v0].items[e_table[v0].num_items-1][2-fg_tmp] = i;
            e_table[v0].items[e_table[v0].num_items-1][1+fg_tmp] = -1;
            ++num_edges;
        }
        else
        {
            if(e_table[v0].items[i_02][1+fg_tmp]>=0) e_table[v0].items[i_02][2-fg_tmp] = i;
            else e_table[v0].items[i_02][1+fg_tmp] = i;
        }
    }

    m->num_edges = num_edges;
    if((m->edges = (MESH_EDGE)malloc(num_edges*sizeof(mesh_edge)))==NULL) mesh_error(MESH_ERR_MALLOC);
    num_edges = 0;
    for(i=0; i<m->num_vertices; ++i)
    {
        for(j=0; j<e_table[i].num_items; ++j)
        {
            m->edges[num_edges].vertices[0] = i;
            m->edges[num_edges].vertices[1] = e_table[i].items[j][0];
            m->edges[num_edges].faces[0] = e_table[i].items[j][1];
            m->edges[num_edges].faces[1] = e_table[i].items[j][2];
            ++num_edges;
        }
        free(e_table[i].items);
    }
    free(e_table);
    m->is_edges = 1;
    return 0;
}

/** \brief Computes vertex adjacent faces of a given mesh
 *
 * \param[in] m Input mesh
 * \return Error code
 *
 */

int mesh_calc_vertex_adjacency(MESH m)
{
    INTDATA i, j;
    if(m==NULL) return 1;
    if(m->is_faces==0) return 2;
    if(m->is_vfaces)
    {
        #pragma omp parallel for shared(m)
        for(i=0; i<m->num_vertices; ++i)
        {
            if(m->vfaces[i].faces!=NULL) free(m->vfaces[i].faces);
        }
        free(m->vfaces);
        m->vfaces = NULL;
    }
    m->vfaces = (MESH_VFACE)malloc(m->num_vertices*sizeof(mesh_vface));
    if(m->vfaces==NULL) mesh_error(MESH_ERR_MALLOC);
    #pragma omp parallel for shared(m)
    for(i=0; i<m->num_vertices; ++i)
    {
        m->vfaces[i].num_faces = 0;
        m->vfaces[i].faces = NULL;
    }
    m->is_vfaces = 1;
    for(i=0; i<m->num_faces; ++i)
    {
        for(j=0; j<m->faces[i].num_vertices; ++j)
        {
            m->vfaces[m->faces[i].vertices[j]].faces = (INTDATA*)realloc(m->vfaces[m->faces[i].vertices[j]].faces, sizeof(INTDATA)*(m->vfaces[m->faces[i].vertices[j]].num_faces+1));
            ++(m->vfaces[m->faces[i].vertices[j]].num_faces);
            m->vfaces[m->faces[i].vertices[j]].faces[m->vfaces[m->faces[i].vertices[j]].num_faces-1] = i;
        }
    }
    return 0;
}

/** \brief Computes face adjacent faces of a given mesh
 *
 * \param[in] m Input mesh
 * \return Error code
 *
 */

int mesh_calc_face_adjacency(MESH m)
{
    INTDATA i;
    INTDATA* vidx = NULL;
    if(m==NULL) return 1;
    if(m->is_faces==0) return 2;
    if(m->is_edges==0) mesh_calc_edges(m);
    if(m->is_ffaces)
    {
        #pragma omp parallel for shared(m)
        for(i=0; i<m->num_faces; ++i)
        {
            if(m->ffaces[i].faces!=NULL) free(m->ffaces[i].faces);
            m->ffaces[i].num_faces = 0;
        }
    }
    else
    {
        m->ffaces = (MESH_FFACE)malloc(m->num_faces*sizeof(mesh_fface));
        if(m->ffaces==NULL) mesh_error(MESH_ERR_MALLOC);
    }
    m->is_ffaces = 1;

    if((vidx = (INTDATA*) malloc((m->num_vertices+1)*sizeof(INTDATA)))==NULL) mesh_error(MESH_ERR_MALLOC);
    memset(vidx, 0, (m->num_vertices+1)*sizeof(INTDATA));
    for(i=0; i<m->num_edges; ++i)
    {
        ++vidx[m->edges[i].vertices[0]+1];
    }

    for(i=0; i<m->num_vertices; ++i)
    {
        vidx[i+1] += vidx[i];
    }

    #pragma omp parallel for shared(m)
    for(i=0; i<m->num_faces; ++i)
    {
        INTDATA k = 0, v0, v1, l, cnt = 0, kk;
        m->ffaces[i].faces = NULL;

        if((m->ffaces[i].faces = (INTDATA*)malloc(sizeof(INTDATA)*m->faces[i].num_vertices))==NULL) mesh_error(MESH_ERR_MALLOC);

        for(k=0; k<m->faces[i].num_vertices; ++k)
        {
            kk = k+1;
            if(kk>=m->faces[i].num_vertices) kk = 0;
            if(m->faces[i].vertices[k]<m->faces[i].vertices[kk])
            {
                v0 = m->faces[i].vertices[k];
                v1 = m->faces[i].vertices[kk];
            }
            else
            {
                v0 = m->faces[i].vertices[kk];
                v1 = m->faces[i].vertices[k]; /* v0<v1 */
            }
            if(vidx[v0]>=0)
            {
                for(l=vidx[v0]; l<vidx[v0+1]; ++l)
                {
                    if(m->edges[l].vertices[1]==v1)
                    {
                        if(m->edges[l].faces[0]==i && m->edges[l].faces[1]>=0)
                        {
                            m->ffaces[i].faces[cnt] = m->edges[l].faces[1];
                            ++cnt;
                            break;
                        }
                        else if(m->edges[l].faces[0]>=0)
                        {
                            m->ffaces[i].faces[cnt] = m->edges[l].faces[0];
                            ++cnt;
                            break;
                        }
                    }
                }
            }
        }
        if((m->ffaces[i].faces = (INTDATA*)realloc(m->ffaces[i].faces, sizeof(INTDATA)*cnt))==NULL) mesh_error(MESH_ERR_MALLOC);
        m->ffaces[i].num_faces = cnt;
    }
    free(vidx);
    return 0;
}


/** \brief Finds an item in an INTDATA structure
 *
 * \param[in] s Input INTDATA structure
 * \param[in] q Query INTDATA
 * \return Index or -1
 *
 */

INTDATA mesh_find(MESH_STRUCT s, INTDATA q)
{
    INTDATA k;
    for(k=0; k<s->num_items; ++k)
    {
        if(s->items[k]==q) return k;
    }
    return -1;
}

/** \brief Finds an item in an INTDATA2 structure
 *
 * \param[in] s Input INTDATA2 structure
 * \param[in] q Query INTDATA2
 * \return Index or -1
 *
 */

INTDATA mesh_find2(MESH_STRUCT2 s, INTDATA q)
{
    INTDATA k;
    for(k=0; k<s->num_items; ++k)
    {
        if(s->items[k][0]==q) return k;
    }
    return -1;
}

/** \brief Finds an item in an INTDATA3 structure
 *
 * \param[in] s Input INTDATA3 structure
 * \param[in] q Query INTDATA3
 * \return Index or -1
 *
 */

INTDATA mesh_find3(MESH_STRUCT3 s, INTDATA q)
{
    INTDATA k;
    for(k=0; k<s->num_items; ++k)
    {
        if(s->items[k][0]==q) return k;
    }
    return -1;
}

/** \brief Upsamples a given mesh
 *
 * \param[in] m Input mesh
 * \param[in] iters Number of iterations
 * \return Error code
 *
 */

int mesh_upsample(MESH m, int iters)
{
    MESH_STRUCT v_table = NULL;
    MESH_FACE new_faces = NULL;
    INTDATA i, k, curr_idx;
    FLOATDATA t;
    INTDATA kv, i_01, i_12, i_20, c_v_0_indcs, c_v_1_indcs, c_v_2_indcs;
    uint8_t is_vfaces, is_ffaces, is_edges;

    if(m==NULL) return 1;
    if(m->is_faces==0) return 2;
    if(m->is_trimesh==0) return 3;
    if(m->is_fcolors!=0) return 4;
    if(m->is_fnormals!=0) return 5;

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
    }
    if(m->is_edges)
    {
        is_edges = 1;
        free(m->edges);
        m->edges = NULL;
        m->num_edges = 0;
    }

    for(k=0; k<iters; ++k)
    {
        new_faces = (MESH_FACE)malloc(4*m->num_faces*sizeof(mesh_face));
        v_table = (MESH_STRUCT)malloc(m->num_vertices*sizeof(mesh_struct));
        if(v_table==NULL ||new_faces==NULL) mesh_error(MESH_ERR_MALLOC);
        for(i=0; i<m->num_vertices; ++i)
        {
            v_table[i].num_items = 0;
            v_table[i].items = NULL;
        }
        kv = m->num_vertices-1;
        for(i=0; i<m->num_faces; ++i)
        {
            i_01 = __mesh_find(&v_table[m->faces[i].vertices[0]], m->faces[i].vertices[1]);

            if(i_01<0)
            {
                ++kv;
                c_v_0_indcs = kv;
                v_table[m->faces[i].vertices[0]].items = (INTDATA*)realloc(v_table[m->faces[i].vertices[0]].items, sizeof(INTDATA)*(v_table[m->faces[i].vertices[0]].num_items+2));
                v_table[m->faces[i].vertices[0]].num_items += 2;
                v_table[m->faces[i].vertices[0]].items[v_table[m->faces[i].vertices[0]].num_items-2] = m->faces[i].vertices[1];
                v_table[m->faces[i].vertices[0]].items[v_table[m->faces[i].vertices[0]].num_items-1] = c_v_0_indcs;

                v_table[m->faces[i].vertices[1]].items = (INTDATA*)realloc(v_table[m->faces[i].vertices[1]].items, sizeof(INTDATA)*(v_table[m->faces[i].vertices[1]].num_items+2));
                v_table[m->faces[i].vertices[1]].num_items += 2;
                v_table[m->faces[i].vertices[1]].items[v_table[m->faces[i].vertices[1]].num_items-2] = m->faces[i].vertices[0];
                v_table[m->faces[i].vertices[1]].items[v_table[m->faces[i].vertices[1]].num_items-1] = c_v_0_indcs;

                m->vertices = (MESH_VERTEX)realloc(m->vertices, sizeof(mesh_vertex)*(kv+1));
                m->vertices[kv].x = 0.5*(m->vertices[m->faces[i].vertices[0]].x+m->vertices[m->faces[i].vertices[1]].x);
                m->vertices[kv].y = 0.5*(m->vertices[m->faces[i].vertices[0]].y+m->vertices[m->faces[i].vertices[1]].y);
                m->vertices[kv].z = 0.5*(m->vertices[m->faces[i].vertices[0]].z+m->vertices[m->faces[i].vertices[1]].z);
                if(m->is_vnormals)
                {
                    m->vnormals = (MESH_NORMAL)realloc(m->vnormals, sizeof(mesh_normal)*(kv+1));
                    m->vnormals[kv].x = m->vnormals[m->faces[i].vertices[0]].x+m->vnormals[m->faces[i].vertices[1]].x;
                    m->vnormals[kv].y = m->vnormals[m->faces[i].vertices[0]].y+m->vnormals[m->faces[i].vertices[1]].y;
                    m->vnormals[kv].z = m->vnormals[m->faces[i].vertices[0]].z+m->vnormals[m->faces[i].vertices[1]].z;
                    t = sqrt(m->vnormals[kv].x*m->vnormals[kv].x+m->vnormals[kv].y*m->vnormals[kv].y+m->vnormals[kv].z*m->vnormals[kv].z);
                    m->vnormals[kv].x /= t;
                    m->vnormals[kv].y /= t;
                    m->vnormals[kv].z /= t;
                }
                if(m->is_vcolors)
                {
                    m->vcolors = (MESH_COLOR)realloc(m->vcolors, sizeof(mesh_color)*(kv+1));
                    m->vcolors[kv].r = 0.5*(m->vcolors[m->faces[i].vertices[0]].r+m->vcolors[m->faces[i].vertices[1]].r);
                    m->vcolors[kv].g = 0.5*(m->vcolors[m->faces[i].vertices[0]].g+m->vcolors[m->faces[i].vertices[1]].g);
                    m->vcolors[kv].b = 0.5*(m->vcolors[m->faces[i].vertices[0]].b+m->vcolors[m->faces[i].vertices[1]].b);
                    m->vcolors[kv].a = 0.5*(m->vcolors[m->faces[i].vertices[0]].a+m->vcolors[m->faces[i].vertices[1]].a);
                }
            }
            else
            {
                c_v_0_indcs = v_table[m->faces[i].vertices[0]].items[i_01+1];
            }

            i_12 = __mesh_find(&v_table[m->faces[i].vertices[1]], m->faces[i].vertices[2]);

            if(i_12<0)
            {
                ++kv;
                c_v_1_indcs = kv;
                v_table[m->faces[i].vertices[1]].items = (INTDATA*)realloc(v_table[m->faces[i].vertices[1]].items, sizeof(INTDATA)*(v_table[m->faces[i].vertices[1]].num_items+2));
                v_table[m->faces[i].vertices[1]].num_items += 2;
                v_table[m->faces[i].vertices[1]].items[v_table[m->faces[i].vertices[1]].num_items-2] = m->faces[i].vertices[2];
                v_table[m->faces[i].vertices[1]].items[v_table[m->faces[i].vertices[1]].num_items-1] = c_v_1_indcs;

                v_table[m->faces[i].vertices[2]].items = (INTDATA*)realloc(v_table[m->faces[i].vertices[2]].items, sizeof(INTDATA)*(v_table[m->faces[i].vertices[2]].num_items+2));
                v_table[m->faces[i].vertices[2]].num_items += 2;
                v_table[m->faces[i].vertices[2]].items[v_table[m->faces[i].vertices[2]].num_items-2] = m->faces[i].vertices[1];
                v_table[m->faces[i].vertices[2]].items[v_table[m->faces[i].vertices[2]].num_items-1] = c_v_1_indcs;

                m->vertices = (MESH_VERTEX)realloc(m->vertices, sizeof(mesh_vertex)*(kv+1));
                m->vertices[kv].x = 0.5*(m->vertices[m->faces[i].vertices[1]].x+m->vertices[m->faces[i].vertices[2]].x);
                m->vertices[kv].y = 0.5*(m->vertices[m->faces[i].vertices[1]].y+m->vertices[m->faces[i].vertices[2]].y);
                m->vertices[kv].z = 0.5*(m->vertices[m->faces[i].vertices[1]].z+m->vertices[m->faces[i].vertices[2]].z);
                if(m->is_vnormals)
                {
                    m->vnormals = (MESH_NORMAL)realloc(m->vnormals, sizeof(mesh_normal)*(kv+1));
                    m->vnormals[kv].x = m->vnormals[m->faces[i].vertices[1]].x+m->vnormals[m->faces[i].vertices[2]].x;
                    m->vnormals[kv].y = m->vnormals[m->faces[i].vertices[1]].y+m->vnormals[m->faces[i].vertices[2]].y;
                    m->vnormals[kv].z = m->vnormals[m->faces[i].vertices[1]].z+m->vnormals[m->faces[i].vertices[2]].z;
                    t = sqrt(m->vnormals[kv].x*m->vnormals[kv].x+m->vnormals[kv].y*m->vnormals[kv].y+m->vnormals[kv].z*m->vnormals[kv].z);
                    m->vnormals[kv].x /= t;
                    m->vnormals[kv].y /= t;
                    m->vnormals[kv].z /= t;
                }
                if(m->is_vcolors)
                {
                    m->vcolors = (MESH_COLOR)realloc(m->vcolors, sizeof(mesh_color)*(kv+1));
                    m->vcolors[kv].r = 0.5*(m->vcolors[m->faces[i].vertices[1]].r+m->vcolors[m->faces[i].vertices[2]].r);
                    m->vcolors[kv].g = 0.5*(m->vcolors[m->faces[i].vertices[1]].g+m->vcolors[m->faces[i].vertices[2]].g);
                    m->vcolors[kv].b = 0.5*(m->vcolors[m->faces[i].vertices[1]].b+m->vcolors[m->faces[i].vertices[2]].b);
                    m->vcolors[kv].a = 0.5*(m->vcolors[m->faces[i].vertices[1]].a+m->vcolors[m->faces[i].vertices[2]].a);
                }
            }
            else
            {
                c_v_1_indcs = v_table[m->faces[i].vertices[1]].items[i_12+1];
            }

            i_20 = __mesh_find(&v_table[m->faces[i].vertices[2]], m->faces[i].vertices[0]);

            if(i_20<0)
            {
                ++kv;
                c_v_2_indcs = kv;
                v_table[m->faces[i].vertices[2]].items = (INTDATA*)realloc(v_table[m->faces[i].vertices[2]].items, sizeof(INTDATA)*(v_table[m->faces[i].vertices[2]].num_items+2));
                v_table[m->faces[i].vertices[2]].num_items += 2;
                v_table[m->faces[i].vertices[2]].items[v_table[m->faces[i].vertices[2]].num_items-2] = m->faces[i].vertices[0];
                v_table[m->faces[i].vertices[2]].items[v_table[m->faces[i].vertices[2]].num_items-1] = c_v_2_indcs;

                v_table[m->faces[i].vertices[0]].items = (INTDATA*)realloc(v_table[m->faces[i].vertices[0]].items, sizeof(INTDATA)*(v_table[m->faces[i].vertices[0]].num_items+2));
                v_table[m->faces[i].vertices[0]].num_items += 2;
                v_table[m->faces[i].vertices[0]].items[v_table[m->faces[i].vertices[0]].num_items-2] = m->faces[i].vertices[2];
                v_table[m->faces[i].vertices[0]].items[v_table[m->faces[i].vertices[0]].num_items-1] = c_v_2_indcs;

                m->vertices = (MESH_VERTEX)realloc(m->vertices, sizeof(mesh_vertex)*(kv+1));
                m->vertices[kv].x = 0.5*(m->vertices[m->faces[i].vertices[2]].x+m->vertices[m->faces[i].vertices[0]].x);
                m->vertices[kv].y = 0.5*(m->vertices[m->faces[i].vertices[2]].y+m->vertices[m->faces[i].vertices[0]].y);
                m->vertices[kv].z = 0.5*(m->vertices[m->faces[i].vertices[2]].z+m->vertices[m->faces[i].vertices[0]].z);
                if(m->is_vnormals)
                {
                    m->vnormals = (MESH_NORMAL)realloc(m->vnormals, sizeof(mesh_normal)*(kv+1));
                    m->vnormals[kv].x = m->vnormals[m->faces[i].vertices[2]].x+m->vnormals[m->faces[i].vertices[0]].x;
                    m->vnormals[kv].y = m->vnormals[m->faces[i].vertices[2]].y+m->vnormals[m->faces[i].vertices[0]].y;
                    m->vnormals[kv].z = m->vnormals[m->faces[i].vertices[2]].z+m->vnormals[m->faces[i].vertices[0]].z;
                    t = sqrt(m->vnormals[kv].x*m->vnormals[kv].x+m->vnormals[kv].y*m->vnormals[kv].y+m->vnormals[kv].z*m->vnormals[kv].z);
                    m->vnormals[kv].x /= t;
                    m->vnormals[kv].y /= t;
                    m->vnormals[kv].z /= t;
                }
                if(m->is_vcolors)
                {
                    m->vcolors = (MESH_COLOR)realloc(m->vcolors, sizeof(mesh_color)*(kv+1));
                    m->vcolors[kv].r = 0.5*(m->vcolors[m->faces[i].vertices[2]].r+m->vcolors[m->faces[i].vertices[0]].r);
                    m->vcolors[kv].g = 0.5*(m->vcolors[m->faces[i].vertices[2]].g+m->vcolors[m->faces[i].vertices[0]].g);
                    m->vcolors[kv].b = 0.5*(m->vcolors[m->faces[i].vertices[2]].b+m->vcolors[m->faces[i].vertices[0]].b);
                    m->vcolors[kv].a = 0.5*(m->vcolors[m->faces[i].vertices[2]].a+m->vcolors[m->faces[i].vertices[0]].a);
                }
            }
            else
            {
                c_v_2_indcs = v_table[m->faces[i].vertices[2]].items[i_20+1];
            }

            curr_idx = i*4;
            new_faces[curr_idx].num_vertices = 3;
            new_faces[curr_idx].vertices = (INTDATA *)malloc(3*sizeof(INTDATA));
            new_faces[curr_idx].vertices[0] = m->faces[i].vertices[0];
            new_faces[curr_idx].vertices[1] = c_v_0_indcs;
            new_faces[curr_idx].vertices[2] = c_v_2_indcs;

            ++curr_idx;
            new_faces[curr_idx].num_vertices = 3;
            new_faces[curr_idx].vertices = (INTDATA *)malloc(3*sizeof(INTDATA));
            new_faces[curr_idx].vertices[0] = c_v_0_indcs;
            new_faces[curr_idx].vertices[1] = m->faces[i].vertices[1];
            new_faces[curr_idx].vertices[2] = c_v_1_indcs;

            ++curr_idx;
            new_faces[curr_idx].num_vertices = 3;
            new_faces[curr_idx].vertices = (INTDATA *)malloc(3*sizeof(INTDATA));
            new_faces[curr_idx].vertices[0] = c_v_1_indcs;
            new_faces[curr_idx].vertices[1] = m->faces[i].vertices[2];
            new_faces[curr_idx].vertices[2] = c_v_2_indcs;

            ++curr_idx;
            new_faces[curr_idx].num_vertices = 3;
            new_faces[curr_idx].vertices = (INTDATA *)malloc(3*sizeof(INTDATA));
            new_faces[curr_idx].vertices[0] = c_v_0_indcs;
            new_faces[curr_idx].vertices[1] = c_v_1_indcs;
            new_faces[curr_idx].vertices[2] = c_v_2_indcs;
        }

        for(i=0; i<m->num_vertices; ++i)
        {
            if(v_table[i].items!=NULL) free(v_table[i].items);
        }
        free(v_table);
        if(m->is_vfaces)
        {
            mesh_calc_vertex_adjacency(m);
        }
        free(m->faces);
        m->faces = new_faces;
        m->num_faces = 4*m->num_faces;
        m->num_vertices = kv+1;
    }
    if(is_vfaces==1) mesh_calc_vertex_adjacency(m);
    if(is_ffaces==1) mesh_calc_face_adjacency(m);
    if(is_edges==1) mesh_calc_edges(m);
    return 0;
}

/** \brief Computes area of a triangle
 *
 * \param[in] a First vertex
 * \param[in] b Second vertex
 * \param[in] c Third vertex
 * \return Area
 *
 */

FLOATDATA mesh_calc_triangle_area(MESH_VERTEX a, MESH_VERTEX b, MESH_VERTEX c)
{
    static mesh_vertex p, q, r;
    static FLOATDATA area;
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

