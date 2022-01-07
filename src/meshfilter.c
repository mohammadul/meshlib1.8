/**
 * @file meshfilter.c
 * @author Sk. Mohammadul Haque
 * @version 1.8.0.0
 * @copyright
 * Copyright (c) 2013-2021 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to different mesh filtering algorithms.
 */

#include "../include/meshlib.h"

/* Fleishman, Shachar, Iddo Drori, and Daniel Cohen-Or.
 * "Bilateral mesh denoising."
 * ACM Transactions on Graphics (TOG). Vol. 22. No. 3. ACM, 2003.
 */

/** \brief Mesh bilateral filter
*
* \param[in] m Input mesh
* \param[in] sigma_c Range standard deviation
* \param[in] sigma_s Spatial standard deviation
* \param[in] niters Number of iterations
* \return Error code
*
*/

int mesh_filter_bilateral(MESH m, FLOATDATA sigma_c, FLOATDATA sigma_s, int niters)
{
    FLOATDATA isigmac, isigmas;
    INTDATA l;
    MESH_VERTEX est_vertices = NULL;
    mesh_calc_vertex_adjacency(m);
    isigmac = 0.5/(sigma_c*sigma_c);
    isigmas = 0.5/(sigma_s*sigma_s);
    est_vertices = (MESH_VERTEX)malloc(m->num_vertices*sizeof(mesh_vertex));
    if(est_vertices==NULL) mesh_error(MESH_ERR_MALLOC);
    for(l=0; l<niters; ++l)
    {
        INTDATA i;
        mesh_calc_vertex_normals(m);
        #pragma omp parallel for
        for(i=0; i<m->num_vertices; ++i)
        {
            INTDATA j, k;
            MESH_FACE cf;
            FLOATDATA t, tx, ty, tz, h, wc, ws;
            FLOATDATA sum = 0, normalizer = 1e-15;
            for(j=0; j<m->vfaces[i].num_faces; ++j)
            {
                cf = m->faces+(m->vfaces[i].faces[j]);
                for(k=0; k<cf->num_vertices; ++k)
                {
                    if(i!=(cf->vertices[k]))
                    {
                        tx = (m->vertices[cf->vertices[k]].x-m->vertices[i].x);
                        ty = (m->vertices[cf->vertices[k]].y-m->vertices[i].y);
                        tz = (m->vertices[cf->vertices[k]].z-m->vertices[i].z);
                        t = tx*tx+ty*ty+tz*tz;
                        h = tx*m->vnormals[i].x + ty*m->vnormals[i].y+tz*m->vnormals[i].z;
                        wc = exp(-t*isigmac);
                        ws = exp(-h*h*isigmas);
                        sum += wc*ws*h;
                        normalizer += wc*ws;
                    }
                }
            }
            sum /= normalizer;
            est_vertices[i].x = m->vertices[i].x + m->vnormals[i].x*sum;
            est_vertices[i].y = m->vertices[i].y + m->vnormals[i].y*sum;
            est_vertices[i].z = m->vertices[i].z + m->vnormals[i].z*sum;
        }
        for(i=0; i<m->num_vertices; ++i)
        {
            m->vertices[i].x = est_vertices[i].x;
            m->vertices[i].y = est_vertices[i].y;
            m->vertices[i].z = est_vertices[i].z;
        }
    }
    free(est_vertices);
    mesh_calc_face_normals(m);
    return 0;
}

/* Field, David A.
 * "Laplacian smoothing and Delaunay triangulations."
 * Communications in applied numerical methods 4.6 (1988): 709-712.
 */

/** \brief Mesh Laplacian filter
 *
 * \param[in] m Input mesh
 * \param[in] r Amount of diffusion
 * \return Error code
 *
 */

int mesh_filter_laplacian(MESH m, FLOATDATA r)
{
    INTDATA i;
    FLOATDATA one_minus_r;
    MESH_VERTEX new_vertices = NULL;
    one_minus_r = 1.0-r;
    if(!m->is_vfaces) mesh_calc_vertex_adjacency(m);
    if((new_vertices = (MESH_VERTEX)malloc(sizeof(mesh_vertex)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    #pragma omp parallel for
    for(i=0; i<m->num_vertices; ++i)
    {
        mesh_vertex s = {0.0, 0.0, 0.0};
        INTDATA t = 0;
        INTDATA j, k;
        for(j=0; j<m->vfaces[i].num_faces; ++j)
        {
            INTDATA mvfij = m->vfaces[i].faces[j];
            for(k=0; k<m->faces[mvfij].num_vertices; ++k)
            {
                INTDATA mvfijvk = m->faces[mvfij].vertices[k];
                if(i==mvfijvk) continue;
                ++t;
                s.x +=m->vertices[mvfijvk].x;
                s.y +=m->vertices[mvfijvk].y;
                s.z +=m->vertices[mvfijvk].z;
            }
        }
        if(t>0)
        {
            s.x /= t;
            s.y /= t;
            s.z /= t;
            new_vertices[i].x = r*s.x + one_minus_r*m->vertices[i].x;
            new_vertices[i].y = r*s.y + one_minus_r*m->vertices[i].y;
            new_vertices[i].z = r*s.z + one_minus_r*m->vertices[i].z;
        }
        else
        {
            new_vertices[i].x = m->vertices[i].x;
            new_vertices[i].y = m->vertices[i].y;
            new_vertices[i].z = m->vertices[i].z;
        }
    }
    free(m->vertices);
    m->vertices = new_vertices;
    return 0;
}

/** \brief Restricted Mesh Laplacian filter
 *
 * \param[in] m Input mesh
 * \param[in] r Amount of diffusion
 * \param[in] ang Minimum angle in degrees to suppress filtering
 * \return Error code
 *
 */

int mesh_filter_laplacian_restricted(MESH m, FLOATDATA r, FLOATDATA ang)
{
    INTDATA i;
    FLOATDATA one_minus_r;
    MESH_VERTEX new_vertices = NULL;
    ang = cos(ang*MESH_PI/180.0);
    one_minus_r = 1.0-r;
    if(!m->is_vfaces) mesh_calc_vertex_adjacency(m);
    if(!m->is_fnormals) mesh_calc_face_normals(m);
    if(!m->is_vnormals) mesh_calc_vertex_normals(m);
    if((new_vertices = (MESH_VERTEX)malloc(sizeof(mesh_vertex)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    #pragma omp parallel for
    for(i=0; i<m->num_vertices; ++i)
    {
        mesh_vertex s = {0.0, 0.0, 0.0};
        INTDATA t = 0;
        INTDATA b, j, k;
        char flg = 0;
        for(j=0; j<m->vfaces[i].num_faces; ++j)
        {
            b = m->vfaces[i].faces[j];
            if((m->fnormals[b].x*m->vnormals[i].x+m->fnormals[b].y*m->vnormals[i].y+m->fnormals[b].z*m->vnormals[i].z)<ang)
            {
                flg = 1;
                break;
            }
        }
        if(flg==0)
        {
            for(j=0; j<m->vfaces[i].num_faces; ++j)
            {
                for(k=0; k<m->faces[m->vfaces[i].faces[j]].num_vertices; ++k)
                {
                    if(i==m->faces[m->vfaces[i].faces[j]].vertices[k]) continue;
                    ++t;
                    s.x +=m->vertices[m->faces[m->vfaces[i].faces[j]].vertices[k]].x;
                    s.y +=m->vertices[m->faces[m->vfaces[i].faces[j]].vertices[k]].y;
                    s.z +=m->vertices[m->faces[m->vfaces[i].faces[j]].vertices[k]].z;
                }
            }
            if(t>0)
            {
                s.x /= t;
                s.y /= t;
                s.z /= t;
                new_vertices[i].x = r*s.x + one_minus_r*m->vertices[i].x;
                new_vertices[i].y = r*s.y + one_minus_r*m->vertices[i].y;
                new_vertices[i].z = r*s.z + one_minus_r*m->vertices[i].z;
            }
            else
            {
                new_vertices[i].x = m->vertices[i].x;
                new_vertices[i].y = m->vertices[i].y;
                new_vertices[i].z = m->vertices[i].z;
            }
        }
        else
        {
            new_vertices[i].x = m->vertices[i].x;
            new_vertices[i].y = m->vertices[i].y;
            new_vertices[i].z = m->vertices[i].z;
        }
    }
    free(m->vertices);
    m->vertices = new_vertices;
    mesh_calc_face_normals(m);
    mesh_calc_vertex_normals(m);
    return 0;
}

/** \brief Mesh Depth Laplacian filter
*
* \param[in] m Input mesh
* \param[in] r Amount of diffusion
* \param[in] vp View-point
* \return Error code
*
*/

int mesh_filter_laplacian_depth(MESH m, FLOATDATA r, MESH_VECTOR3 vp)
{
    INTDATA i;
    MESH_VERTEX new_vertices = NULL;
    if(!m->is_vfaces) mesh_calc_vertex_adjacency(m);
    if((new_vertices = (MESH_VERTEX)malloc(sizeof(mesh_vertex)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    #pragma omp parallel for
    for(i=0; i<m->num_vertices; ++i)
    {
        mesh_vertex s = {0.0, 0.0, 0.0};
        INTDATA t = 0;
        INTDATA j, k;
        for(j=0; j<m->vfaces[i].num_faces; ++j)
        {
            for(k=0; k<m->faces[m->vfaces[i].faces[j]].num_vertices; ++k)
            {
                if(i==m->faces[m->vfaces[i].faces[j]].vertices[k]) continue;
                ++t;
                s.x +=m->vertices[m->faces[m->vfaces[i].faces[j]].vertices[k]].x;
                s.y +=m->vertices[m->faces[m->vfaces[i].faces[j]].vertices[k]].y;
                s.z +=m->vertices[m->faces[m->vfaces[i].faces[j]].vertices[k]].z;
            }
        }
        if(t>0)
        {
            s.x /= t;
            s.y /= t;
            s.z /= t;
            mesh_vector3 delv = {s.x-m->vertices[i].x,s.y-m->vertices[i].y,s.z-m->vertices[i].z};
            mesh_vector3 dr = {vp->x-m->vertices[i].x,vp->y-m->vertices[i].y,vp->z-m->vertices[i].z};
            FLOATDATA mgn = sqrt(dr.x*dr.x+dr.y*dr.y+dr.z*dr.z);
            if(mgn>0.0)
            {
                dr.x /= mgn;
                dr.y /= mgn;
                dr.z /= mgn;
            }
            FLOATDATA sc = dr.x*delv.x+dr.y*delv.y+dr.z*delv.z;
            new_vertices[i].x = m->vertices[i].x+sc*dr.x;
            new_vertices[i].y = m->vertices[i].y+sc*dr.y;
            new_vertices[i].z = m->vertices[i].z+sc*dr.z;
        }
        else
        {
            new_vertices[i].x = m->vertices[i].x;
            new_vertices[i].y = m->vertices[i].y;
            new_vertices[i].z = m->vertices[i].z;
        }
    }
    free(m->vertices);
    m->vertices = new_vertices;
    return 0;
}

/** \brief Mesh \f$ \lambda-\mu \f$ Taubin filter
 *
 * \param[in] m Input mesh
 * \param[in] lambd \f$\lambda\f$ value
 * \param[in] mu \f$\mu\f$ value
 * \param[in] niters Number of iterations
 * \return Error code
 *
 */

int mesh_filter_taubin(MESH m, FLOATDATA lambd, FLOATDATA mu, int niters)
{
    int i;
    if(m==NULL) return -1;
    if(m->is_faces==0) return -2;
    for(i=0; i<niters; ++i)
    {
        mesh_laplacian_filter(m, lambd);
        mesh_laplacian_filter(m, mu);
    }
    return 0;
}

/** \brief Mesh bilateral vertex color filter
*
* \param[in] m Input mesh
* \param[in] sigma_k Color standard deviation
* \param[in] sigma_c Range standard deviation
* \param[in] sigma_s Spatial standard deviation
* \param[in] niters Number of iterations
* \return Error code
*
*/

int mesh_filter_vertex_color_bilateral(MESH m, FLOATDATA sigma_k, FLOATDATA sigma_c, FLOATDATA sigma_s, int niters)
{
    FLOATDATA isigmak, isigmac, isigmas;
    INTDATA l;
    MESH_COLOR est_vcolors = NULL;
    if(m->is_vcolors==0) return -1;
    mesh_calc_vertex_adjacency(m);
    isigmak = 0.5/(sigma_k*sigma_k);
    isigmac = 0.5/(sigma_c*sigma_c);
    isigmas = 0.5/(sigma_s*sigma_s);
    est_vcolors = (MESH_COLOR)malloc(m->num_vertices*sizeof(mesh_color));
    if(est_vcolors==NULL) mesh_error(MESH_ERR_MALLOC);
    for(l=0; l<niters; ++l)
    {
        INTDATA i;
        MESH_COLOR tmp;
        mesh_calc_vertex_normals(m);
        #pragma omp parallel for
        for(i=0; i<m->num_vertices; ++i)
        {
            INTDATA j, k;
            MESH_FACE cf;
            FLOATDATA normalizer = 1.0;
            MESH_COLOR cvcolor = est_vcolors+i, old_cvcolor = m->vcolors+i;
            cvcolor->r = old_cvcolor->r;
            cvcolor->g = old_cvcolor->g;
            cvcolor->b = old_cvcolor->b;
            cvcolor->a = old_cvcolor->a;

            for(j=0; j<m->vfaces[i].num_faces; ++j)
            {
                cf = m->faces+(m->vfaces[i].faces[j]);
                for(k=0; k<cf->num_vertices; ++k)
                {
                    INTDATA kv = cf->vertices[k];
                    if(i!=kv)
                    {
                        FLOATDATA t, tx, ty, tz, kc, cr, cg, cb, h, wkc, wc, ws, w;
                        MESH_COLOR ovcolor = m->vcolors+k;
                        tx = (m->vertices[kv].x-m->vertices[i].x);
                        ty = (m->vertices[kv].y-m->vertices[i].y);
                        tz = (m->vertices[kv].z-m->vertices[i].z);
                        cr = (ovcolor->r-old_cvcolor->r);
                        cg = (ovcolor->g-old_cvcolor->g);
                        cb = (ovcolor->b-old_cvcolor->b);
                        t = tx*tx+ty*ty+tz*tz;
                        kc = cr*cr+cg*cg+cb*cb;
                        h = tx*m->vnormals[i].x + ty*m->vnormals[i].y+tz*m->vnormals[i].z;
                        wkc = exp(-kc*isigmak);
                        wc = exp(-t*isigmac);
                        ws = exp(-h*h*isigmas);
                        w = wkc*wc*ws;
                        cvcolor->r += w*ovcolor->r;
                        cvcolor->g += w*ovcolor->g;
                        cvcolor->b += w*ovcolor->b;
                        cvcolor->a += w*ovcolor->a;
                        normalizer += w;
                    }
                }
            }
            cvcolor->r /= normalizer;
            cvcolor->g /= normalizer;
            cvcolor->b /= normalizer;
            cvcolor->a /= normalizer;
        }
        tmp = m->vcolors;
        m->vcolors = est_vcolors;
        est_vcolors = tmp;
    }
    free(est_vcolors);
    return 0;
}

/** \brief Mesh Laplacian vertex color filter
 *
 * \param[in] m Input mesh
 * \param[in] r Amount of diffusion
 * \return Error code
 *
 */

int mesh_filter_vertex_color_laplacian(MESH m, FLOATDATA r)
{
    INTDATA i;
    FLOATDATA one_minus_r;
    MESH_COLOR new_vcolors = NULL;
    one_minus_r = 1.0-r;
    if(!m->is_vfaces) mesh_calc_vertex_adjacency(m);
    if((new_vcolors = (MESH_COLOR)malloc(m->num_vertices*sizeof(mesh_color)))==NULL) mesh_error(MESH_ERR_MALLOC);
    #pragma omp parallel for
    for(i=0; i<m->num_vertices; ++i)
    {
        mesh_color s = {0.0, 0.0, 0.0, 0.0};
        INTDATA t = 0;
        INTDATA j, k;
        for(j=0; j<m->vfaces[i].num_faces; ++j)
        {
            INTDATA mvfij = m->vfaces[i].faces[j];
            for(k=0; k<m->faces[mvfij].num_vertices; ++k)
            {
                INTDATA mvfijvk = m->faces[mvfij].vertices[k];
                if(i==mvfijvk) continue;
                ++t;
                s.r +=m->vcolors[mvfijvk].r;
                s.g +=m->vcolors[mvfijvk].g;
                s.b +=m->vcolors[mvfijvk].b;
                s.a +=m->vcolors[mvfijvk].a;
            }
        }
        if(t>0)
        {
            s.r /= t;
            s.g /= t;
            s.b /= t;
            s.a /= t;
            new_vcolors[i].r = r*s.r + one_minus_r*m->vcolors[i].r;
            new_vcolors[i].g = r*s.g + one_minus_r*m->vcolors[i].g;
            new_vcolors[i].b = r*s.b + one_minus_r*m->vcolors[i].b;
            new_vcolors[i].a = r*s.a + one_minus_r*m->vcolors[i].a;
        }
        else
        {
            new_vcolors[i].r = m->vcolors[i].r;
            new_vcolors[i].g = m->vcolors[i].g;
            new_vcolors[i].b = m->vcolors[i].b;
            new_vcolors[i].a = m->vcolors[i].a;
        }
    }
    free(m->vcolors);
    m->vcolors = new_vcolors;
    return 0;
}

/** \brief Mesh minimum intensity vertex color filter
 *
 * \param[in] m Input mesh
 * \param[in] niters Number of iterations
 * \return Error code
 *
 */

int mesh_filter_vertex_color_min(MESH m, INTDATA niters)
{
    INTDATA l, li;
    FLOATDATA *vals = NULL, *newvals = NULL;
    MESH_COLOR new_vcolors = NULL;
    if(!m->is_vfaces) mesh_calc_vertex_adjacency(m);
    if((new_vcolors = (MESH_COLOR)malloc(m->num_vertices*sizeof(mesh_color)))==NULL) mesh_error(MESH_ERR_MALLOC);
    if((vals = (FLOATDATA*)malloc(m->num_vertices*sizeof(FLOATDATA)))==NULL) mesh_error(MESH_ERR_MALLOC);
    if((newvals = (FLOATDATA*)malloc(m->num_vertices*sizeof(FLOATDATA)))==NULL) mesh_error(MESH_ERR_MALLOC);
    #pragma omp parallel for
    for(li=0; li<m->num_vertices; ++li)
    {
        MESH_COLOR cvc = m->vcolors+li;
        vals[li] = cvc->r*cvc->r+cvc->g*cvc->g+cvc->b*cvc->b;
    }
    for(l=0; l<niters; ++l)
    {
        INTDATA i;
        MESH_COLOR tmp;
        FLOATDATA* tmpvals;
        #pragma omp parallel for
        for(i=0; i<m->num_vertices; ++i)
        {
            INTDATA smin = i;
            INTDATA j, k;
            for(j=0; j<m->vfaces[i].num_faces; ++j)
            {
                INTDATA mvfij = m->vfaces[i].faces[j];
                for(k=0; k<m->faces[mvfij].num_vertices; ++k)
                {
                    INTDATA mvfijvk = m->faces[mvfij].vertices[k];
                    if(i==mvfijvk) continue;
                    if(vals[mvfijvk]<vals[smin])
                    {
                        smin = mvfijvk;
                    }
                }
            }
            new_vcolors[i].r = m->vcolors[smin].r;
            new_vcolors[i].g = m->vcolors[smin].g;
            new_vcolors[i].b = m->vcolors[smin].b;
            new_vcolors[i].a = m->vcolors[smin].a;
            newvals[i] = vals[smin];
        }
        tmp = m->vcolors;
        m->vcolors = new_vcolors;
        new_vcolors = tmp;

        tmpvals = vals;
        vals = newvals;
        newvals = tmpvals;

    }
    free(vals);
    free(newvals);
    free(new_vcolors);
    return 0;
}

/** \brief Mesh maximum intensity vertex color filter
 *
 * \param[in] m Input mesh
 * \param[in] niters Number of iterations
 * \return Error code
 *
 */

int mesh_filter_vertex_color_max(MESH m, INTDATA niters)
{
    INTDATA l, li;
    FLOATDATA *vals = NULL, *newvals = NULL;
    MESH_COLOR new_vcolors = NULL;
    if(!m->is_vfaces) mesh_calc_vertex_adjacency(m);
    if((new_vcolors = (MESH_COLOR)malloc(m->num_vertices*sizeof(mesh_color)))==NULL) mesh_error(MESH_ERR_MALLOC);
    if((vals = (FLOATDATA*)malloc(m->num_vertices*sizeof(FLOATDATA)))==NULL) mesh_error(MESH_ERR_MALLOC);
    if((newvals = (FLOATDATA*)malloc(m->num_vertices*sizeof(FLOATDATA)))==NULL) mesh_error(MESH_ERR_MALLOC);
    #pragma omp parallel for
    for(li=0; li<m->num_vertices; ++li)
    {
        MESH_COLOR cvc = m->vcolors+li;
        vals[li] = cvc->r*cvc->r+cvc->g*cvc->g+cvc->b*cvc->b;
    }
    for(l=0; l<niters; ++l)
    {
        INTDATA i;
        MESH_COLOR tmp;
        FLOATDATA* tmpvals;
        #pragma omp parallel for
        for(i=0; i<m->num_vertices; ++i)
        {
            INTDATA smax = i;
            INTDATA j, k;
            for(j=0; j<m->vfaces[i].num_faces; ++j)
            {
                INTDATA mvfij = m->vfaces[i].faces[j];
                for(k=0; k<m->faces[mvfij].num_vertices; ++k)
                {
                    INTDATA mvfijvk = m->faces[mvfij].vertices[k];
                    if(i==mvfijvk) continue;
                    if(vals[mvfijvk]>vals[smax])
                    {
                        smax = mvfijvk;
                    }
                }
            }
            new_vcolors[i].r = m->vcolors[smax].r;
            new_vcolors[i].g = m->vcolors[smax].g;
            new_vcolors[i].b = m->vcolors[smax].b;
            new_vcolors[i].a = m->vcolors[smax].a;
            newvals[i] = vals[smax];
        }
        tmp = m->vcolors;
        m->vcolors = new_vcolors;
        new_vcolors = tmp;

        tmpvals = vals;
        vals = newvals;
        newvals = tmpvals;
    }
    free(vals);
    free(newvals);
    free(new_vcolors);
    return 0;
}


