/**
 * @file meshrand.c
 * @author Sk. Mohammadul Haque
 * @version 1.8.0.0
 * @copyright
 * Copyright (c) 2013-2021 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to different mesh random perturbations.
 */

#include "../include/meshlib.h"
#include <stdlib.h>
#include <time.h>


unsigned int MESH_RAND_SEED = 0;
int MESH_SET_RAND_SEED = 0;

void mesh_set_seed(int seed)
{
    if(seed==0)
    {
        MESH_RAND_SEED = (unsigned int)time(NULL)*37519+rand();
        srand((unsigned int)MESH_RAND_SEED);
    }
    else
    {
        MESH_RAND_SEED = seed;
        srand((unsigned int)MESH_RAND_SEED);
    }
    MESH_SET_RAND_SEED = 1;
}

FLOATDATA __mesh_rand(void)
{
    if(!MESH_SET_RAND_SEED) mesh_set_seed(0);
    return ((FLOATDATA)rand())/((FLOATDATA)RAND_MAX+1.0);
}

FLOATDATA __mesh_randn(void)
{
    FLOATDATA tmp0, tmp1;
    if(!MESH_SET_RAND_SEED) mesh_set_seed(0);
    tmp0 = ((FLOATDATA)rand())/((FLOATDATA)RAND_MAX+1.0f);
    tmp1 = ((FLOATDATA)rand())/((FLOATDATA)RAND_MAX+1.0f);
    if(tmp0!=0 ) tmp0 = (FLOATDATA)(sqrt(-2.0f*log(tmp0))*cos(2.0f*3.141592f*tmp1));
    return tmp0;
}

FLOATDATA __mesh_randexp(void)
{
    if(!MESH_SET_RAND_SEED) mesh_set_seed(0);
    return (FLOATDATA)(-log(1.00000000000001f-((FLOATDATA)rand())/((FLOATDATA)RAND_MAX+1.0f)));
}

FLOATDATA __mesh_randfun(FLOATDATA(*fun)(FLOATDATA), FLOATDATA xmin, FLOATDATA xmax)
{
    static FLOATDATA (*Fun)(FLOATDATA) = NULL, YMin, YMax;
    FLOATDATA X, Y;
    int iX;
    if(!MESH_SET_RAND_SEED) mesh_set_seed(0);
    if(fun!=Fun)
    {
        Fun = fun;
        YMin = 0;
        YMax = Fun(xmin);
        for(iX=1; iX<RAND_MAX; ++iX)
        {
            X = xmin+(xmax-xmin)*iX/((FLOATDATA)RAND_MAX+1.0);
            Y = Fun(X);
            YMax = Y>YMax?Y:YMax;
        }
    }
    X = xmin+(xmax-xmin)*((FLOATDATA)rand())/((FLOATDATA)RAND_MAX+1.0);
    Y = YMin+(YMax-YMin)*((FLOATDATA)rand())/((FLOATDATA)RAND_MAX+1.0);;
    return Y <= fun(X) ? X : __mesh_randfun(Fun, xmin, xmax);
}

/** \brief Adds uniform random noise to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \return Error code
*
*/

int mesh_add_noise_uniform(MESH m, FLOATDATA sigma)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    for(INTDATA i=0; i<nv; ++i)
    {
        mv[i].x += sigma*__mesh_rand();
        mv[i].y += sigma*__mesh_rand();
        mv[i].z += sigma*__mesh_rand();
    }
    return 0;
}

/** \brief Adds Gaussian random noise to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \return Error code
*
*/

int mesh_add_noise_gaussian(MESH m, FLOATDATA sigma)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    for(INTDATA i=0; i<nv; ++i)
    {
        mv[i].x += sigma*__mesh_randn();
        mv[i].y += sigma*__mesh_randn();
        mv[i].z += sigma*__mesh_randn();
    }
    return 0;
}

/** \brief Adds exponential random noise to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \return Error code
*
*/

int mesh_add_noise_exp(MESH m, FLOATDATA sigma)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    for(INTDATA i=0; i<nv; ++i)
    {
        mv[i].x += sigma*__mesh_randexp();
        mv[i].y += sigma*__mesh_randexp();
        mv[i].z += sigma*__mesh_randexp();
    }
    return 0;
}

/** \brief Adds random noise given by density function to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \param[in] fun Density function
* \param[in] xmin Lower limit
* \param[in] xmax Upper limit
* \return Error code
*
*/

int mesh_add_noise_func(MESH m, FLOATDATA sigma, FLOATDATA(*fun)(FLOATDATA), FLOATDATA xmin, FLOATDATA xmax)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    for(INTDATA i=0; i<nv; ++i)
    {
        mv[i].x += sigma*__mesh_randfun(fun, xmin, xmax);
        mv[i].y += sigma*__mesh_randfun(fun, xmin, xmax);
        mv[i].z += sigma*__mesh_randfun(fun, xmin, xmax);
    }
    return 0;
}

/** \brief Adds uniform random noise along normals to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \return Error code
*
*/

int mesh_add_noise_uniform_normal(MESH m, FLOATDATA sigma)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    if(!m->is_vnormals)
    {
        int rv = mesh_calc_vertex_normals(m);
        if(rv!=0) return rv;
    }
    {
        const MESH_NORMAL mn = m->vnormals;
        for(INTDATA i=0; i<nv; ++i)
        {
            FLOATDATA r = sigma*__mesh_rand();
            mv[i].x += mn[i].x*r;
            mv[i].y += mn[i].y*r;
            mv[i].z += mn[i].z*r;
        }
    }
    return 0;
}

/** \brief Adds Gaussian random noise along normals to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \return Error code
*
*/

int mesh_add_noise_gaussian_normal(MESH m, FLOATDATA sigma)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    if(!m->is_vnormals)
    {
        int rv = mesh_calc_vertex_normals(m);
        if(rv!=0) return rv;
    }
    {
        const MESH_NORMAL mn = m->vnormals;
        for(INTDATA i=0; i<nv; ++i)
        {
            FLOATDATA r = sigma*__mesh_randn();
            mv[i].x += mn[i].x*r;
            mv[i].y += mn[i].y*r;
            mv[i].z += mn[i].z*r;
        }
    }
    return 0;
}

/** \brief Adds exponential random noise along normals to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \return Error code
*
*/

int mesh_add_noise_exp_normal(MESH m, FLOATDATA sigma)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    if(!m->is_vnormals)
    {
        int rv = mesh_calc_vertex_normals(m);
        if(rv!=0) return rv;
    }
    {
        const MESH_NORMAL mn = m->vnormals;
        for(INTDATA i=0; i<nv; ++i)
        {
            FLOATDATA r = sigma*__mesh_randexp();
            mv[i].x += mn[i].x*r;
            mv[i].y += mn[i].y*r;
            mv[i].z += mn[i].z*r;
        }
    }
    return 0;
}

/** \brief Adds random noise given by density function along normals to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \param[in] fun Density function
* \param[in] xmin Lower limit
* \param[in] xmax Upper limit
* \return Error code
*
*/

int mesh_add_noise_func_normal(MESH m, FLOATDATA sigma, FLOATDATA(*fun)(FLOATDATA), FLOATDATA xmin, FLOATDATA xmax)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    if(!m->is_vnormals)
    {
        int rv = mesh_calc_vertex_normals(m);
        if(rv!=0) return rv;
    }
    {
        const MESH_NORMAL mn = m->vnormals;
        for(INTDATA i=0; i<nv; ++i)
        {
            FLOATDATA r = sigma*__mesh_randfun(fun, xmin, xmax);
            mv[i].x += mn[i].x*r;
            mv[i].y += mn[i].y*r;
            mv[i].z += mn[i].z*r;
        }
    }
    return 0;
}

/** \brief Adds uniform random noise along tangent planes to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \return Error code
*
*/

int mesh_add_noise_uniform_tangent(MESH m, FLOATDATA sigma)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    if(!m->is_vnormals)
    {
        int rv = mesh_calc_vertex_normals(m);
        if(rv!=0) return rv;
    }
    {
        const MESH_NORMAL mn = m->vnormals;
        for(INTDATA i=0; i<nv; ++i)
        {
            FLOATDATA rx = sigma*__mesh_rand();
            FLOATDATA ry = sigma*__mesh_rand();
            FLOATDATA rz = sigma*__mesh_rand();
            MESH_NORMAL cvn = mn+i;

            FLOATDATA axx = (1.0-cvn->x*cvn->x);
            FLOATDATA axy = -(cvn->x*cvn->y);
            FLOATDATA axz = -(cvn->x*cvn->z);
            FLOATDATA ayy = (1.0-cvn->y*cvn->y);
            FLOATDATA ayz = -(cvn->y*cvn->z);
            FLOATDATA azz = (1.0-cvn->z*cvn->z);
            mv[i].x += axx*rx+axy*ry+axz*rz;
            mv[i].y += axy*rx+ayy*ry+ayz*rz;
            mv[i].z += axz*rx+ayz*ry+azz*rz;
        }
    }
    return 0;
}

/** \brief Adds Gaussian random noise along tangent planes to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \return Error code
*
*/

int mesh_add_noise_gaussian_tangent(MESH m, FLOATDATA sigma)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    if(!m->is_vnormals)
    {
        int rv = mesh_calc_vertex_normals(m);
        if(rv!=0) return rv;
    }
    {
        const MESH_NORMAL mn = m->vnormals;
        for(INTDATA i=0; i<nv; ++i)
        {
            FLOATDATA rx = sigma*__mesh_randn();
            FLOATDATA ry = sigma*__mesh_randn();
            FLOATDATA rz = sigma*__mesh_randn();
            MESH_NORMAL cvn = mn+i;

            FLOATDATA axx = (1.0-cvn->x*cvn->x);
            FLOATDATA axy = -(cvn->x*cvn->y);
            FLOATDATA axz = -(cvn->x*cvn->z);
            FLOATDATA ayy = (1.0-cvn->y*cvn->y);
            FLOATDATA ayz = -(cvn->y*cvn->z);
            FLOATDATA azz = (1.0-cvn->z*cvn->z);
            mv[i].x += axx*rx+axy*ry+axz*rz;
            mv[i].y += axy*rx+ayy*ry+ayz*rz;
            mv[i].z += axz*rx+ayz*ry+azz*rz;
        }
    }
    return 0;
}

/** \brief Adds exponential random noise along tangent planes to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \return Error code
*
*/

int mesh_add_noise_exp_tangent(MESH m, FLOATDATA sigma)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    if(!m->is_vnormals)
    {
        int rv = mesh_calc_vertex_normals(m);
        if(rv!=0) return rv;
    }
    {
        const MESH_NORMAL mn = m->vnormals;
        for(INTDATA i=0; i<nv; ++i)
        {
            FLOATDATA rx = sigma*__mesh_randexp();
            FLOATDATA ry = sigma*__mesh_randexp();
            FLOATDATA rz = sigma*__mesh_randexp();
            MESH_NORMAL cvn = mn+i;

            FLOATDATA axx = (1.0-cvn->x*cvn->x);
            FLOATDATA axy = -(cvn->x*cvn->y);
            FLOATDATA axz = -(cvn->x*cvn->z);
            FLOATDATA ayy = (1.0-cvn->y*cvn->y);
            FLOATDATA ayz = -(cvn->y*cvn->z);
            FLOATDATA azz = (1.0-cvn->z*cvn->z);
            mv[i].x += axx*rx+axy*ry+axz*rz;
            mv[i].y += axy*rx+ayy*ry+ayz*rz;
            mv[i].z += axz*rx+ayz*ry+azz*rz;
        }
    }
    return 0;
}

/** \brief Adds random noise given by density function along tangent planes to a mesh
*
* \param[in] m Input mesh
* \param[in] sigma Standard deviation
* \param[in] fun Density function
* \param[in] xmin Lower limit
* \param[in] xmax Upper limit
* \return Error code
*
*/

int mesh_add_noise_func_tangent(MESH m, FLOATDATA sigma, FLOATDATA(*fun)(FLOATDATA), FLOATDATA xmin, FLOATDATA xmax)
{
    const INTDATA nv = m->num_vertices;
    const MESH_VERTEX mv = m->vertices;
    if(!m->is_vnormals)
    {
        int rv = mesh_calc_vertex_normals(m);
        if(rv!=0) return rv;
    }
    {
        const MESH_NORMAL mn = m->vnormals;
        for(INTDATA i=0; i<nv; ++i)
        {
            FLOATDATA rx = sigma*__mesh_randfun(fun, xmin, xmax);
            FLOATDATA ry = sigma*__mesh_randfun(fun, xmin, xmax);
            FLOATDATA rz = sigma*__mesh_randfun(fun, xmin, xmax);
            MESH_NORMAL cvn = mn+i;

            FLOATDATA axx = (1.0-cvn->x*cvn->x);
            FLOATDATA axy = -(cvn->x*cvn->y);
            FLOATDATA axz = -(cvn->x*cvn->z);
            FLOATDATA ayy = (1.0-cvn->y*cvn->y);
            FLOATDATA ayz = -(cvn->y*cvn->z);
            FLOATDATA azz = (1.0-cvn->z*cvn->z);
            mv[i].x += axx*rx+axy*ry+axz*rz;
            mv[i].y += axy*rx+ayy*ry+ayz*rz;
            mv[i].z += axz*rx+ayz*ry+azz*rz;
        }
    }
    return 0;
}






