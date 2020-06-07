/**
 * @file meshdraw.c
 * @author Sk. Mohammadul Haque
 * @version 1.4.2.0
 * @copyright
 * Copyright (c) 2013, 2014, 2015, 2016 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to mesh drawing in OpenGL.
 */

#include "../include/meshlib.h"
#if defined(_WIN32) || defined(WIN32)
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>

/** \brief Draws a given mesh in OpenGL context in flat shading
 *
 * \param[in] m Input mesh
 * \return NULL
 *
 */

void mesh_draw_mesh(MESH m)
{
    INTDATA i, j;
    GLfloat currcolor[4];
    if(!m->is_fnormals) mesh_calc_face_normals(m);
    if(m->is_trimesh)
    {
        if(m->is_fcolors)
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            glBegin(GL_TRIANGLES);
            for(i=0; i<m->num_faces; ++i)
            {
                currcolor[0] = m->fcolors[i].r;
                currcolor[1] = m->fcolors[i].g;
                currcolor[2] = m->fcolors[i].b;
                currcolor[3] = m->fcolors[i].a;
                glColor3fv(currcolor);
                glNormal3f(m->fnormals[i].x, m->fnormals[i].y, m->fnormals[i].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[2]].x, m->vertices[m->faces[i].vertices[2]].y, m->vertices[m->faces[i].vertices[2]].z);
            }
            glEnd();
            glDisable(GL_COLOR_MATERIAL);
        }
        else if(m->is_vcolors)
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            glBegin(GL_TRIANGLES);
            for(i=0; i<m->num_faces; ++i)
            {
                currcolor[0] = 0.33333*(m->vcolors[m->faces[i].vertices[0]].r+m->vcolors[m->faces[i].vertices[1]].r+m->vcolors[m->faces[i].vertices[2]].r);
                currcolor[1] = 0.33333*(m->vcolors[m->faces[i].vertices[0]].g+m->vcolors[m->faces[i].vertices[1]].g+m->vcolors[m->faces[i].vertices[2]].g);
                currcolor[2] = 0.33333*(m->vcolors[m->faces[i].vertices[0]].b+m->vcolors[m->faces[i].vertices[1]].b+m->vcolors[m->faces[i].vertices[2]].b);
                currcolor[3] = 0.33333*(m->vcolors[m->faces[i].vertices[0]].a+m->vcolors[m->faces[i].vertices[1]].a+m->vcolors[m->faces[i].vertices[2]].a);
                glColor3fv(currcolor);
                glNormal3f(m->fnormals[i].x, m->fnormals[i].y, m->fnormals[i].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[2]].x, m->vertices[m->faces[i].vertices[2]].y, m->vertices[m->faces[i].vertices[2]].z);
            }
            glEnd();
            glDisable(GL_COLOR_MATERIAL);
        }
        else
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            currcolor[0] = 1.0;
            currcolor[1] = 1.0;
            currcolor[2] = 1.0;
            currcolor[3] = 1.0;
            glColor3fv(currcolor);
            glBegin(GL_TRIANGLES);
            for(i=0; i<m->num_faces; ++i)
            {
                glNormal3f(m->fnormals[i].x, m->fnormals[i].y, m->fnormals[i].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[2]].x, m->vertices[m->faces[i].vertices[2]].y, m->vertices[m->faces[i].vertices[2]].z);
            }
            glEnd();
            glDisable(GL_COLOR_MATERIAL);
        }
    }
    else if(m->num_faces>0)
    {
        if(m->is_fcolors)
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            for(i=0; i<m->num_faces; ++i)
            {
                currcolor[0] = m->fcolors[i].r;
                currcolor[1] = m->fcolors[i].g;
                currcolor[2] = m->fcolors[i].b;
                currcolor[3] = m->fcolors[i].a;
                glColor3fv(currcolor);
                glBegin(GL_TRIANGLE_FAN);
                glNormal3f(m->fnormals[i].x, m->fnormals[i].y, m->fnormals[i].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                for(j=2; j<m->faces[i].num_vertices; ++j)
                {
                    glVertex3f(m->vertices[m->faces[i].vertices[j]].x, m->vertices[m->faces[i].vertices[j]].y, m->vertices[m->faces[i].vertices[j]].z);
                }
                glEnd();
            }
            glDisable(GL_COLOR_MATERIAL);
        }
        if(m->is_vcolors)
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            for(i=0; i<m->num_faces; ++i)
            {
                currcolor[0] = 0.0;
                currcolor[1] = 0.0;
                currcolor[2] = 0.0;
                currcolor[3] = 0.0;
                for(j=0; j<m->faces[i].num_vertices; ++j)
                {
                    INTDATA idx = m->faces[i].vertices[j];
                    currcolor[0] += m->vcolors[idx].r;
                    currcolor[1] += m->vcolors[idx].g;
                    currcolor[2] += m->vcolors[idx].b;
                    currcolor[3] += m->vcolors[idx].a;
                }
                currcolor[0] /= m->faces[i].num_vertices;
                currcolor[1] /= m->faces[i].num_vertices;
                currcolor[2] /= m->faces[i].num_vertices;
                currcolor[3] /= m->faces[i].num_vertices;
                glColor3fv(currcolor);
                glBegin(GL_TRIANGLE_FAN);
                glNormal3f(m->fnormals[i].x, m->fnormals[i].y, m->fnormals[i].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                for(j=2; j<m->faces[i].num_vertices; ++j)
                {
                    glVertex3f(m->vertices[m->faces[i].vertices[j]].x, m->vertices[m->faces[i].vertices[j]].y, m->vertices[m->faces[i].vertices[j]].z);
                }
                glEnd();
            }
            glDisable(GL_COLOR_MATERIAL);
        }
        else
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            currcolor[0] = 1.0;
            currcolor[1] = 1.0;
            currcolor[2] = 1.0;
            currcolor[3] = 1.0;
            glColor3fv(currcolor);
            for(i=0; i<m->num_faces; ++i)
            {
                glBegin(GL_TRIANGLE_FAN);
                glNormal3f(m->fnormals[i].x, m->fnormals[i].y, m->fnormals[i].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                for(j=2; j<m->faces[i].num_vertices; ++j)
                {
                    glVertex3f(m->vertices[m->faces[i].vertices[j]].x, m->vertices[m->faces[i].vertices[j]].y, m->vertices[m->faces[i].vertices[j]].z);
                }
                glEnd();
            }
            glDisable(GL_COLOR_MATERIAL);
        }
    }
}


/** \brief Draws a given mesh in OpenGL context in smoothing shading
 *
 * \param[in] m Input mesh
 * \return NULL
 *
 */

void mesh_draw_mesh_smooth(MESH m)
{
    INTDATA i, j;
    GLfloat currcolor[4];
    if(!m->is_vnormals) mesh_calc_vertex_normals(m);
    if(m->is_trimesh)
    {
        if(m->is_fcolors)
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            glBegin(GL_TRIANGLES);
            for(i=0; i<m->num_faces; ++i)
            {
                currcolor[0] = m->fcolors[i].r;
                currcolor[1] = m->fcolors[i].g;
                currcolor[2] = m->fcolors[i].b;
                currcolor[3] = m->fcolors[i].a;
                glColor3fv(currcolor);
                glNormal3f(m->vnormals[m->faces[i].vertices[0]].x, m->vnormals[m->faces[i].vertices[0]].y, m->vnormals[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glNormal3f(m->vnormals[m->faces[i].vertices[1]].x, m->vnormals[m->faces[i].vertices[1]].y, m->vnormals[m->faces[i].vertices[1]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                glNormal3f(m->vnormals[m->faces[i].vertices[2]].x, m->vnormals[m->faces[i].vertices[2]].y, m->vnormals[m->faces[i].vertices[2]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[2]].x, m->vertices[m->faces[i].vertices[2]].y, m->vertices[m->faces[i].vertices[2]].z);
            }
            glEnd();
            glDisable(GL_COLOR_MATERIAL);
        }
        else if(m->is_vcolors)
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            glBegin(GL_TRIANGLES);
            for(i=0; i<m->num_faces; ++i)
            {
                currcolor[0] = 0.33333*(m->vcolors[m->faces[i].vertices[0]].r+m->vcolors[m->faces[i].vertices[1]].r+m->vcolors[m->faces[i].vertices[2]].r);
                currcolor[1] = 0.33333*(m->vcolors[m->faces[i].vertices[0]].g+m->vcolors[m->faces[i].vertices[1]].g+m->vcolors[m->faces[i].vertices[2]].g);
                currcolor[2] = 0.33333*(m->vcolors[m->faces[i].vertices[0]].b+m->vcolors[m->faces[i].vertices[1]].b+m->vcolors[m->faces[i].vertices[2]].b);
                currcolor[3] = 0.33333*(m->vcolors[m->faces[i].vertices[0]].a+m->vcolors[m->faces[i].vertices[1]].a+m->vcolors[m->faces[i].vertices[2]].a);
                glColor3fv(currcolor);
                glNormal3f(m->vnormals[m->faces[i].vertices[0]].x, m->vnormals[m->faces[i].vertices[0]].y, m->vnormals[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glNormal3f(m->vnormals[m->faces[i].vertices[1]].x, m->vnormals[m->faces[i].vertices[1]].y, m->vnormals[m->faces[i].vertices[1]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                glNormal3f(m->vnormals[m->faces[i].vertices[2]].x, m->vnormals[m->faces[i].vertices[2]].y, m->vnormals[m->faces[i].vertices[2]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[2]].x, m->vertices[m->faces[i].vertices[2]].y, m->vertices[m->faces[i].vertices[2]].z);
            }
            glEnd();
            glDisable(GL_COLOR_MATERIAL);
        }
        else
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            currcolor[0] = 1.0;
            currcolor[1] = 1.0;
            currcolor[2] = 1.0;
            currcolor[3] = 1.0;
            glColor3fv(currcolor);
            glBegin(GL_TRIANGLES);
            for(i=0; i<m->num_faces; ++i)
            {
                glNormal3f(m->vnormals[m->faces[i].vertices[0]].x, m->vnormals[m->faces[i].vertices[0]].y, m->vnormals[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glNormal3f(m->vnormals[m->faces[i].vertices[1]].x, m->vnormals[m->faces[i].vertices[1]].y, m->vnormals[m->faces[i].vertices[1]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                glNormal3f(m->vnormals[m->faces[i].vertices[2]].x, m->vnormals[m->faces[i].vertices[2]].y, m->vnormals[m->faces[i].vertices[2]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[2]].x, m->vertices[m->faces[i].vertices[2]].y, m->vertices[m->faces[i].vertices[2]].z);
            }
            glEnd();
            glDisable(GL_COLOR_MATERIAL);
        }
    }
    else if(m->num_faces>0)
    {
        if(m->is_fcolors)
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            for(i=0; i<m->num_faces; ++i)
            {
                currcolor[0] = m->fcolors[i].r;
                currcolor[1] = m->fcolors[i].g;
                currcolor[2] = m->fcolors[i].b;
                currcolor[3] = m->fcolors[i].a;
                glColor3fv(currcolor);
                glBegin(GL_TRIANGLE_FAN);
                glNormal3f(m->vnormals[m->faces[i].vertices[0]].x, m->vnormals[m->faces[i].vertices[0]].y, m->vnormals[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glNormal3f(m->vnormals[m->faces[i].vertices[1]].x, m->vnormals[m->faces[i].vertices[1]].y, m->vnormals[m->faces[i].vertices[1]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                for(j=2; j<m->faces[i].num_vertices; ++j)
                {
                    glNormal3f(m->vnormals[m->faces[i].vertices[j]].x, m->vnormals[m->faces[i].vertices[j]].y, m->vnormals[m->faces[i].vertices[j]].z);
                    glVertex3f(m->vertices[m->faces[i].vertices[j]].x, m->vertices[m->faces[i].vertices[j]].y, m->vertices[m->faces[i].vertices[j]].z);
                }
                glEnd();
            }
            glDisable(GL_COLOR_MATERIAL);
        }
        if(m->is_vcolors)
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            for(i=0; i<m->num_faces; ++i)
            {
                glBegin(GL_TRIANGLE_FAN);
                for(j=0; j<m->faces[i].num_vertices; ++j)
                {
                    INTDATA idx = m->faces[i].vertices[j];
                    currcolor[0] = m->vcolors[idx].r;
                    currcolor[1] = m->vcolors[idx].g;
                    currcolor[2] = m->vcolors[idx].b;
                    currcolor[3] = m->vcolors[idx].a;
                    glColor3fv(currcolor);
                    glNormal3f(m->vnormals[idx].x, m->vnormals[idx].y, m->vnormals[idx].z);
                    glVertex3f(m->vertices[idx].x, m->vertices[idx].y, m->vertices[idx].z);
                }
                glEnd();

            }
            glDisable(GL_COLOR_MATERIAL);
        }
        else
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            currcolor[0] = 1.0;
            currcolor[1] = 1.0;
            currcolor[2] = 1.0;
            currcolor[3] = 1.0;
            glColor3fv(currcolor);
            for(i=0; i<m->num_faces; ++i)
            {
                glBegin(GL_TRIANGLE_FAN);
                glNormal3f(m->vnormals[m->faces[i].vertices[0]].x, m->vnormals[m->faces[i].vertices[0]].y, m->vnormals[m->faces[i].vertices[0]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[0]].x, m->vertices[m->faces[i].vertices[0]].y, m->vertices[m->faces[i].vertices[0]].z);
                glNormal3f(m->vnormals[m->faces[i].vertices[1]].x, m->vnormals[m->faces[i].vertices[1]].y, m->vnormals[m->faces[i].vertices[1]].z);
                glVertex3f(m->vertices[m->faces[i].vertices[1]].x, m->vertices[m->faces[i].vertices[1]].y, m->vertices[m->faces[i].vertices[1]].z);
                for(j=2; j<m->faces[i].num_vertices; ++j)
                {
                    glNormal3f(m->vnormals[m->faces[i].vertices[j]].x, m->vnormals[m->faces[i].vertices[j]].y, m->vnormals[m->faces[i].vertices[j]].z);
                    glVertex3f(m->vertices[m->faces[i].vertices[j]].x, m->vertices[m->faces[i].vertices[j]].y, m->vertices[m->faces[i].vertices[j]].z);
                }
                glEnd();
            }
            glDisable(GL_COLOR_MATERIAL);
        }
    }
}


/** \brief Draws a given mesh in OpenGL context as pointcloud
 *
 * \param[in] m Input mesh
 * \return NULL
 *
 */

void mesh_draw_point_cloud(MESH m)
{
    INTDATA i;
    GLfloat currcolor[4];
    if(!m->is_vnormals) mesh_calc_vertex_normals(m);
    if(m->is_vcolors)
    {
        if(m->is_vnormals)
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            glBegin(GL_POINTS);
            for(i=0; i<m->num_vertices; ++i)
            {
                currcolor[0] = m->vcolors[i].r;
                currcolor[1] = m->vcolors[i].g;
                currcolor[2] = m->vcolors[i].b;
                currcolor[3] = m->vcolors[i].a;

                glColor3fv(currcolor);
                glNormal3f(m->vnormals[i].x, m->vnormals[i].y, m->vnormals[i].z);
                glVertex3f(m->vertices[i].x, m->vertices[i].y, m->vertices[i].z);
            }
            glEnd();
            glDisable(GL_COLOR_MATERIAL);
        }
        else
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            glBegin(GL_POINTS);
            for(i=0; i<m->num_vertices; ++i)
            {
                currcolor[0] = m->vcolors[i].r;
                currcolor[1] = m->vcolors[i].g;
                currcolor[2] = m->vcolors[i].b;
                currcolor[3] = m->vcolors[i].a;

                glColor3fv(currcolor);
                glVertex3f(m->vertices[i].x, m->vertices[i].y, m->vertices[i].z);
            }
            glEnd();
            glDisable(GL_COLOR_MATERIAL);
        }
    }
    else
    {
        if(m->is_vnormals)
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            glBegin(GL_POINTS);
            for(i=0; i<m->num_vertices; ++i)
            {
                currcolor[0] = 1.0;
                currcolor[1] = 1.0;
                currcolor[2] = 1.0;
                currcolor[3] = 1.0;
                glColor3fv(currcolor);
                glNormal3f(m->vnormals[i].x, m->vnormals[i].y, m->vnormals[i].z);
                glVertex3f(m->vertices[i].x, m->vertices[i].y, m->vertices[i].z);
            }
            glEnd();
            glDisable(GL_COLOR_MATERIAL);
        }
        else
        {
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
            glBegin(GL_POINTS);
            for(i=0; i<m->num_vertices; ++i)
            {
                currcolor[0] = 1.0;
                currcolor[1] = 1.0;
                currcolor[2] = 1.0;
                currcolor[3] = 1.0;
                glColor3fv(currcolor);
                glVertex3f(m->vertices[i].x, m->vertices[i].y, m->vertices[i].z);
            }
            glEnd();
            glDisable(GL_COLOR_MATERIAL);
        }
    }
}
