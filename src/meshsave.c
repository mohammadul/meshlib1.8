/**
 * @file meshsave.c
 * @author Sk. Mohammadul Haque
 * @version 1.8.0.0
 * @copyright
 * Copyright (c) 2013-2021 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to saving different mesh file types.
 */
#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#ifndef _MSC_VER
#include <inttypes.h>
#else
#define PRId64 "I64d"
#endif
#include "../include/meshlib.h"

/** \brief Saves a mesh to an OFF/PLY/ASC/XYZ/BIN/OBJ file
 *
 * \param[in] m Input mesh
 * \param[in] fname Output filename
 * \return Error code
 *
 */

int mesh_save_file(MESH m, const char* fname)
{
    const char *ext = strrchr(fname, '.');
    if(strcmpi(ext,".off")==0) return mesh_write_off(m, fname);
    else if(strcmpi(ext,".ply")==0) return mesh_save_ply(m, fname);
    else if(strcmpi(ext,".asc")==0) return mesh_save_xyz(m, fname);
    else if(strcmpi(ext,".xyz")==0) return mesh_save_xyz(m, fname);
    else if(strcmpi(ext,".bin")==0) return mesh_save_bin(m, fname);
    else if(strcmpi(ext,".obj")==0) return mesh_save_obj(m, fname);
    return -1;
}

/** \brief Saves a mesh to an OFF file
 *
 * \param[in] m Input mesh
 * \param[in] fname Output filename
 * \return Error code
 *
 */

int mesh_save_off(MESH m, const char* fname)
{
    INTDATA i, j;
    if(m->is_vertices)
    {
        FILEPOINTER fp = NULL;
        if((fp = fopen(fname,"wb"))==NULL) mesh_error(MESH_ERR_FNOTOPEN);
        if(m->is_vcolors && !m->is_vnormals)
        {
            fprintf(fp, "COFF\n");
#if MESH_INTDATA_TYPE==0
            fprintf(fp, "%d %d 0\n", m->num_vertices, m->num_faces);
#else
            fprintf(fp, "%"PRId64" %"PRId64" 0\n", m->num_vertices, m->num_faces);
#endif
            for(i=0; i<m->num_vertices; ++i)
            {
                fprintf(fp, "%f %f %f %f %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z, m->vcolors[i].r,  m->vcolors[i].g, m->vcolors[i].b, m->vcolors[i].a);
            }
        }
        else if(m->is_vnormals && !m->is_vcolors)
        {
            fprintf(fp, "NOFF\n");
#if MESH_INTDATA_TYPE==0
            fprintf(fp, "%d %d 0\n", m->num_vertices, m->num_faces);
#else
            fprintf(fp, "%"PRId64" %"PRId64" 0\n", m->num_vertices, m->num_faces);
#endif
            for(i=0; i<m->num_vertices; ++i)
            {
                fprintf(fp, "%f %f %f %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z, m->vnormals[i].x,  m->vnormals[i].y, m->vnormals[i].z);
            }
        }
        else if(m->is_vcolors && m->is_vnormals)
        {
            fprintf(fp, "NCOFF\n");
#if MESH_INTDATA_TYPE==0
            fprintf(fp, "%d %d 0\n", m->num_vertices, m->num_faces);
#else
            fprintf(fp, "%"PRId64" %"PRId64" 0\n", m->num_vertices, m->num_faces);
#endif
            for(i=0; i<m->num_vertices; ++i)
            {
                fprintf(fp, "%f %f %f %f %f %f %f %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z, m->vnormals[i].x,  m->vnormals[i].y, m->vnormals[i].z, m->vcolors[i].r,  m->vcolors[i].g, m->vcolors[i].b, m->vcolors[i].a);
            }
        }
        else
        {
            fprintf(fp, "OFF\n");
#if MESH_INTDATA_TYPE==0
            fprintf(fp, "%d %d 0\n", m->num_vertices, m->num_faces);
#else
            fprintf(fp, "%"PRId64" %"PRId64" 0\n", m->num_vertices, m->num_faces);
#endif
            for(i=0; i<m->num_vertices; ++i)
            {
                fprintf(fp, "%f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z);
            }
        }
        if(m->is_faces)
        {
            if(m->is_fcolors)
            {
                for(i=0; i<m->num_faces; ++i)
                {
#if MESH_INTDATA_TYPE==0
                    fprintf(fp, "%d", m->faces[i].num_vertices);
#else
                    fprintf(fp, "%"PRId64"", m->faces[i].num_vertices);
#endif
                    for(j=0; j<m->faces[i].num_vertices; ++j)
                    {
#if MESH_INTDATA_TYPE==0
                        fprintf(fp, " %d", m->faces[i].vertices[j]);
#else
                        fprintf(fp, " %"PRId64"", m->faces[i].vertices[j]);
#endif
                    }
                    fprintf(fp, " %f %f %f %f\n", m->fcolors[i].r, m->fcolors[i].g, m->fcolors[i].b, m->fcolors[i].a);
                }
            }
            else
            {
                for(i=0; i<m->num_faces; ++i)
                {
#if MESH_INTDATA_TYPE==0
                    fprintf(fp, "%d", m->faces[i].num_vertices);
#else
                    fprintf(fp, "%"PRId64"", m->faces[i].num_vertices);
#endif
                    for(j=0; j<m->faces[i].num_vertices; ++j)
                    {
#if MESH_INTDATA_TYPE==0
                        fprintf(fp, " %d", m->faces[i].vertices[j]);
#else
                        fprintf(fp, " %"PRId64"", m->faces[i].vertices[j]);
#endif
                    }
                    fprintf(fp, "\n");
                }
            }
        }
        fclose(fp);
    }
    return 0;
}

/** \brief Saves a mesh to an XYZ file
 *
 * \param[in] m Input mesh
 * \param[in] fname Output filename
 * \return Error code
 *
 */

int mesh_save_xyz(MESH m, const char* fname)
{
    INTDATA i;
    if(m->is_vertices)
    {
        FILEPOINTER fp = NULL;
        if((fp = fopen(fname,"wb"))==NULL) mesh_error(MESH_ERR_FNOTOPEN);
        for(i=0; i<m->num_vertices; ++i)
        {
            fprintf(fp, "%f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z);
        }
        fclose(fp);
    }
    return 0;
}

/** \brief Saves a mesh to a PLY file
 *
 * \param[in] m Input mesh
 * \param[in] fname Output filename
 * \return Error code
 *
 */

int mesh_save_ply(MESH m, const char* fname)
{
    INTDATA i, j;
    if(m->is_vertices)
    {
        FILEPOINTER fp = NULL;
        if((fp = fopen(fname,"wb"))==NULL) mesh_error(MESH_ERR_FNOTOPEN);
        fprintf(fp, "ply\n");
        fprintf(fp, "format ascii 1.0\n");
        if(m->is_vscalars)
        {
            if(m->is_vcolors && !m->is_vnormals)
            {
#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element vertex %d\n", m->num_vertices);
#else
                fprintf(fp, "element vertex %"PRId64"\n", m->num_vertices);
#endif
                fprintf(fp, "property float x\n");
                fprintf(fp, "property float y\n");
                fprintf(fp, "property float z\n");
                fprintf(fp, "property float red\n");
                fprintf(fp, "property float green\n");
                fprintf(fp, "property float blue\n");
                fprintf(fp, "property float alpha\n");
                fprintf(fp, "property float scalar_intensity\n");

#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element face %d\n", m->num_faces);
#else
                fprintf(fp, "element face %"PRId64"\n", m->num_faces);
#endif
                fprintf(fp, "property list uchar int vertex_index\n");
                fprintf(fp, "end_header\n");

                for(i=0; i<m->num_vertices; ++i)
                {
                    fprintf(fp, "%f %f %f %f %f %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z, m->vcolors[i].r,  m->vcolors[i].g, m->vcolors[i].b, m->vcolors[i].a, m->vscalars[i]);
                }
            }
            else if(m->is_vnormals && !m->is_vcolors)
            {
#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element vertex %d\n", m->num_vertices);
#else
                fprintf(fp, "element vertex %"PRId64"\n", m->num_vertices);
#endif
                fprintf(fp, "property float x\n");
                fprintf(fp, "property float y\n");
                fprintf(fp, "property float z\n");
                fprintf(fp, "property float nx\n");
                fprintf(fp, "property float ny\n");
                fprintf(fp, "property float nz\n");
                fprintf(fp, "property float scalar_intensity\n");

#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element face %d\n", m->num_faces);
#else
                fprintf(fp, "element face %"PRId64"\n", m->num_faces);
#endif
                fprintf(fp, "property list uchar int vertex_index\n");
                fprintf(fp, "end_header\n");


                for(i=0; i<m->num_vertices; ++i)
                {
                    fprintf(fp, "%f %f %f %f %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z, m->vnormals[i].x,  m->vnormals[i].y, m->vnormals[i].z, m->vscalars[i]);
                }
            }
            else if(m->is_vnormals && m->is_vcolors)
            {
#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element vertex %d\n", m->num_vertices);
#else
                fprintf(fp, "element vertex %"PRId64"\n", m->num_vertices);
#endif
                fprintf(fp, "property float x\n");
                fprintf(fp, "property float y\n");
                fprintf(fp, "property float z\n");
                fprintf(fp, "property float nx\n");
                fprintf(fp, "property float ny\n");
                fprintf(fp, "property float nz\n");
                fprintf(fp, "property float red\n");
                fprintf(fp, "property float green\n");
                fprintf(fp, "property float blue\n");
                fprintf(fp, "property float alpha\n");
                fprintf(fp, "property float scalar_intensity\n");

#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element face %d\n", m->num_faces);
#else
                fprintf(fp, "element face %"PRId64"\n", m->num_faces);
#endif
                fprintf(fp, "property list uchar int vertex_index\n");
                fprintf(fp, "end_header\n");


                for(i=0; i<m->num_vertices; ++i)
                {
                    fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z, m->vnormals[i].x,  m->vnormals[i].y, m->vnormals[i].z, m->vcolors[i].r,  m->vcolors[i].g, m->vcolors[i].b, m->vcolors[i].a, m->vscalars[i]);
                }
            }
            else
            {
#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element vertex %d\n", m->num_vertices);
#else
                fprintf(fp, "element vertex %"PRId64"\n", m->num_vertices);
#endif
                fprintf(fp, "property float x\n");
                fprintf(fp, "property float y\n");
                fprintf(fp, "property float z\n");
                fprintf(fp, "property float scalar_intensity\n");

#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element face %d\n", m->num_faces);
#else
                fprintf(fp, "element face %"PRId64"\n", m->num_faces);
#endif
                fprintf(fp, "property list uchar int vertex_index\n");
                fprintf(fp, "end_header\n");

                for(i=0; i<m->num_vertices; ++i)
                {
                    fprintf(fp, "%f %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z, m->vscalars[i]);
                }
            }

        }
        else
        {
            if(m->is_vcolors && !m->is_vnormals)
            {
#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element vertex %d\n", m->num_vertices);
#else
                fprintf(fp, "element vertex %"PRId64"\n", m->num_vertices);
#endif
                fprintf(fp, "property float x\n");
                fprintf(fp, "property float y\n");
                fprintf(fp, "property float z\n");
                fprintf(fp, "property float red\n");
                fprintf(fp, "property float green\n");
                fprintf(fp, "property float blue\n");
                fprintf(fp, "property float alpha\n");

#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element face %d\n", m->num_faces);
#else
                fprintf(fp, "element face %"PRId64"\n", m->num_faces);
#endif
                fprintf(fp, "property list uchar int vertex_index\n");
                fprintf(fp, "end_header\n");

                for(i=0; i<m->num_vertices; ++i)
                {
                    fprintf(fp, "%f %f %f %f %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z, m->vcolors[i].r,  m->vcolors[i].g, m->vcolors[i].b, m->vcolors[i].a);
                }
            }
            else if(m->is_vnormals && !m->is_vcolors)
            {
#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element vertex %d\n", m->num_vertices);
#else
                fprintf(fp, "element vertex %"PRId64"\n", m->num_vertices);
#endif
                fprintf(fp, "property float x\n");
                fprintf(fp, "property float y\n");
                fprintf(fp, "property float z\n");
                fprintf(fp, "property float nx\n");
                fprintf(fp, "property float ny\n");
                fprintf(fp, "property float nz\n");

#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element face %d\n", m->num_faces);
#else
                fprintf(fp, "element face %"PRId64"\n", m->num_faces);
#endif
                fprintf(fp, "property list uchar int vertex_index\n");
                fprintf(fp, "end_header\n");


                for(i=0; i<m->num_vertices; ++i)
                {
                    fprintf(fp, "%f %f %f %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z, m->vnormals[i].x,  m->vnormals[i].y, m->vnormals[i].z);
                }
            }
            else if(m->is_vnormals && m->is_vcolors)
            {
#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element vertex %d\n", m->num_vertices);
#else
                fprintf(fp, "element vertex %"PRId64"\n", m->num_vertices);
#endif
                fprintf(fp, "property float x\n");
                fprintf(fp, "property float y\n");
                fprintf(fp, "property float z\n");
                fprintf(fp, "property float nx\n");
                fprintf(fp, "property float ny\n");
                fprintf(fp, "property float nz\n");
                fprintf(fp, "property float red\n");
                fprintf(fp, "property float green\n");
                fprintf(fp, "property float blue\n");
                fprintf(fp, "property float alpha\n");

#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element face %d\n", m->num_faces);
#else
                fprintf(fp, "element face %"PRId64"\n", m->num_faces);
#endif
                fprintf(fp, "property list uchar int vertex_index\n");
                fprintf(fp, "end_header\n");


                for(i=0; i<m->num_vertices; ++i)
                {
                    fprintf(fp, "%f %f %f %f %f %f %f %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z, m->vnormals[i].x,  m->vnormals[i].y, m->vnormals[i].z, m->vcolors[i].r,  m->vcolors[i].g, m->vcolors[i].b, m->vcolors[i].a);
                }
            }
            else
            {
#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element vertex %d\n", m->num_vertices);
#else
                fprintf(fp, "element vertex %"PRId64"\n", m->num_vertices);
#endif
                fprintf(fp, "property float x\n");
                fprintf(fp, "property float y\n");
                fprintf(fp, "property float z\n");

#if MESH_INTDATA_TYPE==0
                fprintf(fp, "element face %d\n", m->num_faces);
#else
                fprintf(fp, "element face %"PRId64"\n", m->num_faces);
#endif
                fprintf(fp, "property list uchar int vertex_index\n");
                fprintf(fp, "end_header\n");

                for(i=0; i<m->num_vertices; ++i)
                {
                    fprintf(fp, "%f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z);
                }
            }
        }
        if(m->is_faces)
            for(i=0; i<m->num_faces; ++i)
            {
#if MESH_INTDATA_TYPE==0
                fprintf(fp, "%d", m->faces[i].num_vertices);
#else
                fprintf(fp, "%"PRId64"", m->faces[i].num_vertices);
#endif
                for(j=0; j<m->faces[i].num_vertices; ++j)
                {
#if MESH_INTDATA_TYPE==0
                    fprintf(fp, " %d", m->faces[i].vertices[j]);
#else
                    fprintf(fp, " %"PRId64"", m->faces[i].vertices[j]);
#endif
                }
                fprintf(fp, "\n");
            }
        fclose(fp);
    }
    return 0;
}

/** \brief Saves a mesh to a BINv1 file
 *
 * \param[in] m Input mesh
 * \param[in] fname Output filename
 * \return Error code
 *
 */

int mesh_save_bin(MESH m, const char* fname)
{
    if(m->is_vertices)
    {
        FILEPOINTER fp = NULL;
        INTDATA i;
        uint32_t buff32;
        uint8_t buff8;
        float xyz[3];
        uint8_t rgb[3];
        if((fp = fopen(fname,"wb"))==NULL) mesh_error(MESH_ERR_FNOTOPEN);
        buff32 = 1;
        fwrite(&buff32, sizeof(uint32_t), 1, fp); /* number of pcds */
        buff32 = m->num_vertices;
        fwrite(&buff32, sizeof(uint32_t), 1, fp); /* number of vertices */

        buff8 = 1+((m->is_vcolors)*2)+((m->is_vnormals)*4)+((m->is_vscalars)*8); /* flags */
        fwrite(&buff8, sizeof(uint8_t), 1, fp); /* number of vertices */

        if(!m->is_vcolors && !m->is_vnormals)
        {
            for(i=0; i<m->num_vertices; ++i)
            {
                const MESH_VERTEX cmv = m->vertices+i;
                xyz[0] = cmv->x;
                xyz[1] = cmv->y;
                xyz[2] = cmv->z;
                fwrite(xyz, sizeof(float), 3, fp);
            }
        }
        else if(m->is_vcolors && !m->is_vnormals)
        {
            for(i=0; i<m->num_vertices; ++i)
            {
                const MESH_VERTEX cmv = m->vertices+i;
                const MESH_COLOR cmvc = m->vcolors+i;
                xyz[0] = cmv->x;
                xyz[1] = cmv->y;
                xyz[2] = cmv->z;
                fwrite(xyz, sizeof(float), 3, fp);
                rgb[0] = cmvc->r;
                rgb[1] = cmvc->g;
                rgb[2] = cmvc->b;
                fwrite(rgb, sizeof(uint8_t), 3, fp);
            }
        }
        else if(!m->is_vcolors && m->is_vnormals)
        {
            for(i=0; i<m->num_vertices; ++i)
            {
                const MESH_VERTEX cmv = m->vertices+i;
                const MESH_NORMAL cmvn = m->vnormals+i;
                xyz[0] = cmv->x;
                xyz[1] = cmv->y;
                xyz[2] = cmv->z;
                fwrite(xyz, sizeof(float), 3, fp);
                xyz[0] = cmvn->x;
                xyz[1] = cmvn->y;
                xyz[2] = cmvn->z;
                fwrite(xyz, sizeof(float), 3, fp);
            }
        }
        else if(m->is_vcolors && m->is_vnormals)
        {
            for(i=0; i<m->num_vertices; ++i)
            {
                const MESH_VERTEX cmv = m->vertices+i;
                const MESH_COLOR cmvc = m->vcolors+i;
                const MESH_NORMAL cmvn = m->vnormals+i;
                xyz[0] = cmv->x;
                xyz[1] = cmv->y;
                xyz[2] = cmv->z;
                fwrite(xyz, sizeof(float), 3, fp);
                rgb[0] = cmvc->r;
                rgb[1] = cmvc->g;
                rgb[2] = cmvc->b;
                fwrite(rgb, sizeof(uint8_t), 3, fp);
                xyz[0] = cmvn->x;
                xyz[1] = cmvn->y;
                xyz[2] = cmvn->z;
                fwrite(xyz, sizeof(float), 3, fp);
            }
        }
        fclose(fp);
    }
    return 0;
}

/** \brief Saves a mesh to an OBJ file
 *
 * \param[in] m Input mesh
 * \param[in] fname Output filename
 * \return Error code
 *
 */

int mesh_save_obj(MESH m, const char* fname)
{
    INTDATA i, j;
    if(m->is_vertices)
    {
        FILEPOINTER fp = NULL;
        if((fp = fopen(fname,"wb"))==NULL) mesh_error(MESH_ERR_FNOTOPEN);
        for(i=0; i<m->num_vertices; ++i)
        {
            fprintf(fp,"v %f %f %f\n", m->vertices[i].x, m->vertices[i].y, m->vertices[i].z);
        }
        if(m->is_vnormals)
        {
            for(i=0; i<m->num_vertices; ++i)
            {
                fprintf(fp,"vn %f %f %f\n", m->vnormals[i].x,  m->vnormals[i].y, m->vnormals[i].z);
            }
        }

        for(i=0; i<m->num_faces; ++i)
        {
            fprintf(fp, "f ");
            for(j=0; j<m->faces[i].num_vertices; ++j)
            {
#if MESH_INTDATA_TYPE==0
                fprintf(fp, " %d", m->faces[i].vertices[j]+1);
#else
                fprintf(fp, " %"PRId64"", m->faces[i].vertices[j]+1);
#endif
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    return 0;
}



