/**
 * @file meshload.c
 * @author Sk. Mohammadul Haque
 * @version 1.4.2.0
 * @copyright
 * Copyright (c) 2013, 2014, 2015, 2016 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to loading different mesh file types.
 */

#include <string.h>
#include "../include/meshlib.h"

/** \brief Reads a mesh from an OFF/PLY/ASC/XYZ file
 *
 * \param[in] fname Input filename
 * \return Output mesh
 *
 */

MESH mesh_load_file(const char* fname)
{
    char *ext = strrchr(fname, '.');
    if(ext!=NULL)
    {
        if(strcmp(ext,".off")==0) return mesh_load_off(fname);
        else if(strcmp(ext,".ply")==0) return mesh_load_ply(fname);
        else if(strcmp(ext,".asc")==0) return mesh_load_xyz(fname);
        else if(strcmp(ext,".xyz")==0) return mesh_load_xyz(fname);
    }
    mesh_error(MESH_ERR_FNOTOPEN);
    return NULL;
}

/** \brief Reads a mesh from an OFF file
 *
 * \param[in] fname Input filename
 * \return Output mesh
 *
 */

MESH mesh_load_off(const char* fname)
{
    FILEPOINTER fp = NULL;
    MESH m = NULL;
    if((fp = fopen(fname,"rb"))==NULL) mesh_error(MESH_ERR_FNOTOPEN);
    m = mesh_create_mesh_new();
    m = __mesh_parse_off_header(m, fp);
    if(m->num_vertices>0) m = __mesh_parse_off_vertices(m, fp);
    if(m->num_faces>0 && m->is_vertices) m = __mesh_parse_off_faces(m, fp);
    if(m->is_vertices) m->is_loaded = 1;
    fclose(fp);
    return m;
}

/** \cond HIDDEN_SYMBOLS */

MESH __mesh_parse_off_header(MESH m, FILEPOINTER fp)
{
    char dummy[16];
    int flag;
    do
    {
        flag = mesh_read_word_only(fp, dummy, 16);
    }
    while(flag==3);
    if(flag>0) return m;
    if(strcmp(dummy, "OFF")==0) m->origin_type = MESH_ORIGIN_TYPE_OFF;
    else if(strcmp(dummy, "NOFF")==0) m->origin_type = MESH_ORIGIN_TYPE_NOFF;
    else if(strcmp(dummy, "COFF")==0) m->origin_type = MESH_ORIGIN_TYPE_COFF;
    else if(strcmp(dummy, "NCOFF")==0) m->origin_type = MESH_ORIGIN_TYPE_NCOFF;
    else return m;
    m->is_loaded = 0;
    do
    {
        flag = mesh_read_word_only(fp, dummy, 16);
    }
    while(flag==3);
    if(flag>0) return m;
    m->num_vertices = strtol(dummy, NULL, 0);

    if(mesh_read_word_only(fp, dummy, 16)>0) return m;
    m->num_faces = strtol(dummy, NULL, 0);

    if(mesh_read_word_only(fp, dummy, 16)>0) return m;
    /* ignore edges */
    return m;
}

MESH __mesh_parse_off_vertices(MESH m, FILEPOINTER fp)
{
    INTDATA i;
    if(m->is_vertices) free(m->vertices);
    if((m->vertices = (MESH_VERTEX)malloc(sizeof(mesh_vertex)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    m->is_vertices = 1;

    switch(m->origin_type)
    {
    case MESH_ORIGIN_TYPE_OFF:
        for(i=0; i<m->num_vertices; ++i)
        {
#if MESH_FLOATDATA_TYPE==0
            if(fscanf(fp, " %f %f %f \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z)!=3)
#else
            if(fscanf(fp, " %lf %lf %lf \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z)!=3)
#endif
            {
                free(m->vertices);
                m->is_vertices = 0;
                return m;
            }
        }
        break;

    case MESH_ORIGIN_TYPE_COFF:
        if((m->vcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        m->is_vcolors = 1;
        for(i=0; i<m->num_vertices; ++i)
        {
#if MESH_FLOATDATA_TYPE==0
            if(fscanf(fp, " %f %f %f %f %f %f %f \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vcolors[i].r, &m->vcolors[i].g, &m->vcolors[i].b, &m->vcolors[i].a)!=7)
#else
            if(fscanf(fp, " %lf %lf %lf %lf %lf %lf %lf \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vcolors[i].r, &m->vcolors[i].g, &m->vcolors[i].b, &m->vcolors[i].a)!=7)
#endif
            {
                free(m->vertices);
                free(m->vcolors);
                m->is_vertices = 0;
                m->is_vcolors = 0;
                return m;
            }
        }
        break;

    case MESH_ORIGIN_TYPE_NOFF:
        if((m->vnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        m->is_vnormals = 1;
        for(i=0; i<m->num_vertices; ++i)
        {
#if MESH_FLOATDATA_TYPE==0
            if(fscanf(fp, " %f %f %f %f %f %f \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vnormals[i].x, &m->vnormals[i].y, &m->vnormals[i].z)!=6)
#else
            if(fscanf(fp, " %lf %lf %lf %lf %lf %lf \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vnormals[i].x, &m->vnormals[i].y, &m->vnormals[i].z)!=6)
#endif
            {
                free(m->vertices);
                free(m->vnormals);
                m->is_vertices = 0;
                m->is_vnormals = 0;
                return m;
            }
        }
        break;

    case MESH_ORIGIN_TYPE_NCOFF:
        if((m->vnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        if((m->vcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        m->is_vnormals = 1;
        m->is_vcolors = 1;
        for(i=0; i<m->num_vertices; ++i)
        {
#if MESH_FLOATDATA_TYPE==0
            if(fscanf(fp, " %f %f %f %f %f %f %f %f %f %f \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vnormals[i].x, &m->vnormals[i].y, &m->vnormals[i].z, &m->vcolors[i].r, &m->vcolors[i].g, &m->vcolors[i].b, &m->vcolors[i].a)!=10)
#else
            if(fscanf(fp, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vnormals[i].x, &m->vnormals[i].y, &m->vnormals[i].z, &m->vcolors[i].r, &m->vcolors[i].g, &m->vcolors[i].b, &m->vcolors[i].a)!=10)
#endif
            {
                free(m->vertices);
                free(m->vcolors);
                m->is_vertices = 0;
                m->is_vcolors = 0;
                return m;
            }
        }
        break;

    }

    return m;
}

MESH __mesh_parse_off_faces(MESH m, FILEPOINTER fp)
{
    INTDATA i, j, nv;
    char dummy[32];
    long int currpos;
    int nwrds = 0, flag = 0;
    INTDATA nverts = 0;
    if(m->is_faces) free(m->faces);
    if((m->faces = (MESH_FACE)malloc(sizeof(mesh_face)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
    m->is_faces = 1;
    m->is_trimesh = 1;
    currpos = ftell(fp);
    mesh_count_words_in_line(fp, &nwrds);
    fseek(fp, currpos, SEEK_SET);
#if MESH_INTDATA_TYPE==0
    if(fscanf(fp, "%d", &nverts)!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#else
    if(fscanf(fp, "%ld", &nverts)!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#endif
    fseek(fp, currpos, SEEK_SET);
    if(nwrds==(nverts+1))
    {
        m->is_fcolors = 0;
        for(i=0; i<m->num_faces; ++i)
        {
            do
            {
                flag = mesh_read_word_only(fp, dummy, 32);
            }
            while(flag==3);
            if(flag>0)
            {
                free(m->faces);
                m->is_faces = 0;
                return m;
            }
            nv = strtol(dummy, NULL, 0);
            m->faces[i].num_vertices = nv;
            if(nv!=3) m->is_trimesh = 0;
            if((m->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*nv))==NULL) mesh_error(MESH_ERR_MALLOC);
            for(j=0; j<nv; ++j)
            {
#if MESH_INTDATA_TYPE==0
                if(fscanf(fp, "%d", &(m->faces[i].vertices[j]))!=1)
#else
                if(fscanf(fp, "%ld", &(m->faces[i].vertices[j]))!=1)
#endif
                {
                    free(m->faces);
                    m->is_faces = 0;
                    return m;
                }
            }
        }
    }
    else if(nwrds==(nverts+2))
    {
        if(m->is_fcolors) free(m->fcolors);
        if((m->fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        m->is_faces = 1;
        m->is_fcolors = 1;
        for(i=0; i<m->num_faces; ++i)
        {
            do
            {
                flag = mesh_read_word_only(fp, dummy, 32);
            }
            while(flag==3);
            if(flag>0)
            {
                free(m->faces);
                m->is_faces = 0;
                return m;
            }
            nv = strtol(dummy, NULL, 0);
            m->faces[i].num_vertices = nv;
            if(nv!=3) m->is_trimesh = 0;
            if((m->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*nv))==NULL) mesh_error(MESH_ERR_MALLOC);
            for(j=0; j<nv; ++j)
            {
#if MESH_INTDATA_TYPE==0
                if(fscanf(fp, "%d", &(m->faces[i].vertices[j]))!=1)
#else
                if(fscanf(fp, "%ld", &(m->faces[i].vertices[j]))!=1)
#endif
                {
                    mesh_error(MESH_ERR_SIZE_MISMATCH);
                }
            }
#if MESH_FLOATDATA_TYPE==0
            if(fscanf(fp, "%f", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#else
            if(fscanf(fp, "%lf", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#endif
            m->fcolors[i].g = m->fcolors[i].r;
            m->fcolors[i].b = m->fcolors[i].r;
            m->fcolors[i].a = 1.0f;
        }
    }
    else if(nwrds==(nverts+4))
    {
        if(m->is_fcolors) free(m->fcolors);
        if((m->fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        m->is_faces = 1;
        m->is_fcolors = 1;
        for(i=0; i<m->num_faces; ++i)
        {
            do
            {
                flag = mesh_read_word_only(fp, dummy, 32);
            }
            while(flag==3);
            if(flag>0)
            {
                free(m->faces);
                m->is_faces = 0;
                return m;
            }
            nv = strtol(dummy, NULL, 0);
            m->faces[i].num_vertices = nv;
            if(nv!=3) m->is_trimesh = 0;
            if((m->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*nv))==NULL) mesh_error(MESH_ERR_MALLOC);
            for(j=0; j<nv; ++j)
            {
#if MESH_INTDATA_TYPE==0
                if(fscanf(fp, "%d", &(m->faces[i].vertices[j]))!=1)
#else
                if(fscanf(fp, "%ld", &(m->faces[i].vertices[j]))!=1)
#endif
                {
                    mesh_error(MESH_ERR_SIZE_MISMATCH);
                }
            }
#if MESH_FLOATDATA_TYPE==0
            if(fscanf(fp, "%f", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
            if(fscanf(fp, "%f", &(m->fcolors[i].g))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
            if(fscanf(fp, "%f", &(m->fcolors[i].b))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#else
            if(fscanf(fp, "%lf", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
            if(fscanf(fp, "%lf", &(m->fcolors[i].g))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
            if(fscanf(fp, "%lf", &(m->fcolors[i].b))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#endif
            m->fcolors[i].a = 1.0f;
        }
    }
    else
    {
        if(m->is_fcolors) free(m->fcolors);
        if((m->fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        m->is_faces = 1;
        m->is_fcolors = 1;
        for(i=0; i<m->num_faces; ++i)
        {
            do
            {
                flag = mesh_read_word_only(fp, dummy, 32);
            }
            while(flag==3);
            if(flag>0)
            {
                free(m->faces);
                m->is_faces = 0;
                return m;
            }
            nv = strtol(dummy, NULL, 0);
            m->faces[i].num_vertices = nv;
            if(nv!=3) m->is_trimesh = 0;
            if((m->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*nv))==NULL) mesh_error(MESH_ERR_MALLOC);
            for(j=0; j<nv; ++j)
            {
#if MESH_INTDATA_TYPE==0
                if(fscanf(fp, "%d", &(m->faces[i].vertices[j]))!=1)
#else
                if(fscanf(fp, "%ld", &(m->faces[i].vertices[j]))!=1)
#endif
                {
                    mesh_error(MESH_ERR_SIZE_MISMATCH);
                }
            }
#if MESH_FLOATDATA_TYPE==0
            if(fscanf(fp, "%f", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
            if(fscanf(fp, "%f", &(m->fcolors[i].g))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
            if(fscanf(fp, "%f", &(m->fcolors[i].b))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
            if(fscanf(fp, "%f", &(m->fcolors[i].a))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#else
            if(fscanf(fp, "%lf", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
            if(fscanf(fp, "%lf", &(m->fcolors[i].g))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
            if(fscanf(fp, "%lf", &(m->fcolors[i].b))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
            if(fscanf(fp, "%lf", &(m->fcolors[i].a))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#endif
        }
    }
    return m;
}

/** \endcond */

/** \brief Read a mesh from an ASC/XYZ file
 *
 * \param[in] fname Input filename
 * \return Output mesh
 *
 */

MESH mesh_load_xyz(const char* fname)
{
    FILEPOINTER fp = NULL;
    MESH m = NULL;
    if((fp = fopen(fname,"rb"))==NULL) mesh_error(MESH_ERR_FNOTOPEN);
    m = mesh_create_mesh_new();
    m = __mesh_parse_xyz_data(m, fp);
    if(m->is_vertices) m->is_loaded = 1;
    fclose(fp);
    return m;
}

/** \cond HIDDEN_SYMBOLS */

MESH __mesh_parse_xyz_data(MESH m, FILEPOINTER fp)
{
    INTDATA n1 = 0, n2 = 0, i, flag = 0, tmp = 0, k = 0;
    FLOATDATA in_value = 0;
    char c_word[100];

    while(!flag)
    {
        k = mesh_isnumeric(fp);
        if(k==1)
        {
            flag = mesh_count_words_in_line(fp, &tmp);
            if(tmp==3) ++n2; /*improved skip non-data line*/
            else if(tmp>0) mesh_error(MESH_ERR_SIZE_MISMATCH);
        }
        else if(k==0)
        {
            flag = mesh_count_words_in_line(fp, &tmp);
        }
        else flag = 1;
    }
    n1 = tmp;
    if(n1==0 || n2==0) mesh_error(MESH_ERR_SIZE_MISMATCH);
    rewind(fp);

    m->num_vertices = n2;
    if(m->is_vertices) free(m->vertices);
    if((m->vertices = (MESH_VERTEX)malloc(sizeof(mesh_vertex)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    m->is_vertices = 1;
    for (i=0; i<n2; ++i)
    {
        k = mesh_isnumeric(fp);
        if(k==0)
        {
            do
            {
                flag = mesh_count_words_in_line(fp, &tmp);
                k = mesh_isnumeric(fp);
            }
            while(k==0);
        }
        do
        {
            flag = mesh_read_word_only(fp, c_word, 100);
        }
        while(flag==3);
        in_value = (float)strtod(c_word, NULL);
        in_value = strtod(c_word, NULL);
        m->vertices[i].x = in_value;

        mesh_read_word_only(fp, c_word, 100);
        in_value = (float)strtod(c_word, NULL);
        in_value = strtod(c_word, NULL);
        m->vertices[i].y = in_value;

        mesh_read_word_only(fp, c_word, 100);
        in_value = (float)strtod(c_word, NULL);
        in_value = strtod(c_word, NULL);
        m->vertices[i].z = in_value;
    }
    return m;
}

struct plydsizet
{
    char* str;
    int t;
};
typedef struct plydsizet plysizet;

plysizet plydsizes[] = {{"char",1}, /* 0 */
    {"uchar",1},
    {"short",2}, /* 2 */
    {"ushort",2},
    {"int",4}, /* 4 */
    {"uint",4},
    {"float",4}, /* 6 */
    {"float32",4},
    {"double",8}
}; /* 8 */

static int __mesh_checkdtype(char* str)
{
    int i;
    for(i=0; i<9; ++i)
    {
        if(strcmp(plydsizes[i].str, str)==0) return i;
    }
    return -1;
}

static INTDATA __mesh_convert_format_int_0(int dtype, char* in)
{
    switch(dtype)
    {
    case 0:
        return *in;
    case 1:
        return *((unsigned char *)in);
    case 2:
        return *((short *)in);
    case 3:
        return *((unsigned short *)in);
    case 4:
        return *((int *)in);
    case 5:
        return *((unsigned int *)in);
    default:
        mesh_error(MESH_ERR_UNKNOWN);
    }
    return -1;
}

static INTDATA __mesh_convert_format_int_1(int dtype, char* in)
{
    switch(dtype)
    {
    case 0:
        return *((char *)in);
    case 1:
        return *((unsigned char *)in);
    case 2:
        return (short)in[0]|(short)(in[1]<<8);
    case 3:
        return (unsigned short)in[0]|((unsigned short)in[1]<<8);
    case 4:
        return (int)in[0]|((int)in[1]<<8)|((int)in[2]<<16)|((int)in[3]<<24);
    case 5:
        return (unsigned int)in[0]|((unsigned int)in[1]<<8)|((unsigned int)in[2]<<16)|((unsigned int)in[3]<<24);
    default:
        mesh_error(MESH_ERR_UNKNOWN);
    }
    return -1;
}

static __inline float __mesh_convertf(char* in)
{
    float out;
    char* pout = (char*) &out;
    pout[0] = in[3];
    pout[1] = in[2];
    pout[2] = in[1];
    pout[3] = in[0];
    return out;
}

static __inline double __mesh_convertd(char* in)
{
    double out;
    char* pout = (char*) &out;
    pout[0] = in[7];
    pout[1] = in[6];
    pout[2] = in[5];
    pout[3] = in[4];
    pout[4] = in[3];
    pout[5] = in[2];
    pout[6] = in[1];
    pout[7] = in[0];
    return out;
}

static FLOATDATA __mesh_convert_format_float_0(int dtype, char* in)
{
    switch(dtype)
    {
    case 6:
    case 7:
        return *((float *)in);
    case 8:
        return *((double *)in);
    default:
        mesh_error(MESH_ERR_UNKNOWN);
    }
    return -1;
}

static FLOATDATA __mesh_convert_format_float_1(int dtype, char* in)
{
    switch(dtype)
    {
    case 6:
    case 7:
        return __mesh_convertf(in);
    case 8:
        return __mesh_convertd(in);
    default:
        mesh_error(MESH_ERR_UNKNOWN);
    }
    return -1;
}

static void __mesh_convert_format_mixed(int* dtypes, int*cnte, int l, char* in, INTDATA* outi, FLOATDATA* outf, int cnt, int type)
{
    int i, j, k;
    int ci = 0, cf = 0, c = 0;
    if(type==0)
    {
        for(i=0; i<cnt; ++i)
        {
            for(j=0; j<l; ++j)
            {
                for(k=0; k<cnte[j]; ++k)
                {
                    if(dtypes[j]<6)
                    {
                        outi[ci++] = __mesh_convert_format_int_0(dtypes[j], in+c);
                    }
                    else
                    {
                        outf[cf++] = __mesh_convert_format_float_0(dtypes[j], in+c);
                    }
                    c += plydsizes[dtypes[j]].t;

                }
            }

        }
    }
    else
    {
        for(i=0; i<cnt; ++i)
        {
            for(j=0; j<l; ++j)
            {
                for(k=0; k<cnte[j]; ++k)
                {
                    if(dtypes[j]<6)
                    {
                        outi[ci++] = __mesh_convert_format_int_1(dtypes[j], in+c);
                    }
                    else
                    {
                        outf[cf++] = __mesh_convert_format_float_1(dtypes[j], in+c);
                    }
                    c += plydsizes[dtypes[j]].t;

                }
            }
        }
    }
}


#define __v 0
#define __vn 1
#define __vc 2
#define __ft 3
#define __f 4
#define __fnt 5
#define __fn 6
#define __fct 7
#define __fct1 8
#define __endian 9
#define __doffset 10

/** \endcond */

/** \brief Reads a mesh from a PLY file
 *
 * \param[in] fname Input filename
 * \return Output mesh
 *
 */

MESH mesh_load_ply(const char* fname)
{
    FILEPOINTER fp = NULL;
    unsigned int i = 1;
    char *c = (char*)&i;
    int dtypes[20] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0}; /* v, vn, vc, ft, f, fnt, fn, fct, fc */
    MESH m = NULL;
    if((fp = fopen(fname,"rb"))==NULL) mesh_error(MESH_ERR_FNOTOPEN);
    m = mesh_create_mesh_new();
    /* check this machine endianness */
    if(*c) dtypes[__endian] = MESH_ORIGIN_TYPE_PLY_BINARY_LITTLE_ENDIAN;
    else dtypes[__endian] = MESH_ORIGIN_TYPE_PLY_BINARY_BIG_ENDIAN;
    m = __mesh_parse_ply_header(m, fp, dtypes);
    m = __mesh_parse_ply_body(m, fp, dtypes);
    if(m->is_vertices) m->is_loaded = 1;
    else
    {
        mesh_free_mesh(m);
        m = NULL;
    }
    fclose(fp);
    return m;
}

/** \cond HIDDEN_SYMBOLS */

MESH __mesh_parse_ply_header(MESH m, FILEPOINTER fp, int* dtypes)
{
    char dummy[32], dummy2[32];
    int flag, element_done = 0;
    do
    {
        flag = mesh_read_word_only(fp, dummy, 32);
        if(strcmp(dummy, "comment")==0)
        {
            mesh_skip_line(fp);
            flag = 3;
        }
    }
    while(flag==3);
    if(flag>0) return m;
    if(strcmp(dummy, "ply")==0) m->origin_type = MESH_ORIGIN_TYPE_PLY_ASCII;
    else return m;
    m->is_loaded = 0;
    do
    {
        flag = mesh_read_word_only(fp, dummy, 32);
        if(strcmp(dummy, "comment")==0)
        {
            mesh_skip_line(fp);
            flag = 3;
        }
    }
    while(flag==3);
    if(flag>0) return m;
    if(strcmp(dummy, "format")==0)
    {
        flag = mesh_read_word_only(fp, dummy, 32);
        if(strcmp(dummy, "ascii")==0) m->origin_type = MESH_ORIGIN_TYPE_PLY_ASCII;
        else if(strcmp(dummy, "binary_little_endian")==0) m->origin_type = MESH_ORIGIN_TYPE_PLY_BINARY_LITTLE_ENDIAN;
        else if(strcmp(dummy, "binary_big_endian")==0) m->origin_type = MESH_ORIGIN_TYPE_PLY_BINARY_BIG_ENDIAN;
        else mesh_error(MESH_ERR_UNKNOWN);
        mesh_skip_line(fp);
    }
    else return m;
    do
    {
        flag = mesh_read_word_only(fp, dummy, 32);
        if(strcmp(dummy, "comment")==0)
        {
            mesh_skip_line(fp);
            flag = 3;
        }
    }
    while(flag==3);
    if(strcmp(dummy, "element")==0)
    {
        mesh_read_word_only(fp, dummy, 32);
        if(strcmp(dummy, "vertex")==0)
        {
            if(mesh_read_word_only(fp, dummy, 32)>0) return m;
            m->num_vertices = strtol(dummy, NULL, 0);
            element_done = 0;
            do
            {
                do
                {
                    flag = mesh_read_word_only(fp, dummy, 32);
                    if(strcmp(dummy, "element")==0) break;
                    if(strcmp(dummy, "comment")==0)
                    {
                        mesh_skip_line(fp);
                        flag = 3;
                    }
                }
                while(flag==3); /* v, vn, vc, ft, f, fnt, fn, fct, fc */
                if(strcmp(dummy, "element")==0) break;
                if(strcmp(dummy, "property")==0)
                {
                    mesh_read_word_only(fp, dummy2, 32);
                    mesh_read_word_only(fp, dummy, 32);
                    if(strcmp(dummy, "red")==0)
                    {
                        m->is_vcolors = 1;
                        dtypes[__vc] = __mesh_checkdtype(dummy2);
                    }
                    else if(strcmp(dummy, "nx")==0)
                    {
                        m->is_vnormals = 1;
                        dtypes[__vn] = __mesh_checkdtype(dummy2);
                    }
                    else if(strcmp(dummy, "x")==0)
                    {
                        dtypes[__v] = __mesh_checkdtype(dummy2);
                    }
                }
                else element_done = 1;
            }
            while(element_done==0);
        }
    }
    if(strcmp(dummy, "element")==0)
    {
        mesh_read_word_only(fp, dummy, 32);
        if(strcmp(dummy, "face")==0)
        {
            if(mesh_read_word_only(fp, dummy, 32)>0) return m;
            m->num_faces = strtol(dummy, NULL, 0);
            element_done = 0;
            do
            {
                do
                {
                    flag = mesh_read_word_only(fp, dummy, 32);
                    if(strcmp(dummy, "comment")==0)
                    {
                        mesh_skip_line(fp);
                        flag = 3;
                    }
                }
                while(flag==3); /* v, vn, vc, ft, f, fnt, fn, fct, fc */
                if(strcmp(dummy, "property")==0)
                {
                    mesh_read_word_only(fp, dummy, 32);
                    if(strcmp(dummy, "list")==0)
                    {
                        int tmp1, tmp2;
                        mesh_read_word_only(fp, dummy2, 32);
                        tmp1 = __mesh_checkdtype(dummy2);
                        mesh_read_word_only(fp, dummy2, 32);
                        tmp2 = __mesh_checkdtype(dummy2);
                        mesh_read_word_only(fp, dummy, 32);
                        if((strcmp(dummy, "vertex_indices")==0)||(strcmp(dummy, "vertex_index")==0))
                        {
                            dtypes[__ft] = tmp1;
                            dtypes[__f] = tmp2;
                        }
                        else
                        {
                            m->num_faces = 0;
                            return m;
                        }
                        /* if(strcmp(dummy, "rgb")==0) m->is_fcolors = 1; */
                        /* if(strcmp(dummy, "nx")==0) m->is_fnormals = 1; */
                    }
                    else mesh_skip_line(fp);
                }
                else element_done = 1;
            }
            while(element_done==0);
        }
    }
    while(strcmp(dummy, "end_header")!=0)
    {
        flag = mesh_read_word_only(fp, dummy, 32);
        if(flag==1) mesh_error(MESH_ERR_UNKNOWN);
    }
    mesh_skip_line(fp);
    return m;
}

MESH __mesh_parse_ply_body(MESH m, FILEPOINTER fp, int* dtypes)
{
    __mesh_parse_ply_vertices(m, fp, dtypes);
    __mesh_parse_ply_faces(m, fp, dtypes);
    return m;
}

MESH __mesh_parse_ply_vertices(MESH m, FILEPOINTER fp, int *dtypes)
{
    INTDATA i, j, j2, tmp[16]; /* v, vn, vc, ft, f, fnt, fn, fct, fc */
    char* buff;
    INTDATA *obuff2;
    FLOATDATA* obuff1;
    if(m->is_vertices) free(m->vertices);
    if((m->vertices = (MESH_VERTEX)malloc(sizeof(mesh_vertex)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
    m->is_vertices = 1;

    if(!m->is_vcolors && !m->is_vnormals)
    {
        if(m->origin_type==MESH_ORIGIN_TYPE_PLY_ASCII)
        {
            for(i=0; i<m->num_vertices; ++i)
            {
#if MESH_FLOATDATA_TYPE==0
                if(fscanf(fp, " %f %f %f \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z)!=3)
#else
                if(fscanf(fp, " %lf %lf %lf \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z)!=3)
#endif
                {
                    free(m->vertices);
                    m->num_vertices = 0;
                    m->is_vertices = 0;
                    return m;
                }
            }
        }
        else
        {
            INTDATA csz = m->num_vertices*plydsizes[dtypes[__v]].t*3;
            if((buff = (char*)malloc(csz))==NULL) mesh_error(MESH_ERR_MALLOC);
            if(fread(buff, 1, csz, fp)!=(size_t)csz)
            {
                free(m->vertices);
                free(buff);
                m->is_vertices = 0;
                return m;
            }
            tmp[0] = 3;
            if((obuff1 = (FLOATDATA*)malloc(m->num_vertices*sizeof(FLOATDATA)*3))==NULL) mesh_error(MESH_ERR_MALLOC);
            if(m->origin_type==dtypes[__endian]) __mesh_convert_format_mixed(&dtypes[__v], tmp, 1, buff, NULL, obuff1, m->num_vertices, 0);
            else __mesh_convert_format_mixed(&dtypes[__v], tmp, 1, buff, NULL, obuff1, m->num_vertices, 1);
            j = 0;
            for(i=0; i<m->num_vertices; ++i, j+=3)
            {
                m->vertices[i].x = obuff1[j];
                m->vertices[i].y = obuff1[j+1];
                m->vertices[i].z = obuff1[j+2];
            }
            free(buff);
            free(obuff1);
        }
    }
    else if(m->is_vcolors && !m->is_vnormals)
    {
        if((m->vcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        if(m->origin_type==MESH_ORIGIN_TYPE_PLY_ASCII)
        {
            for(i=0; i<m->num_vertices; ++i)
            {
#if MESH_FLOATDATA_TYPE==0
                if(fscanf(fp, " %f %f %f %f %f %f %f \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vcolors[i].r, &m->vcolors[i].g, &m->vcolors[i].b, &m->vcolors[i].a)!=7)
#else
                if(fscanf(fp, " %lf %lf %lf %lf %lf %lf %lf \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vcolors[i].r, &m->vcolors[i].g, &m->vcolors[i].b, &m->vcolors[i].a)!=7)
#endif
                {
                    free(m->vertices);
                    free(m->vcolors);
                    m->is_vertices = 0;
                    m->is_vcolors = 0;
                    return m;
                }
            }
        }
        else
        {
            INTDATA csz = m->num_vertices*(plydsizes[dtypes[__v]].t*3+plydsizes[dtypes[__vc]].t*4);
            if((buff = (char*)malloc(csz))==NULL) mesh_error(MESH_ERR_MALLOC);
            if(fread(buff, 1, csz, fp)!=(size_t)csz)
            {
                free(m->vertices);
                free(m->vcolors);
                free(buff);
                m->is_vertices = 0;
                m->is_vcolors = 0;
                return m;
            }
            tmp[0] = dtypes[0];
            tmp[1] = dtypes[2];
            tmp[2] = 3;
            tmp[3] = 4; /* check if 3 or 4 does this refer to the number of color components, then include alpha also ==4 */
            if(dtypes[__vc]<6)
            {
                if((obuff1 = (FLOATDATA*)malloc(m->num_vertices*sizeof(FLOATDATA)*3))==NULL) mesh_error(MESH_ERR_MALLOC);
                if((obuff2 = (INTDATA*)malloc(m->num_vertices*sizeof(INTDATA)*4))==NULL) mesh_error(MESH_ERR_MALLOC);
                if(m->origin_type==dtypes[__endian]) __mesh_convert_format_mixed(tmp, tmp+2, 2, buff, obuff2, obuff1, m->num_vertices, 0);
                else __mesh_convert_format_mixed(tmp, tmp+2, 2, buff, obuff2, obuff1, m->num_vertices, 1);
                j = 0;
                j2 = 0;
                for(i=0; i<m->num_vertices; ++i, j+=3, j2+=4)
                {
                    m->vertices[i].x = obuff1[j];
                    m->vertices[i].y = obuff1[j+1];
                    m->vertices[i].z = obuff1[j+2];
                    m->vcolors[i].r = obuff2[j2];
                    m->vcolors[i].g = obuff2[j2+1];
                    m->vcolors[i].b = obuff2[j2+2];
                    m->vcolors[i].a = obuff2[j2+3];
                }
                free(obuff1);
                free(obuff2);
            }
            else
            {
                if((obuff1 = (FLOATDATA*)malloc(m->num_vertices*sizeof(FLOATDATA)*7))==NULL) mesh_error(MESH_ERR_MALLOC);
                if(m->origin_type==dtypes[__endian]) __mesh_convert_format_mixed(tmp, tmp+2, 2, buff, NULL, obuff1, m->num_vertices, 0);
                else __mesh_convert_format_mixed(tmp, tmp+2, 2, buff, NULL, obuff1, m->num_vertices, 1);
                j = 0;
                for(i=0; i<m->num_vertices; ++i, j+=6)
                {
                    m->vertices[i].x = obuff1[j];
                    m->vertices[i].y = obuff1[j+1];
                    m->vertices[i].z = obuff1[j+2];
                    m->vcolors[i].r = obuff1[j+3];
                    m->vcolors[i].g = obuff1[j+4];
                    m->vcolors[i].b = obuff1[j+5];
                    m->vcolors[i].a = obuff1[j+6];
                }
                free(obuff1);
            }
            free(buff);
        }
    }
    else if(!m->is_vcolors && m->is_vnormals)
    {
        if((m->vnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        if(m->origin_type==MESH_ORIGIN_TYPE_PLY_ASCII)
        {
            for(i=0; i<m->num_vertices; ++i)
            {
#if MESH_FLOATDATA_TYPE==0
                if(fscanf(fp, " %f %f %f %f %f %f \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vnormals[i].x, &m->vnormals[i].y, &m->vnormals[i].z)!=6)
#else
                if(fscanf(fp, " %lf %lf %lf %lf %lf %lf \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vnormals[i].x, &m->vnormals[i].y, &m->vnormals[i].z)!=6)
#endif
                {
                    free(m->vertices);
                    free(m->vnormals);
                    m->is_vertices = 0;
                    m->is_vnormals = 0;
                    return m;
                }
            }
        }
        else
        {
            INTDATA csz = m->num_vertices*(plydsizes[dtypes[__v]].t+plydsizes[dtypes[__vn]].t)*3;
            if((buff = (char*)malloc(csz))==NULL) mesh_error(MESH_ERR_MALLOC);
            if(fread(buff, 1, csz, fp)!=(size_t)csz)
            {
                free(m->vertices);
                free(m->vnormals);
                free(buff);
                m->is_vertices = 0;
                m->is_vnormals = 0;
                return m;
            }
            tmp[0] = dtypes[__v];
            tmp[1] = dtypes[__vn];
            tmp[2] = 3;
            tmp[3] = 3;
            if((obuff1 = (FLOATDATA*)malloc(m->num_vertices*sizeof(FLOATDATA)*6))==NULL) mesh_error(MESH_ERR_MALLOC);
            if(m->origin_type==dtypes[__endian]) __mesh_convert_format_mixed(tmp, tmp+2, 2, buff, NULL, obuff1, m->num_vertices, 0);
            else __mesh_convert_format_mixed(tmp, tmp+2, 2, buff, NULL, obuff1, m->num_vertices, 1);
            j = 0;
            for(i=0; i<m->num_vertices; ++i, j+=6)
            {
                m->vertices[i].x = obuff1[j];
                m->vertices[i].y = obuff1[j+1];
                m->vertices[i].z = obuff1[j+2];
                m->vnormals[i].x = obuff1[j+3];
                m->vnormals[i].y = obuff1[j+4];
                m->vnormals[i].z = obuff1[j+5];
            }
            free(buff);
            free(obuff1);
        }
    }
    else if(m->is_vcolors && m->is_vnormals)
    {
        if((m->vcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        if((m->vnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_vertices)))==NULL) mesh_error(MESH_ERR_MALLOC);
        if(m->origin_type==MESH_ORIGIN_TYPE_PLY_ASCII)
        {
            for(i=0; i<m->num_vertices; ++i)
            {
#if MESH_FLOATDATA_TYPE==0
                if(fscanf(fp, " %f %f %f %f %f %f %f %f %f %f \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vnormals[i].x, &m->vnormals[i].y, &m->vnormals[i].z, &m->vcolors[i].r, &m->vcolors[i].g, &m->vcolors[i].b, &m->vcolors[i].a)!=10)
#else
                if(fscanf(fp, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", &m->vertices[i].x, &m->vertices[i].y, &m->vertices[i].z, &m->vnormals[i].x, &m->vnormals[i].y, &m->vnormals[i].z, &m->vcolors[i].r, &m->vcolors[i].g, &m->vcolors[i].b, &m->vcolors[i].a)!=10)
#endif
                {
                    free(m->vertices);
                    free(m->vnormals);
                    free(m->vcolors);
                    m->is_vertices = 0;
                    m->is_vnormals = 0;
                    m->is_vcolors = 0;
                    return m;
                }
            }
        }
        else
        {
            INTDATA csz = m->num_vertices*(plydsizes[dtypes[__v]].t*3+plydsizes[dtypes[__vn]].t*3+plydsizes[dtypes[__vc]].t*4);
            if((buff = (char*)malloc(csz))==NULL) mesh_error(MESH_ERR_MALLOC);
            if(fread(buff, 1, csz, fp)!=(size_t)csz)
            {
                free(m->vertices);
                free(m->vnormals);
                free(m->vcolors);
                free(buff);
                m->is_vertices = 0;
                m->is_vnormals = 0;
                m->is_vcolors = 0;
                return m;
            }
            tmp[0] = dtypes[__v];
            tmp[1] = dtypes[__vn];
            tmp[2] = dtypes[__vc];
            tmp[3] = 3;
            tmp[4] = 3;
            tmp[5] = 4;
            if(dtypes[__vc]<6)
            {
                if((obuff1 = (FLOATDATA*)malloc(m->num_vertices*sizeof(FLOATDATA)*6))==NULL) mesh_error(MESH_ERR_MALLOC);
                if((obuff2 = (INTDATA*)malloc(m->num_vertices*sizeof(INTDATA)*4))==NULL) mesh_error(MESH_ERR_MALLOC);
                if(m->origin_type==dtypes[__endian]) __mesh_convert_format_mixed(tmp, tmp+3, 3, buff, obuff2, obuff1, m->num_vertices, 0);
                else __mesh_convert_format_mixed(tmp, tmp+3, 3, buff, obuff2, obuff1, m->num_vertices, 1);
                j = 0;
                j2 = 0;
                for(i=0; i<m->num_vertices; ++i, j+=6, j2+=4)
                {
                    m->vertices[i].x = obuff1[j];
                    m->vertices[i].y = obuff1[j+1];
                    m->vertices[i].z = obuff1[j+2];
                    m->vnormals[i].x = obuff1[j+3];
                    m->vnormals[i].y = obuff1[j+4];
                    m->vnormals[i].z = obuff1[j+5];
                    m->vcolors[i].r = obuff2[j2];
                    m->vcolors[i].g = obuff2[j2+1];
                    m->vcolors[i].b = obuff2[j2+2];
                    m->vcolors[i].a = obuff2[j2+3];
                }
                free(obuff1);
                free(obuff2);
            }
            else
            {
                if((obuff1 = (FLOATDATA*)malloc(m->num_vertices*sizeof(FLOATDATA)*10))==NULL) mesh_error(MESH_ERR_MALLOC);
                if(m->origin_type==dtypes[__endian]) __mesh_convert_format_mixed(tmp, tmp+3, 3, buff, NULL, obuff1, m->num_vertices, 0);
                else __mesh_convert_format_mixed(tmp, tmp+3, 3, buff, NULL, obuff1, m->num_vertices, 1);
                j = 0;
                for(i=0; i<m->num_vertices; ++i, j+=10)
                {
                    m->vertices[i].x = obuff1[j];
                    m->vertices[i].y = obuff1[j+1];
                    m->vertices[i].z = obuff1[j+2];
                    m->vnormals[i].x = obuff1[j+3];
                    m->vnormals[i].y = obuff1[j+4];
                    m->vnormals[i].z = obuff1[j+5];
                    m->vcolors[i].r = obuff1[j+6];
                    m->vcolors[i].g = obuff1[j+7];
                    m->vcolors[i].b = obuff1[j+8];
                    m->vcolors[i].a = obuff1[j+9];
                }
                free(obuff1);
            }
            free(buff);
        }
    }
    return m;
}

MESH __mesh_parse_ply_faces(MESH m, FILEPOINTER fp, int *dtypes)
{
    INTDATA i, j, nv;
    char dummy[32];
    long int currpos;
    int nwrds = 0, flag = 0;
    INTDATA nverts = 0;
    if(m->is_faces) free(m->faces);
    m->is_faces = 0;
    if(m->num_faces!=0)
    {
        if((m->faces = (MESH_FACE)malloc(sizeof(mesh_face)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
        m->is_faces = 1;
        m->is_trimesh = 1;
        if(m->origin_type==MESH_ORIGIN_TYPE_PLY_ASCII)
        {
            currpos = ftell(fp);
            mesh_count_words_in_line(fp, &nwrds);
            fseek(fp, currpos, SEEK_SET);
#if MESH_INTDATA_TYPE==0
            if(fscanf(fp, "%d", &nverts)!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#else
            if(fscanf(fp, "%ld", &nverts)!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#endif
            fseek(fp, currpos, SEEK_SET);
            if(nwrds==(nverts+1) && !m->is_fcolors && !m->is_fnormals)
            {
                m->is_faces = 1;
                for(i=0; i<m->num_faces; ++i)
                {
                    do
                    {
                        flag = mesh_read_word_only(fp, dummy, 32);
                    }
                    while(flag==3);
                    if(flag>0)
                    {
                        free(m->faces);
                        m->is_faces = 0;
                        return m;
                    }
                    nv = strtol(dummy, NULL, 0);
                    m->faces[i].num_vertices = nv;
                    if(nv!=3) m->is_trimesh = 0;
                    if((m->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*nv))==NULL) mesh_error(MESH_ERR_MALLOC);
                    for(j=0; j<nv; ++j)
                    {
#if MESH_INTDATA_TYPE==0
                        if(fscanf(fp, "%d", &(m->faces[i].vertices[j]))!=1)
#else
                        if(fscanf(fp, "%ld", &(m->faces[i].vertices[j]))!=1)
#endif
                        {
                            free(m->faces);
                            m->is_faces = 0;
                            return m;
                        }
                    }
                }
            }
            else if(nwrds==(nverts+2) && m->is_fcolors && !m->is_fnormals)
            {
                if(m->is_fcolors) free(m->fcolors);
                if((m->fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
                m->is_faces = 1;
                for(i=0; i<m->num_faces; ++i)
                {
                    do
                    {
                        flag = mesh_read_word_only(fp, dummy, 32);
                    }
                    while(flag==3);
                    if(flag>0)
                    {
                        free(m->faces);
                        m->is_faces = 0;
                        return m;
                    }
                    nv = strtol(dummy, NULL, 0);
                    m->faces[i].num_vertices = nv;
                    if(nv!=3) m->is_trimesh = 0;
                    if((m->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*nv))==NULL) mesh_error(MESH_ERR_MALLOC);
                    for(j=0; j<nv; ++j)
                    {
#if MESH_INTDATA_TYPE==0
                        if(fscanf(fp, "%d", &(m->faces[i].vertices[j]))!=1)
#else
                        if(fscanf(fp, "%ld", &(m->faces[i].vertices[j]))!=1)
#endif
                        {
                            mesh_error(MESH_ERR_SIZE_MISMATCH);
                        }
                    }
#if MESH_FLOATDATA_TYPE==0
                    if(fscanf(fp, "%f", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#else
                    if(fscanf(fp, "%lf", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#endif
                    m->fcolors[i].g = m->fcolors[i].r;
                    m->fcolors[i].b = m->fcolors[i].r;
                    m->fcolors[i].a = 1.0f;
                }
            }
            else if(nwrds==(nverts+4) && m->is_fcolors && !m->is_fnormals)
            {
                if(m->is_fcolors) free(m->fcolors);
                if((m->fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
                m->is_faces = 1;
                for(i=0; i<m->num_faces; ++i)
                {
                    do
                    {
                        flag = mesh_read_word_only(fp, dummy, 32);
                    }
                    while(flag==3);
                    if(flag>0)
                    {
                        free(m->faces);
                        m->is_faces = 0;
                        return m;
                    }
                    nv = strtol(dummy, NULL, 0);
                    m->faces[i].num_vertices = nv;
                    if(nv!=3) m->is_trimesh = 0;
                    if((m->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*nv))==NULL) mesh_error(MESH_ERR_MALLOC);
                    for(j=0; j<nv; ++j)
                    {
#if MESH_INTDATA_TYPE==0
                        if(fscanf(fp, "%d", &(m->faces[i].vertices[j]))!=1)
#else
                        if(fscanf(fp, "%ld", &(m->faces[i].vertices[j]))!=1)
#endif
                        {
                            mesh_error(MESH_ERR_SIZE_MISMATCH);
                        }
                    }
#if MESH_FLOATDATA_TYPE==0
                    if(fscanf(fp, "%f", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%f", &(m->fcolors[i].g))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%f", &(m->fcolors[i].b))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#else
                    if(fscanf(fp, "%lf", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%lf", &(m->fcolors[i].g))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%lf", &(m->fcolors[i].b))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#endif
                    m->fcolors[i].a = 1.0f;
                }
            }
            else if(nwrds==(nverts+5) && m->is_fcolors && !m->is_fnormals)
            {
                if(m->is_fcolors) free(m->fcolors);
                if((m->fcolors = (MESH_COLOR)malloc(sizeof(mesh_color)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
                m->is_faces = 1;
                for(i=0; i<m->num_faces; ++i)
                {
                    do
                    {
                        flag = mesh_read_word_only(fp, dummy, 32);
                    }
                    while(flag==3);
                    if(flag>0)
                    {
                        free(m->faces);
                        m->is_faces = 0;
                        return m;
                    }
                    nv = strtol(dummy, NULL, 0);
                    m->faces[i].num_vertices = nv;
                    if(nv!=3) m->is_trimesh = 0;
                    if((m->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*nv))==NULL) mesh_error(MESH_ERR_MALLOC);
                    for(j=0; j<nv; ++j)
                    {
#if MESH_INTDATA_TYPE==0
                        if(fscanf(fp, "%d", &(m->faces[i].vertices[j]))!=1)
#else
                        if(fscanf(fp, "%ld", &(m->faces[i].vertices[j]))!=1)
#endif
                        {
                            mesh_error(MESH_ERR_SIZE_MISMATCH);
                        }
                    }
#if MESH_FLOATDATA_TYPE==0
                    if(fscanf(fp, "%f", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%f", &(m->fcolors[i].g))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%f", &(m->fcolors[i].b))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%f", &(m->fcolors[i].a))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#else
                    if(fscanf(fp, "%lf", &(m->fcolors[i].r))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%lf", &(m->fcolors[i].g))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%lf", &(m->fcolors[i].b))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%lf", &(m->fcolors[i].a))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#endif
                }
            }
            else if(nwrds==(nverts+4) && !m->is_fcolors && m->is_fnormals)
            {
                if(m->is_fnormals) free(m->fnormals);
                if((m->fnormals = (MESH_NORMAL)malloc(sizeof(mesh_normal)*(m->num_faces)))==NULL) mesh_error(MESH_ERR_MALLOC);
                m->is_faces = 1;
                for(i=0; i<m->num_faces; ++i)
                {
                    do
                    {
                        flag = mesh_read_word_only(fp, dummy, 32);
                    }
                    while(flag==3);
                    if(flag>0)
                    {
                        free(m->faces);
                        m->is_faces = 0;
                        return m;
                    }
                    nv = strtol(dummy, NULL, 0);
                    m->faces[i].num_vertices = nv;
                    if(nv!=3) m->is_trimesh = 0;
                    if((m->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*nv))==NULL) mesh_error(MESH_ERR_MALLOC);
                    for(j=0; j<nv; ++j)
                    {
#if MESH_INTDATA_TYPE==0
                        if(fscanf(fp, "%d", &(m->faces[i].vertices[j]))!=1)
#else
                        if(fscanf(fp, "%ld", &(m->faces[i].vertices[j]))!=1)
#endif
                        {
                            mesh_error(MESH_ERR_SIZE_MISMATCH);
                        }
                    }
#if MESH_FLOATDATA_TYPE==0
                    if(fscanf(fp, "%f", &(m->fnormals[i].x))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%f", &(m->fnormals[i].y))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%f", &(m->fnormals[i].z))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#else
                    if(fscanf(fp, "%lf", &(m->fnormals[i].x))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%lf", &(m->fnormals[i].y))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
                    if(fscanf(fp, "%lf", &(m->fnormals[i].z))!=1) mesh_error(MESH_ERR_SIZE_MISMATCH);
#endif
                }
            }
        }
        else
        {
            char* buff;
            INTDATA csz1, csz2, tmp[4];
            csz1 =  plydsizes[dtypes[__ft]].t;

            if((buff = (char*)malloc(256))==NULL) mesh_error(MESH_ERR_MALLOC);
            for(i=0; i<m->num_faces; ++i)
            {
                if(fread(buff, 1, csz1, fp)!=(size_t)csz1) /* num of vertices */
                {
                    for(j=0; j<i; ++j) free(m->faces[j].vertices);
                    free(m->faces);
                    m->faces = NULL;
                    m->is_faces = 0;
                    free(buff);
                    return m;
                }
                tmp[0] = dtypes[__ft];
                tmp[1] = 1;
                if(m->origin_type==dtypes[__endian]) __mesh_convert_format_mixed(tmp, tmp+1, 1, buff, &nv, NULL, 1, 0);
                else __mesh_convert_format_mixed(tmp, tmp+1, 1, buff, &nv, NULL, 1, 1);
                csz2 = plydsizes[dtypes[__f]].t*nv;
                m->faces[i].num_vertices = nv;
                if(nv!=3) m->is_trimesh = 0;
                if((m->faces[i].vertices = (INTDATA *)malloc(sizeof(INTDATA)*nv))==NULL) mesh_error(MESH_ERR_MALLOC);
                if(fread(buff, 1, csz2, fp)!=(size_t)csz2) /* vertex indices */
                {
                    for(j=0; j<=i; ++j) free(m->faces[j].vertices);
                    free(m->faces);
                    m->faces = NULL;
                    m->is_faces = 0;
                    free(buff);
                    return m;
                }
                tmp[0] = dtypes[__f];
                tmp[1] = nv;
                if(m->origin_type==dtypes[__endian]) __mesh_convert_format_mixed(tmp, tmp+1, 1, buff, m->faces[i].vertices, NULL, 1, 0);
                else __mesh_convert_format_mixed(tmp, tmp+1, 1, buff, m->faces[i].vertices, NULL, 1, 1);
            }
            free(buff);
        }
    }
    return m;
}

/** \endcond */
