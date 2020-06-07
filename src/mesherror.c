/**
 * @file mesherror.c
 * @author Sk. Mohammadul Haque
 * @version 1.4.2.0
 * @copyright
 * Copyright (c) 2013, 2014, 2015, 2016 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to handling errors.
 */

#include "../include/meshlib.h"

/** \brief Displays error message and exits
 *
 * \param[in] type Error type (MESH_ERR_MALLOC/MESH_ERR_SIZE_MISMATCH/MESH_ERR_FNOTOPEN)
 * \return NULL
 *
 */

void mesh_error(int type)
{
    switch(type)
    {
    case MESH_ERR_MALLOC:
        fprintf(stderr, "mesh: malloc error\n");
        break;

    case MESH_ERR_SIZE_MISMATCH:
        fprintf(stderr, "mesh: dimension size mismatch\n");
        break;

    case MESH_ERR_FNOTOPEN:
        fprintf(stderr, "mesh: can not open file\n");
        break;

    case MESH_ERR_INCOMPATIBLE:
        fprintf(stderr, "mesh: incompatible data\n");
        break;

    default:
        fprintf(stderr, "mesh: unknown error\n");
        break;
    }
    exit(-1);
}

