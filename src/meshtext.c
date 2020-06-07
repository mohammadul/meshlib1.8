/**
 * @file meshtext.c
 * @author Sk. Mohammadul Haque
 * @version 1.4.2.0
 * @copyright
 * Copyright (c) 2013, 2014, 2015, 2016 Sk. Mohammadul Haque.
 * @brief This file contains functions pertaining to different text routines.
 */

#include <string.h>
#include "../include/meshlib.h"

/** \brief Checks if numeric or not
 *
 * \param[in] fp Pointer to input file
 * \return 1 for numeric/ else - for non-numeric
 *
 */

int mesh_isnumeric(FILEPOINTER fp)
{
    char flag = 0;
    int ch;
    while((ch = getc(fp))!=EOF)
    {
        if(ch!=' ')
        {
            if((ch>=48) &&(ch<=57))
            {
                ungetc(ch, fp);
                if(flag==1) ungetc('.', fp);
                if(flag==2) ungetc('-', fp);
                if(flag==3)
                {
                    ungetc('.', fp);
                    ungetc('-', fp);
                }
                return 1;
            }
            else
            {
                if(ch=='.' )
                {
                    if(flag==2) flag = 3;
                    else flag = 1;
                }
                else if(ch=='-') flag = 2;

                else
                {
                    ungetc(ch, fp);
                    return 0;
                }
            }
        }
        else
        {
            if(flag>0) /* inserted condition for initial space in line */
            {
                ungetc(ch, fp);
                return 0;
            }
        }
    }
    if(ch==EOF)
    {
        ungetc(ch, fp); /*edited just now */

        return -1;
    }
    return 0;
}

/** \brief Points to the next word
 *
 * \param[in] fp Pointer to input file
 * \return Status 0 - Normal/ 1- EOF
 *
 */

int mesh_go_next_word(FILEPOINTER fp)
{
    int flag = 0, ch;
    while((flag<2) && ((ch = getc(fp))!=EOF))
    {
       if ((ch =='\v') || (ch =='\r') || (ch =='\n') || (ch =='\t') || (ch ==' ') || (ch =='\f') || (ch==',') || (ch=='!') || (ch=='(') || (ch==')') || (ch=='{') || (ch=='}') || (ch=='[') || (ch==']'))
{
            if(flag==0) flag = 1;
        }
        else if(flag==1) flag = 2;
    }
    if(ch!=EOF)
    {
        ungetc(ch, fp);
        return 0;
    }
    else return 1;
}

/** \brief Counts number of words in the current line
 *
 * \param[in] fp Pointer to input file
 * \param[out] count Count
 * \return Status 0 - Normal/ 1- EOF
 *
 */

int mesh_count_words_in_line(FILEPOINTER fp, int *count)
{
    int flag = -1, ch;
    *count = 0;
    while((flag<3) && ((ch = (char)getc(fp))!=EOF))
    {
        if((ch=='\v') || (ch=='\n'))
        {
            if (flag==0)
            {
                ++(*count);
                flag = 3;
            }
            if(flag==-1) flag = 4;/*  line included to handle empty line */
            else flag = 2;
        }
        if ((ch =='\v') || (ch =='\r') || (ch =='\n') || (ch =='\t') || (ch ==' ') || (ch =='\f') || (ch==',') || (ch=='!') || (ch=='(') || (ch==')') || (ch=='{') || (ch=='}') || (ch=='[') || (ch==']'))
    {
            if(flag==0)/* changed from  <=0 to ==0 to skip initial space */
            {
                flag = 1;
                ++(*count);
            }
        }
        else if(flag==1) flag = 0;
        else if(flag==2) flag = 3;
        else flag = 0;
    }
    if(flag!=-1 && flag!=4) ungetc(ch, fp);
    if(ch==EOF)
    {
        if(flag==0) ++(*count);
        return 1;
    }
    else return 0;
}

/** \brief Reads current word and moves to the next word
 *
 * \param[in] fp Pointer to input file
 * \param[out] c_word Variable to store the word
 * \param[in] sz Maximum size to read
 * \return Status 0 - Normal/ 1- EOF
 *
 */

int mesh_read_word(FILEPOINTER fp, char *c_word, int sz)
{
    int flag = 0, t = 0, comment_flag = 0;
    int ch;
    while((flag<3) && ((ch = getc(fp))!=EOF))/*no need for state 3 to be corrected*/
    {
        if ((ch =='\v') || (ch =='\r') || (ch =='\n') || (ch =='\t') || (ch ==' ') || (ch =='\f') || (ch==',') || (ch=='!') || (ch=='(') || (ch==')') || (ch=='{') || (ch=='}') || (ch=='[') || (ch==']'))
        {
            if(flag!=0) flag = 2;
        }
        else if(flag<2)
        {
            flag = 1;
            if(ch=='#')
            {
                flag = 4; /* off comment */
                comment_flag = 1;
                /* skip remaining part of the line */
                while((comment_flag==1) && ((ch = getc(fp))!=EOF))/*no need for state 3 to be corrected*/
                {
                    if(ch=='\n')
                    {
                        comment_flag = 0;
                    }
                }
                if(ch!=EOF)
                {
                    ungetc(ch, fp);
                    return 3;
                }
                else return 3;
            }
            c_word[t] = ch;
            if((t+1)>=sz)
            {
                c_word[t+1] = '\0'; /* prevent buffer overflow */
                ungetc(ch, fp);
                return 2;
            }
            else ++t;
        }
        else if(flag==2) flag = 3; /* reached next word */ /*to be corrected for deleting state 3*/
    }
    c_word[t] = '\0';
    if(ch!=EOF)
    {
        ungetc(ch, fp);
        return 0;
    }
    else return 1;
}

/** \brief Reads current word withot moving to the next word
 *
 * \param[in] fp Pointer to input file
 * \param[out] c_word Variable to store the word
 * \param[in] sz Maximum size to read
 * \return Status 0 - Normal/ 1- EOF
 *
 */

int mesh_read_word_only(FILEPOINTER fp, char *c_word, int sz)
{
    int flag = 0, t = 0, comment_flag = 0, ch;
    while((flag<2) && ((ch = getc(fp))!=EOF))/*no need for state 3 to be corrected*/
    {
        if ((ch =='\v') || (ch =='\r') || (ch =='\n') || (ch =='\t') || (ch ==' ') || (ch =='\f') || (ch==',') || (ch=='!') || (ch=='(') || (ch==')') || (ch=='{') || (ch=='}') || (ch=='[') || (ch==']'))
        {
            if(flag!=0) flag = 2;
        }
        else if(flag<2)
        {
            flag = 1;
            if(ch=='#')
            {
                flag = 4; /* off comment */
                comment_flag = 1;
                /* skip remaining part of the line */
                while((comment_flag==1) && ((ch = getc(fp))!=EOF))/*no need for state 3 to be corrected*/
                {
                    if(ch=='\n')
                    {
                        comment_flag = 0;
                    }
                }
                if(ch!=EOF)
                {
                    ungetc(ch, fp);
                    return 3;
                }
                else return 3;
            }
            c_word[t] = ch;
            if((t+1)>=sz)
            {
                c_word[t+1] = '\0'; /* prevent buffer overflow */
                ungetc(ch, fp);
                return 2;
            }
            else ++t;
        }
    }
    c_word[t] = '\0';
    if(ch!=EOF)
    {
        ungetc(ch, fp);
        return 0;
    }
    else return 1;
}

/** \brief Skips to next line
 *
 * \param[in] fp Pointer to input file
 * \return Status 0 - Normal/ 1- EOF
 *
 */

int mesh_skip_line(FILEPOINTER fp)
{
    char comment_flag = 1;
    int ch;
    /* skip remaining part of the line */
    while((comment_flag==1) && ((ch = getc(fp))!=EOF))/*no need for state 3 to be corrected*/
    {
        if (ch=='\r' || ch=='\n')
        {
            comment_flag = 0;
        }
    }
    if(ch!=EOF) return 0;
    else return 1;
}

