#include <stdio.h>

void msg(const char *msg)
{
    fprintf(stdout, "\nMESSAGE:\n%s\n", msg);
}
void error_msg(const char *msg, const char *file, int line)
{
    fprintf(stderr, "\a\nERROR in file %s line %i:\n%s\n", file, line, msg);
    getchar();
}

void print_mat(const double *M, int row, int col, const char *c)
{
    int i, j;
    printf("\n%s:\n", c);
    for(i = 0; i < row; i++)
    {
        for(j = 0; j < col; j++) printf("%12.4g ", M[i * col + j]);
        printf("\n");
    }
}
