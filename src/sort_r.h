/* Header for sort_r.c */

#ifndef SORT_R_H_
#define SORT_R_H_

#include <stdlib.h>

 // structs
typedef struct sort_r_data{
  void *arg;
  int (*compar)(const void *a1, const void *a2, void *aarg);
} sort_r_data;




// functions
int sort_r_arg_swap(void *s, const void *aa, const void *bb);
void sort_r(void *base, size_t nel, size_t width, int (*compar)(const void *a1, const void *a2, void *aarg), void *arg);





#endif // SORT_R_H_