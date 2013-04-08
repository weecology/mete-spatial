/*
Purpose: to compute the separation order distance between 
each unique pairwise site comparison. This is carried out 
in conjunction with the R function dist_bisect().
*/ 

void get_sep_dist(int *vec_coords, int *i_bisect, int *N, int *sep_dist){
  int i, j, icount, diff, sep ;
  icount = 0 ;
  for (i = 1; i < *N; i++) {
    for (j = (i + 1); j < *N + 1; j++) {
      diff = 0 ;
      sep = 0 ;
      while (diff == 0) {
        sep++ ; 
        diff = vec_coords[(i - 1) * *i_bisect + sep - 1] -
               vec_coords[(j - 1) * *i_bisect + sep - 1] ; 
      }
      sep_dist[icount] = sep ; 
      icount++ ; 
    }
  }
}

