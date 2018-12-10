#define N_feature 17
/* Features for each atom consists of
 * 1 total number of nearest neighbor
 * 2 average distance of nearest neighbor
 * 3 minimum distance of nearest neighbor
 * 4 maximum distance of nearest neighbor
 * 5 average distance between nearest neighbor
 * 6 minimum distance between nearest neighbor
 * 7 maximum distance between nearest neighbor
 * 8 average  number of neighobrs of nearest neighbor
 * 9 average distance of neighobrs of nearest neighbor
 * 10 minimum distance of neighobrs of nearest neighbor
 * 11 maximum distance of neighobrs of nearest neighbor
 * 12-13 number and distance  of neighbor between 3 to 4
 * 14-15 number and distance of neighbor between 4 to 5
 * 16-17 number and distance of neighbor between 5 to 6
 * common neighbors
 */
extern int cal_feature(a_coodrinates *atoms,int centeratom,int **nlist,float *boxmd,float *halfboxmd);
extern int cal_distance(int iatom, int jatomp,a_coodrinates *atoms,float *boxmd,float *halfboxmd,float *retval);
extern int write_features(a_coodrinates *input_atoms,char *filename,int Natoms,float *box);
