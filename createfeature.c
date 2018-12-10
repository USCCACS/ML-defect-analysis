#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atom_property.h"
#include "createfeature.h"

int cal_feature(a_coodrinates *atoms,int centeratom,int **nlist,float *boxmd,float *halfboxmd){
    int iatom,jatom,katom,jj=0,count=0;
    int first_nn[max_neighbour];
    float distance[4],first_nn_avg=0.0,first_nn_min=10000.0,first_nn_max=0.0;
    float nn_nn_avg=0.0,nn_nn_min=10000.0,nn_nn_max=0.0;
    float sec_nn=0.0,sec_nn_avg=0.0,sec_nn_min=10000.0,sec_nn_max=0.0;
    float navg_34=0,navg_45=0,navg_56=0,ravg_34=0,ravg_45=0,ravg_56=0;
    float *features;

    features=(float *)malloc(N_feature*sizeof(float));
    memset(features,0.0,N_feature*sizeof(float));
    memset(first_nn,0,max_neighbour*sizeof(int));

    for (int ii=1;ii<=nlist[centeratom][0];ii++){
        jatom=nlist[centeratom][ii];
        cal_distance(centeratom,jatom,atoms,boxmd,halfboxmd,distance);
        if(distance[3] <= 3.0) {
          first_nn[0]+=1;
          first_nn[first_nn[0]]=jatom;
          first_nn_avg+=distance[3];
          if(first_nn_min > distance[3]) first_nn_min = distance[3];
          if(first_nn_max < distance[3]) first_nn_max = distance[3];
        }
        else if (distance[3] > 3.0 && distance[3] <= 4.0){
                navg_34 +=1.0;
                ravg_34 +=distance[3];
        }else if (distance[3] > 4.0 && distance[3] <= 5.0){
                navg_45 +=1.0;
                ravg_45 +=distance[3];
        }else if (distance[3] > 5.0 && distance[3] <= 6.0){
                navg_56 +=1.0;
                ravg_56 +=distance[3];
        }
    }
    first_nn_avg/=(float)first_nn[0];
//    printf("%d %f %f %f\n",first_nn[0],first_nn_avg,first_nn_min,first_nn_max);
    features[jj++]=first_nn[0];
    features[jj++]=first_nn_avg;
    features[jj++]=first_nn_min;
    features[jj++]=first_nn_max;

    //compute statistics among neighbor atoms
    for (int ii=1; ii<first_nn[0]; ii++){
        for (int kk=ii+1; kk<=first_nn[0]; kk++){
            count++;
            iatom=first_nn[ii];
            katom=first_nn[kk];
            cal_distance(iatom,katom,atoms,boxmd,halfboxmd,distance);
            nn_nn_avg+=distance[3];
            if(nn_nn_min > distance[3]) nn_nn_min = distance[3];
            if(nn_nn_max < distance[3]) nn_nn_max = distance[3];
        }
    }
    nn_nn_avg/=(float)(count);
//    printf("%d %f %f %f\n",count,nn_nn_avg,nn_nn_min,nn_nn_max);
    features[jj++]=nn_nn_avg;
    features[jj++]=nn_nn_min;
    features[jj++]=nn_nn_max;

    // compute statistics of the nearest neighors
    for (int ii=1; ii<first_nn[0]; ii++){
        iatom=first_nn[ii];
        for (int kk=1; kk<=nlist[iatom][0];kk++){
            katom=nlist[iatom][kk];
            cal_distance(iatom,katom,atoms,boxmd,halfboxmd,distance);
            if (distance[3] <= 3.0) {
                sec_nn += 1;
                sec_nn_avg += distance[3];
                if(sec_nn_min > distance[3]) sec_nn_min = distance[3];
                if(sec_nn_max < distance[3]) sec_nn_max = distance[3];
            }
        }
    }
    sec_nn_avg /= sec_nn;
    sec_nn /= first_nn[0];
//    printf("%f %f %f %f\n",sec_nn,sec_nn_avg,sec_nn_min,sec_nn_max);
    features[jj++]=sec_nn;
    features[jj++]=sec_nn_avg;
    features[jj++]=sec_nn_min;
    features[jj++]=sec_nn_max;
    // higher order features
    features[jj++]= navg_34;
    if(navg_34 == 0) navg_34 = 1.0;
    features[jj++]= ravg_34/navg_34;
    features[jj++]= navg_45;
    if(navg_45 == 0) navg_45 = 1.0;
    features[jj++]= ravg_45/navg_45;
    features[jj++]= navg_56;
    if(navg_56 == 0) navg_56 = 1.0;
    features[jj++]= ravg_56/navg_56;
//    exit(0);
    atoms[centeratom].feature=features;
    return 0;
}

int cal_distance(int iatom, int jatom,a_coodrinates *atoms,float *boxmd,float *halfboxmd,float *retval){
    float rsq=0.0;
    float dr[3]={0.0,0.0,0.0};

    for(int kk=0;kk<3;kk++){
        dr[kk] = atoms[iatom].loc[kk] - atoms[jatom].loc[kk];
        if (dr[kk] > halfboxmd[kk]) {dr[kk] -= boxmd[kk];}
        if (dr[kk] < -halfboxmd[kk]) {dr[kk] += boxmd[kk];}
        rsq+=dr[kk]*dr[kk];
        retval[kk]=dr[kk];
    }
    retval[3]=sqrt(rsq);
    return(0);
}

int write_features(a_coodrinates *input_atoms,char *filename,int Natoms,float *box){
    FILE *fp;
    int total = N_feature;
  //  int bluk=0,surface=0,defect=0;
    fp=fopen(filename,"w");
    fprintf(fp,"%d \n",Natoms);
    fprintf(fp,"%d \n",total);
    fprintf(fp,"%12.6f %12.6f %12.6f \n",box[0],box[1],box[2]);

    for (int i=0; i<Natoms; i++){
        fprintf(fp,"%s \t %12.6f \t %12.6f \t %12.6f \t %10d \t",input_atoms[i].aname,
                        input_atoms[i].loc[0],input_atoms[i].loc[1],input_atoms[i].loc[2],input_atoms[i].label);
        for (int j=0;j<total;j++){
            fprintf(fp,"%12.6f \t",input_atoms[i].feature[j]);
        }
        fprintf(fp,"\n");
    }
    return(0);
}

int main(int argc, char* argv[]){
    char *filename,*outputfile;
    a_systeminfo mdatom_info;
    a_coodrinates *input_atoms;

    filename=argv[1];
    outputfile = argv[2];
    read_input(filename,&mdatom_info,&input_atoms);
    makelinkedlist(input_atoms,mdatom_info.Natoms,mdatom_info.cellsize,mdatom_info.ng,mdatom_info.llst,mdatom_info.lshd);
    makeneighboutlist(input_atoms,mdatom_info.Natoms,mdatom_info.boxmd,mdatom_info.halfboxmd,mdatom_info.cellsize,
                    mdatom_info.ng,mdatom_info.rcutoffsq,mdatom_info.llst,mdatom_info.lshd,&(mdatom_info.neigh_info));
    fprintf(stdout,"Total number of atoms %10d \n",mdatom_info.Natoms);
    fprintf(stdout,"Box size %12.6f %12.6f %12.6f \n",mdatom_info.boxmd[0],mdatom_info.boxmd[1],mdatom_info.boxmd[2]);

    //compute feature vector

    for(int i=0;i<mdatom_info.Natoms;i++){
        cal_feature(input_atoms,i,mdatom_info.neigh_info,mdatom_info.boxmd,mdatom_info.halfboxmd);
    }
    write_features(input_atoms,outputfile,mdatom_info.Natoms,mdatom_info.boxmd);
    return(0);
}
