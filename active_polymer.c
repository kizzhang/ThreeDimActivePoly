  #include <math.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include "ran2.h"
  #include "myparam.h"
  #include "confinementcheck.h"
  #include "update_position.h"
  #include "initialize_polymer.h"
  #include "update_v_vec.h"
  #include "adopt_new_pos.h"
  #include "buffer_vecs.h"
  #include"initialize_polymer.h"

  // Dai's variable
  double knottable[Nknottypetotal][3];

  FILE *fp_knottype, *fp_knotsize;

  // Functions
  void Reflection_from_sphere(double* pt_in, double* next_pt, double* reflection);
  int check_in_sphere(double x, double y, double z){
    if(sqrt(pow(x,2) + pow(y,2) + pow(z,2)) > radius){
      return 0;
    }
    return 1;
  }
  double Dot_product(double* vec1, double* vec2);
  void initialize_polymer(double x[NB][3],  double dir_vec[NB][3], double angle[NB][2]);
  void update_position(double* dx, double* dy, double* dz, double x[NB][3], double dir_vec[NB][3], long randx_mat[NB][dd]);
  void buffer_vecs(double dir_vec[NB][3], double angle[NB][2], double dir_buffer[NB][3], double a_buffer[NB][2]);
  void disp_diffuse(double next_pt[NB][3], double x[NB][3], double* dx, double* dy, double* dz);
  void update_v_vec(double dir_vec[NB][3], double dir_buffer[NB][3], double angle[NB][2], double abuffer[NB][2], long randv_mat[NB][2]);;
  void adopt_new_pos(double x[NB][3], double next_pt[NB][3]);
  void populate_randx_mat(long randx_mat[NB][dd]);
  void populate_randv_mat(long randv_mat[NB][2]);
  //Dai's functions
  int    get_knottype_of_entire_chain(double x[NBmax][dd], int Lchain); 
  int    get_knottype_of_subchain(double x[NBmax][dd], int ileft, int iright);
  void   get_knotcore(double x[NBmax][dd], int Lchain, int knottype, int result[3]) ;
  void   get_knotcore_circular(double x[NBmax][dd], int Lchain, int knottype_init, int knotsize[3]);
  void   get_knotcore_faster(double x[NBmax][dd], int Lchain, int knottype, int result[3]);

  static double x[1][NB][dd];
  static double dir_vec[1][NB][dd];
  static double angle[1][NB][2]; // angle[0] = phi ; angle[1] = theta 
  double theta, phi;
  static double dtheta[NB], dphi[NB];
  static double dx[NB], dy[NB], dz[NB];
  static double next_pt[NB][3];
  static double dir_buffer[NB][3], a_buffer[NB][2];

  int  main(){
        FILE *walk_traj;
        char* filename[80];
        int intersect_bool;
        int num_bounce;
        int walk_step;
        FILE *fp_knotsize;
        FILE *fp_knottype;

        // Dai's variable
        int knottype, knotsize[3], Ncross;
        double small_value = 0.001;
        long dum = -100;
        long *dvar;
        int i,j;
        int irep;
        long randx_mat[NB][dd], randv_mat[NB][2]; // array that records random values for particles
        
        FILE *fp_knottable;

        dvar=&dum;
        
        // Open knot table
        if((fp_knottable = fopen("table_knot_Apoly_all.txt","r")) == NULL ) { fprintf(stderr,"cannot find table_knot_Apoly_all.txt\n"); exit(-1); }
        for(i=0;i<Nknottypetotal;i++) fscanf(fp_knottable,"%lf %lf %lf\n",&knottable[i][0],&knottable[i][1],&knottable[i][2]);
        fclose(fp_knottable);

        // Create file that records knot 
        fp_knotsize  = fopen("result_knotsize.txt","w");
        fp_knottype  = fopen("result_knottype.txt","w");

        // Create file that records trajectory
        walk_traj  = fopen("walk_trajectory.txt","w");

        // MAIN CODE: TEST RUNS FOR irep TIMES: WALK FOR NB STEPS
        for(irep=0; irep < Nstep; irep++){
          // Create the file for recording trajectories
          fprintf(walk_traj,"%d\n%d's trial\n",NB,irep);
          num_bounce = 0;

          // Create matrices for recording random values
          populate_randx_mat(randx_mat);
          populate_randv_mat(randv_mat);
          
          // Initializing position, angular direction, and velocity 
          initialize_polymer(x[0], dir_vec[0], angle[0]);
          
          /*fprintf(walk_traj,"0\n");
          for(i=0; i<NB; i++){
            fprintf(walk_traj,"%d\t%10lf\t%10lf\t%10lf\n",i,x[0][i][0],x[0][i][1],x[0][i][2]); 
          }*/

          // NOW WALK
          for(walk_step=1;walk_step<TS;walk_step++){
            printf("\n\nTHIS IS TIME STEP AT %d", walk_step);
            // Translational diffusion
            update_position(dx, dy, dz, x[0], dir_vec[0], randx_mat);

            // Save the current velocity and direction; if not reflected, next point's velocity is v_buffer + diffusion
            buffer_vecs(dir_vec[0], angle[0], dir_buffer, a_buffer);

            // Particle moves to next_pt
            disp_diffuse(next_pt, x[0], dx, dy, dz);
          
            /*
            // Check if the particle hits the confinement boundary
            if(check_in_sphere(next_pt[0], next_pt[1], next_pt[2]) == 0){
              Reflection_from_sphere(x[walk_step-1], next_pt, speed_dir);   
              ++num_bounce;                                              // return the direction of reflected velocity vector,
              
              dir_buffer[0] = speed_dir[0]; dir_buffer[1] = speed_dir[1]; dir_buffer[2] = speed_dir[2]; 
              v_buffer[0] = walk_speed * speed_dir[0]; v_buffer[1] = walk_speed * speed_dir[1]; v_buffer[2] = walk_speed * speed_dir[2];  // and transform into reflected velocity
              a_buffer[0] = atan(speed_dir[1]/speed_dir[0]); a_buffer[1] = acos(speed_dir[2]); // angles are changed toooooo
            }*/ 

            // Next point specified
            adopt_new_pos(x[0], next_pt);
            
            
            // Speed angle now is subject to rotational diffusion: Strategy is first to bounce back then diffuse
            update_v_vec(dir_vec[0], dir_buffer, angle[0], a_buffer, randv_mat);
           
            // Record the positions into the trajectory file
            
            if(walk_step>=0){
              fprintf(walk_traj,"%d\n",walk_step);
              for(i=0;i<NB;i++){
                fprintf(walk_traj,"%d\t%10lf\t%10lf\t%10lf\n",i,x[0][i][0],x[0][i][1],x[0][i][2]);
              }
            }
          }
        //printf("Walkstep is %lf", walk_step);
        // Save the last position
        fprintf(walk_traj,"bounced for %d times\n\n",num_bounce);
        //printf("FINISHED!\n");
        // Begin identify the knot of the final configuration
      /*  knottype = get_knottype_of_entire_chain(x[TS-1],NB);
        printf("knottype = %lf\n", knottype);
        Ncross = 99; // unknown knotype
        for(i=0;i<Nknottypetotal;i++) {
          if(fabs( knottype - knottable[i][1]) < small_value) {
              Ncross = (int)(knottable[i][0] + 0.1);
              break;
          }
        }
        //printf("FINISHED!\n");
        fprintf(fp_knottype,"%7d %6d %3d %3d\n",irep,knottype,Ncross, num_bounce);

        if((knottype != 999999) && (knottype !=0)){
          get_knotcore_circular(x[TS-1], NB, knottype, knotsize);
          fprintf(fp_knotsize,"%7d  %6d  %5d %5d %5d %3d\n",irep, knottype, knotsize[0], knotsize[1], knotsize[2], num_bounce);
        }*/
      }
      //printf("HHH");
      fclose(walk_traj);
      fclose(fp_knotsize);
      fclose(fp_knottype);

      
      return 0;
  }

