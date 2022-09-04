  #include <math.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include "ran2.h"
  #include "myparam.h"
  #include "confinementcheck.h"
  #include "checklineintersect.h"
  #include "update_position.h"
  #include "initialize_polymer.h"

  // Dai's variable
  double knottable[Nknottypetotal][3];

  FILE *fp_knottype, *fp_knotsize;

  // Functions
  void Reflection_from_sphere(double* pt_in, double* next_pt, double* reflection);
  void update_v_vec(double* dir_vec, double* dir_buffer, double* v_vec, double* buffer, double* angle, double* a_buffer);
  int check_in_sphere(double x, double y, double z){
    if(sqrt(pow(x,2) + pow(y,2) + pow(z,2)) > radius){
      return 0;
    }
    return 1;
  }
  double Dot_product(double* vec1, double* vec2);
  int check_line_intersect(double* point1, double* point2, double* point3, double* point4);
  void initialize_polymer(double** x, double** v_vec, double** dir_vec, double** angle, long* dvar);
  void update_position(double* dx, double* dy, double* dz, double** x, double** v_vec, long* dvar);

  //Dai's functions
  int    get_knottype_of_entire_chain(double x[NBmax][dd], int Lchain); 
  int    get_knottype_of_subchain(double x[NBmax][dd], int ileft, int iright);
  void   get_knotcore(double x[NBmax][dd], int Lchain, int knottype, int result[3]) ;
  void   get_knotcore_circular(double x[NBmax][dd], int Lchain, int knottype_init, int knotsize[3] );
  void   get_knotcore_faster(double x[NBmax][dd], int Lchain, int knottype, int result[3]);


  int  main(){
        double x[TS][NB][dd];
        double v_vec[TS][NB][dd], dir_vec[TS][NB][dd];
        double angle[TS][NB][2]; // angle[0] = phi ; angle[1] = theta 
        double theta, phi;
        double dtheta[NB], dphi{NB};
        double dx[NB], dy[NB], dz[NB];
        double pt_in[3];
        double next_pt[NB][3];
        double speed_dir[3];
        double intersection_pt[3];
        double v_buffer[NB][3], dir_buffer[NB][3], a_buffer[NB][3];
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
        /*
        for(irep =0; irep< 1307; irep++){
          ran2(dvar,'u');
          ran2(dvar,'u');
          for(i=1; i <NB;i++){
            ran2(dvar, 'g');ran2(dvar, 'g');ran2(dvar, 'g');
            ran2(dvar, 'g');ran2(dvar, 'g');ran2(dvar, 'g');
          }
        }*/

        // MAIN CODE: TEST RUNS FOR irep TIMES: WALK FOR NB STEPS
        for(irep=0; irep < Nstep; irep++) {
          // Create the file for recording trajectories
          fprintf(walk_traj,"%d\n%d's trial\n",NB,irep);
          num_bounce = 0;

          // Initializing posiiton, angular direction, and velocity 
          initialize_polymer(x[0], v_vec[0], dir_vec[0], angle[0], dvar);
        
          // NOW WALK
          for(walk_step=1;walk_step<NB;walk_step++){
            // Translational diffusion
            update_position(dx, dy, dz, x[walk_step-1],v_vec[walk_step-1], dvar);

            // Save the current velocity and direction; if not reflected, next point's velocity is v_buffer + diffusion
            buffer_vecs(dir_vec[walk_step-1], v_vec[walk_step-1], angle[walk_step-1], dir_buffer, v_buffer, a_buffer);

            // Particle moves to next_pt
            disp_diffuse(next_pt, x[walk_step-1], dx, dy, dz);
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
            adopt_new_pos(x, next_pt);

            // Speed angle now is subject to rotational diffusion: Strategy is first to bounce back then diffuse
            update_v_vec(dir_vec[walk_step], dir_buffer, v_vec[walk_step], v_buffer, angle[walk_step], a_buffer);
  
            // Record the positions into the trajectory file
            fprintf(walk_traj,"%d\t%10lf\t%10lf\t%10lf\n",irep,x[walk_step-1][0],x[walk_step-1][1],x[walk_step-1][2]);
        }

        // Check if trajectories intersect with each other: No one wants a piece of knot to be passing directly on each other
        for(i=0;i<NB;i++){
          for(j=i+1;j<NB;j++){
            intersect_bool = check_line_intersect(x[i],x[i+1],x[j],x[j+1]);
          
            if(intersect_bool==1){
              printf("ERROR: line-line interesection at i = %d and j = %d\n\n", i,j);
              exit(-1);
            }
          }
        }

        // Save the last position
        fprintf(walk_traj,"%d\t%10lf\t%10lf\t%10lf\n",irep,x[NB-1][0],x[NB-1][1],x[NB-1][2]);  
        fprintf(walk_traj,"bounced for %d times\n\n",num_bounce);

        // Begin identify the knot
        knottype = get_knottype_of_entire_chain(x,NB);
    
        Ncross = 99; // unknown knotype
        for(i=0;i<Nknottypetotal;i++) {
          if(fabs( knottype - knottable[i][1]) < small_value) {
              Ncross = (int)(knottable[i][0] + 0.1);
              break;
          }
        }
       

        fprintf(fp_knottype,"%7d %6d %3d %3d\n",irep,knottype,Ncross, num_bounce);

        if((knottype != 999999) && (knottype !=0)){
          get_knotcore_circular(x, NB, knottype, knotsize);
          fprintf(fp_knotsize,"%7d  %6d  %5d %5d %5d %3d\n",irep, knottype, knotsize[0], knotsize[1], knotsize[2], num_bounce);
        }
      }
  
      fclose(walk_traj);
      fclose(fp_knotsize);
      fclose(fp_knottype);

      
      return 0;
  }

