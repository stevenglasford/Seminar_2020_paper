//----------------------------------------
//- Biotrophic fungi differential game, dimension 2
//----------------------------------------
#include "data_default.h"	// do not modify for basic examples
//----------------------------------------
const char    NAME[]            = "data_user_biotroph_fungi_game.h, April 2017 (I. Yegorov)";
const int     DIM               = 2;   //- space dimension
//- 0: local Hnum function ; 1/2: Hnum defined with dynamics and 
//     distributed_cost 
const int     COMMANDS          = 2;   functions
//const int     OPTIM             = MAXMIN;
//----------------------------
//- DF method parameters:
//----------------------------
const int     OPTIM             = MAXMIN;
const int     METHOD            = MFD;
const int     TYPE_SCHEME       = ENO2;
const int     TYPE_RK           = RK2;
//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 1000000;
      double                   T= 60.0;         //- Terminal time
//----------------------------
//- discretization parameters:
//----------------------------
//- if DT=0 the program will compute a DT based on space steps.
const double  DT                = 5e-3;
      double  CFL               = 0.5;
const int     NN=750;
      int     ND[DIM]           = {NN, NN};
      double  XMIN[DIM]         = {0.0, 0.0};
      //- bound of the domain
      double  XMAX[DIM]         = {1.5, 1.5};
//- periodic mesh (1:periodic, 0:otherwise)
const int     PERIODIC[DIM]     = { 0, 0 };
//- 0 : xi = center of cell ; 1 : xi contains the boundary
const int     MESH              = 1;
const int     cDIM              = 1;
const int     NCD[cDIM]         = { 2 };
const double  UMIN[cDIM]        = { 0.};
const double  UMAX[cDIM]        = { 1.};
const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 2 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};
//-----------------------------
//- BOUNDARY : 0 (Void, for FD only) or 1 (Dirichlet, using g_border,
//             for SL/FD) or 2 (Vxx=0, for SL/FD)
//-----------------------------
const int BOUNDARY=0;
const double  VBORD             = 0.0;
//-  Dirichlet boundary condition:
double g_border(double t, const double* x){
  return VBORD;
}

//- Mixed Neumann bc :  vx=g(t,x,v) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}

//-----------------------------
//- mainloop parameters
//-----------------------------
const int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 0;
const int     COMPUTE_TOPT      = 0;
const int     TOPT_TYPE         = 0;       //- 0= min time , 1= exit time

//- 1 to save "VFn.dat" every  SAVE_VFALL_STEP iterations
const int     SAVE_VF_ALL       = 1;       
const int     SAVE_VF_ALL_STEP  = 50;
const int     SAVE_VF_FINAL     = 1;

//- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR       = 0;       
const int     CHECK_ERROR_STEP  = 100;

//-------------------
//- "coupe" : used to make some cut into some plane ==> results in files
//            "coupe.dat"/"coupeex.dat" (if COUPE_DIM not equal to {0,...,0}) 
//-------------------
//- put 1 to keep the dimensions in the final savings.
const int     COUPE_DIMS[DIM] = {1 ,1 }; 

//- if COUPE_DIM[d]=0 then make precise the value of the cut in direction x_d.
const double  COUPE_VALS[DIM] = {0.,1.}; 

//---------------
//- initial data
//---------------

inline double   v0(const double * x) {

	return 0.0;

}

//------------------------------------------------------------------------------------------------------
//- dynamics and distributed cost functions used in the case of COMMANDS=0 or 1 (assumes only 1 control)
//------------------------------------------------------------------------------------------------------
inline void dynamics(const double* x, C u, double t, double* res){ return; }
inline double distributed_cost(const double* arg, C u, double t){ return 0.; }


const double  new_k = 1.0 / 6.0,
			  new_alpha = 0.2,
			  new_beta = 0.1,
			  gamma_rate = 0.06,
			  mu_rate = 0.03;
	
const double  n_1 = 9.0,
			  n_2 = 1.0;


//---------------------------------------------------------------
//- if COMMANDS=2 (2 player games): dynamics and distributed cost
//---------------------------------------------------------------

inline void dynamics2(const double* x, C u, C u2, double t, double* res) {
     double  new_M_1 = x[0],
             new_M_2 = x[1],
             u_1 = u[0],
             u_2 = u2[0];
     if (new_M_1 < 0.0)
         new_M_1 = 0.0;
     if (new_M_2 < 0.0)
         new_M_2 = 0.0;
     double  compet_term = 1.0 / (1.0 + new_beta * (n_1 * new_M_1...
                                  + n_2 * new_M_2));
     res[0] = new_M_1 * ((1.0 - u_1) * compet_term * new_alpha / ...
                         (new_M_1 + new_k) - gamma_rate);
     res[1] = new_M_2 * ((1.0 - u_2) * compet_term * new_alpha / ...
                         (new_M_2 + new_k) - gamma_rate);
}

inline double distributed_cost2(const double* x, C u, C u2, double t) {
     double new_M_1 = x[0],
            new_M_2 = x[1],
            u_1 = u[0],
            u_2 = u2[0];
     if (new_M_1 < 0.0)
            new_M_1 = 0.0;
     if (new_M_2 < 0.0)
            new_M_2 = 0.0;
     double compet_term = 1.0 / (1.0 + new_beta * (n_1 * new_M_1 + n_2 ...
                                                   * new_M_2));
     return compet_term * new_alpha * exp(-mu_rate * (T - t)) *
            (u_2 * new_M_2  / (new_M_2 + new_k) - u_1 ...
             * new_M_1 / (new_M_1 + new_k));
}


//---------------------------
//- Exact solution (if known)
//---------------------------
inline double Vex(double t, const double* x){ return 0.; }

//-----------------------------------------------
//- PARAMETERS : compute trajectory 
//-----------------------------------------------
// method of trajectory reconstruction and starting time
int           TRAJ_METHOD       = 0; 
double        time_TRAJ_START   = 0.00;

//const int     TRAJPT     = 0;  
//const double  initialpoint[TRAJPT*DIM] = {}; //- Zero initial point
const int     TRAJPT     = 1;  
const double  initialpoint[TRAJPT*DIM] = {-1.,0.0}; //- 1 initial point 

// stopping criteria
//- to stop traj reconstruction when val(x) <=min (val=topt here)
double        min_TRAJ_STOP     = 0.00;    
//- to stop traj reconstruction when val(x) >=max (val=topt here)
double        max_TRAJ_STOP     = 100.00;  
//- to stop traj reconstruction when g_target(x)<=0
int           TARGET_STOP       = 1;       
inline double g_target(double t, const double* x){ return v0(x);}  

// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
inline void u2_adverse(double t, const double* x, C& u){}


//-----------------------------------------------
//- MORE ADVANCED PARAMETERS
//-----------------------------------------------

//- Case COMMANDS=0: numerical Hamiltonian function Hnum.
inline void compute_Hconst(double* aMAX, double t){};
inline double Hnum(const double t, const double* x, const double vi, ...
		   const double* dv){ return 0.; };

//--------------------------------
//- PARAMETERS FOR ADVANCED USERS
//--------------------------------
//- if 1 then starts computation from data in VF.dat and not from the 
//  v0 function
const int     EXTERNALV0        = 0;           
//- to precompute mesh coordinates (faster but needs more memory).
const int     PRECOMPUTE_COORDS = 1;           

//--------------------
//- OBSTACLE g
//--------------------
const int OBSTACLE=0;
inline double g_obstacle(double t, const double* x){ return 0.; }

//--------------------------------
//- border: number of ghost cells in each direction (default is 
//  {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {2,2};

//--------------------
//- Special parameter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL
