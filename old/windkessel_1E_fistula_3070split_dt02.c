#include "udf.h"
#include "unsteady.h"

//define IDOUT1 15

int p_var;								/* Initialize profile macro variable */
int no_face;							/* Initialize no. of faces integer variable*/
int thread_out_artery_ID = 6;			/* Initialize current time step artery flow rate */
int thread_out_vein_ID = 8;				/* Initialize current time step artery flow rate */

real den = 1060;						/* Blood density */
real dt = 0.02;							/* time step size */
real NV_VEC(A);							/* Initialize vector A for area calculations */
real Area = 0.;							/* Initialize Area for area calculations */
real t_np1 = -1;						/* Halt time stepping unitl pre-transient macro used */

/* Artery outlet parameters and variables*/
int node_id_artery;						/* Initialize node id for artery */

real Rp_Ar = 2.231840472085634e+09;		/* Artery Proximal Resistance (Near)*/

real Q_n_Ar = 0.;						/* Initialize current time step artery flow rate */
real Q_np1_Ar = 0.;						/* Initialize previous time step artery flow rate */
real P_n_Ar = 6.125e+03;				/* Initialize current time step artery flow rate */
real P_np1_Ar = 6.125e+03;				/* Initialize previous time step artery flow rate */

/* Vein outlet parameters and variables*/
int node_id_vein;						/* Initialize node id for vein */

real Rp_Ve = 9.206341947353240e+08;		/* Vein Proximal Resistance (Near)*/

real Q_n_Ve = 0.;						/* Initialize current time step vein flow rate */
real Q_np1_Ve = 0.;						/* Initialize previous time step vein flow rate */
real P_n_Ve = 6.025e+03;				/* Initialize current time step vein flow rate */
real P_np1_Ve = 6.025e+03;				/* Initialize previous time step vein flow rate */

 /***********************************************************************
	Macro that finds which threads are on which nodes and checking the area of the artery and vein outlets
 ************************************************************************/

DEFINE_ON_DEMAND(pre_initialize)	
{
	Domain *d_e = Get_Domain(1);
	Thread *f_thread_artery = Lookup_Thread(d_e,thread_out_artery_ID);
	Thread *f_thread_vein = Lookup_Thread(d_e,thread_out_vein_ID);

	face_t f;
	node_id_artery = 111;																	/* to ensure artery node id dosn't auto to 0 */
	node_id_vein = 111;																		/* to ensure vein node id dosn't auto to 0 */

	if (myid == 0)
	{
		Message("WINDKESSEL 1D UDF\n");
	}

	/* Find on which node the artery outlet is and calculate it's total area */
	no_face=0;
	Area = 0;

    begin_f_loop(f,f_thread_artery)
	{
		no_face++;
		F_AREA(A,f,f_thread_artery);
		Area += NV_MAG(A);
	}
    end_f_loop(f,f_thread_artery)

	if (no_face > 0)
	{
		node_id_artery = myid;
		Message("ARTERY \n - Initial Pressure (%d)Pa \n", P_n_Ar);
		Message(" - Outlet is on node (%d)\n", myid);
		Message(" - Outlet is made up of (%e) faces with area (%8.8e)\n", no_face , Area);
	}

	/* Find on which node the vein outlet is and calculate it's total area */
	no_face = 0;
	Area = 0;

    begin_f_loop(f,f_thread_vein)
	{
		no_face++;
		F_AREA(A,f,f_thread_vein);
		Area += NV_MAG(A);
	}
    end_f_loop(f,f_thread_vein)

	if (no_face > 0)
	{
		node_id_vein = myid;
		Message("VEIN \n - Initial Pressure (%d)Pa \n", P_n_Ve);
		Message(" - Outlet is on node (%d)\n", myid);
		Message(" - Outlet is made up of (%e) faces with area (%8.8e)\n", no_face , Area);
	}
} 

 /********************************************************************
   Macro that calculates the initialized flow rates of the artery and vein outlets and sets the time to the current time so that other ADJUST and EXECUTE_AT_END Macros are used
 *********************************************************************/
DEFINE_ON_DEMAND(pre_transient)
{
	/* Initialized artery outlet flow rate calculation */
	if (myid == node_id_artery)
	{
		Domain *d_e = Get_Domain(1);
		Thread *f_thread_artery = Lookup_Thread(d_e,thread_out_artery_ID);

		face_t f;

		Q_n_Ar = 0.;

		begin_f_loop(f,f_thread_artery)
		{
			Q_n_Ar += F_FLUX(f,f_thread_artery)/den;
		}
		end_f_loop(f,f_thread_artery)

		Q_np1_Ar = Q_n_Ar;

		t_np1 = CURRENT_TIME;

		Message("Node %d = myid = %d\n", node_id_artery, myid);
		Message("CURRENT TIME:			Initial Massflow 1:				P_np1_Ar:\n");
		Message("    %e		%16.16e    %16.16e \n", t_np1, Q_n_Ar, P_np1_Ar);
	}

	/* Initialized vein outlet flow rate calculation */
	if (myid == node_id_vein)
	{
		Domain *d_e = Get_Domain(1);
		Thread *f_thread_vein = Lookup_Thread(d_e,thread_out_vein_ID);

		face_t f;

		Q_n_Ve = 0.;

		begin_f_loop(f,f_thread_vein)
		{
			Q_n_Ve += F_FLUX(f,f_thread_vein)/den;
		}
		end_f_loop(f,f_thread_vein)

		Q_np1_Ve = Q_n_Ve;

		t_np1 = CURRENT_TIME;

		Message("Node %d = myid = %d\n", node_id_vein, myid);
		Message("CURRENT TIME:			Initial Massflow 1:				P_np1_Ve:\n");
		Message("    %e		%16.16e    %16.16e \n", t_np1, Q_n_Ve, P_np1_Ve);
	}
}

 /********************************************************************
   Macro for adjusting the pressure at the artery and vein outlets based on the Windkessel algorithm using the Crank–Nicolson method
 *********************************************************************/

DEFINE_ADJUST(my_adjust,d)
{
	/* Artery outlet pressure adjust */
	if (myid == node_id_artery && t_np1 != -1)
	{
		Thread *f_thread_artery = Lookup_Thread(d,thread_out_artery_ID); 

		face_t f;
		Q_np1_Ar = 0;

		begin_f_loop(f,f_thread_artery)
		{
			Q_np1_Ar += F_FLUX(f,f_thread_artery)/den;
		}
		end_f_loop(f,f_thread_artery)
		
		//Message("Q_np1_Ar:		    P_np1_Ar:\n");
		P_np1_Ar = Rp_Ar*(Q_np1_Ar);
		//Message("%16.16e    %16.16e\n",Q_np1_Ar, P_np1_Ar);
	}

	/* Vein outlet pressure adjust */
	if (myid == node_id_vein && t_np1 != -1)
	{
		Thread *f_thread_vein = Lookup_Thread(d,thread_out_vein_ID); 

		face_t f;
		Q_np1_Ve = 0;

		begin_f_loop(f,f_thread_vein)
		{
			Q_np1_Ve += F_FLUX(f,f_thread_vein)/den;
		}
		end_f_loop(f,f_thread_vein)
		
		//Message("Q_np1_Ve:		    P_np1_Ve:\n");
		P_np1_Ve = Rp_Ve*(Q_np1_Ve);
		//Message("%16.16e    %16.16e\n",Q_np1_Ve, P_np1_Ve);
	}
} 

 /********************************************************************
   Macro for applying the pressure at the artery outlet as calculated in the ADJUST Macro
 *********************************************************************/
DEFINE_PROFILE(pressure_out_artery,t,i)
{
	if (myid == node_id_artery)
	{
		face_t f;
		p_var = i; //error compiling when this is above the previous statement

		begin_f_loop(f,t)
			{
				F_PROFILE(f,t,i) = P_np1_Ar; /* Applies pressure to face */
			}
		end_f_loop(f,t)
		//Message("%d - Pressure out 1: %e	p_var: %d\n", node_id ,P_np1, p_var);
	}
}

 /********************************************************************
   Macro for applying the pressure at the vein outlet as calculated in the ADJUST Macro
 *********************************************************************/
DEFINE_PROFILE(pressure_out_vein,t,i)
{
	if (myid == node_id_vein)
	{
		face_t f;
		p_var = i; //error compiling when this is above the previous statement

		begin_f_loop(f,t)
			{
				F_PROFILE(f,t,i) = P_np1_Ve; /* Applies pressure to face */
			}
		end_f_loop(f,t)
		//Message("%d - Pressure out 1: %e	p_var: %d\n", node_id ,P_np1, p_var);
	}
}

 /********************************************************************
   Macro for setting the previous flow rates, pressures and h to current values once the coupling iterations for the current time step have finished
 *********************************************************************/

DEFINE_EXECUTE_AT_END(execute_at_end)
{
	if (myid == node_id_artery && t_np1 != -1)
	{
		Q_n_Ar = Q_np1_Ar;

		Message("t_np1_Ar:			  Q_n_Ar:			P_n_Ar:\n");
		Message("%e		%8.8e		%8.8e\n", t_np1, Q_n_Ar, P_np1_Ar);
	}

	if (myid == node_id_vein && t_np1 != -1)
	{
		Q_n_Ve = Q_np1_Ve;
		t_np1 = t_np1 + dt;

		Message("t_np1_Ve:			  Q_n_Ve:				P_n_Ve:\n");
		Message("%e 	%8.8e		%8.8e\n", t_np1, Q_n_Ve, P_np1_Ve);
	}
}

 /********************************************************************
   Macro for resetting all the variables and parameters to start again
 *********************************************************************/
DEFINE_ON_DEMAND(restart)
{
den = 1060;
dt = 0.01;
NV_VEC(A);
Area = 0.;
t_np1 = -1;

/* Artery parameters*/

Rp_Ar = 2.231840472085634e+09;		/* Distal Resistance (Far)*/

Q_n_Ar = 0.;
Q_np1_Ar = 0.;
P_n_Ar = 6.125e+03;
P_np1_Ar = 6.125e+03;

/* Vein parameters*/

Rp_Ve = 9.206341947353240e+08;		/* Distal Resistance (Far)*/

Q_n_Ve = 0.;
Q_np1_Ve = 0.;
P_n_Ve = 6.025e+03;
P_np1_Ve = 6.025e+03;
}