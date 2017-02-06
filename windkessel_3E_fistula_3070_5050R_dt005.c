#include "udf.h"
#include "unsteady.h"

int p_var;								/* Initialize profile macro variable */
int no_face;							/* Initialize no. of faces integer variable*/
int thread_out_artery_ID = 30;			/* Initialize current time step artery flow rate */
int thread_out_vein_ID = 33;				/* Initialize current time step artery flow rate */

real den = 1050;						/* Blood density */
real dt = 0.005;							/*Time step size*/
real NV_VEC(A);							/* Initialize vector A for area calculations */
real Area = 0.;							/* Initialize Area for area calculations */
real stage = -1;					/* Halt time stepping unitl pre-transient macro used */
FILE *fp = NULL;
int count_i=0;

/* Artery outlet parameters and variables*/
int node_id_artery;						/* Initialize node id for artery */

real Rp_Ar_init = 2.231840472085634e+09;		/* Artery Proximal Resistance (Near)*/
real Rp_Ar = 1.116E+09;		/* Artery Proximal Resistance (Near)*/
real Rd_Ar = 1.116E+09;		/* Distal Resistance (Far)*/
real C_Ar  = 1.34E-10;
real tau_Ar = 0;

real Q_n_Ar = 0.;
real Q_np1_Ar = 0.;
real P_np1_Ar = 6.125e+03;
						/* Initialize current time step artery flow rate */

real P_inst_Ar = 0;
real h_np1_Ar = 0;
real h_n_Ar = 0;

/* Vein outlet parameters and variables*/
int node_id_vein;						/* Initialize node id for vein */

real Rp_Ve_init = 9.206341947353240e+08;		/* Vein Proximal Resistance (Near)*/
real Rp_Ve = 4.603E+08;		/* Artery Proximal Resistance (Near)*/
real Rd_Ve = 4.603E+08;		/* Distal Resistance (Far)*/
real C_Ve  = 1.74E-10;
real tau_Ve = 0;

real Q_n_Ve = 0.;
real Q_np1_Ve = 0.;
real P_np1_Ve = 6.025e+03;				/* Initialize current time step vein flow rate */


real P_inst_Ve = 0;
real h_np1_Ve = 0;
real h_n_Ve = 0;

 /***********************************************************************
	Macro - finds which threads are on which nodes and checking the area of the artery and vein outlets
 ************************************************************************/

DEFINE_ON_DEMAND(pre_initialize)	
{
	Domain *d_e = Get_Domain(1);
	Thread *f_thread_artery = Lookup_Thread(d_e,thread_out_artery_ID);
	Thread *f_thread_vein = Lookup_Thread(d_e,thread_out_vein_ID);

	face_t f;
	node_id_artery = 111;			/* to ensure artery node id dosn't auto to 0 */	
	node_id_vein = 111;			/* to ensure vein node id dosn't auto to 0*/

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
		Message("ARTERY \n - Starting Pressure (%d)Pa \n", P_np1_Ar);
		Message(" - Outlet is on node (%d)\n", myid);
		Message(" - Outlet is made up of (%d) faces with area (%8.8e)\n", no_face , Area);
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
		Message("VEIN \n - Starting Pressure (%d)Pa \n", P_np1_Ve);
		Message(" - Outlet is on node (%d)\n", myid);
		Message(" - Outlet is made up of (%d) faces with area (%8.8e)\n", no_face , Area);
	}
} 

 /********************************************************************
   Macro - find and set initialsing values
 *********************************************************************/

DEFINE_ON_DEMAND(windkessel1_on)
{
	/* Initialized artery outlet flow rate calculation */
	if (myid == node_id_artery)
	{
		Domain *d_e = Get_Domain(1);
		Thread *f_thread_artery = Lookup_Thread(d_e,thread_out_artery_ID);

		face_t f;
		fp = fopen ("count_iter.txt", "w");
		fclose(fp);
		
		Q_n_Ar = 0.;

		begin_f_loop(f,f_thread_artery)
		{
			Q_n_Ar += F_FLUX(f,f_thread_artery)/den;
		}
		end_f_loop(f,f_thread_artery)
		
		P_np1_Ar = Rp_Ar_init*(Q_n_Ar);
		stage = 0;

		Message("Node %d = myid = %d\n", node_id_artery, myid);
		Message("CURRENT TIME:			Initial Massflow 1:				P_np1_Ar:\n");
		Message("    %e		%16.16e    %16.16e \n", stage, Q_n_Ar, P_np1_Ar);
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

		P_np1_Ve = Rp_Ve_init*(Q_n_Ve);
		stage = 0;

		Message("Node %d = myid = %d\n", node_id_vein, myid);
		Message("CURRENT TIME:			Initial Massflow 1:				P_np1_Ve:\n");
		Message("    %e		%16.16e    %16.16e \n", stage, Q_n_Ve, P_np1_Ve);
	}
}

 /********************************************************************
   Macro - find and set initialsing values
 *********************************************************************/

DEFINE_ON_DEMAND(windkessel3_on)
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
		
		tau_Ar = Rd_Ar*C_Ar;
		P_np1_Ar = Rp_Ar_init*(Q_n_Ar);
		/*h_n_Ar = (exp(-dt/(tau_Ar))*dt*Q_n_Ar+dt*Q_n_Ar)/(2-2*exp(-dt/(tau_Ar))); */
		/*h_np1_Ar = h_n_Ar; */
		stage = 1;

		Message("Node %d = myid = %d\n", node_id_artery, myid);
		Message("CURRENT TIME:			Initial Massflow 1:				P_np1_Ar:\n");
		Message("    %e		%16.16e    %16.16e \n", stage, Q_n_Ar, P_np1_Ar);
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
		
		tau_Ve = Rd_Ve*C_Ve;
		P_np1_Ve = Rp_Ve_init*(Q_n_Ve);
		/*h_n_Ve = (exp(-dt/(tau_Ve))*dt*Q_n_Ve+dt*Q_n_Ve)/(2-2*exp(-dt/(tau_Ve)));*/
		/*h_np1_Ve = h_n_Ve;*/
		stage = 1;

		Message("Node %d = myid = %d\n", node_id_vein, myid);
		Message("CURRENT TIME:			Initial Massflow 1:				P_np1_Ve:\n");
		Message("    %e		%16.16e    %16.16e \n", stage, Q_n_Ve, P_np1_Ve);
	}
}

 /********************************************************************
   Macro for adjusting the pressure at the artery and vein outlets based on the 
	Windkessel algorithm using the CrankNicolson method
 *********************************************************************/

DEFINE_ADJUST(adjust_outlet_pressures,d)
{
	/* Adjusts the artery outlet pressure when initialising */
	if (myid == node_id_artery && stage > -0.5 && stage < 0.5)
	{
		Thread *f_thread_artery = Lookup_Thread(d,thread_out_artery_ID); 

		face_t f;
		Q_n_Ar = 0;

		begin_f_loop(f,f_thread_artery)
		{
			Q_n_Ar += F_FLUX(f,f_thread_artery)/den;
		}
		end_f_loop(f,f_thread_artery)
		
		count_i++;
		P_np1_Ar = Rp_Ar_init*(Q_n_Ar);
	}
	/* Adjusts the artery outlet pressure for 3 element windkessel */
	else if (myid == node_id_artery && stage > 0.5 && stage < 1.5)
	{
		Thread *f_thread_artery = Lookup_Thread(d,thread_out_artery_ID); 

		face_t f;
		Q_np1_Ar = 0;

		begin_f_loop(f,f_thread_artery)
		{
			Q_np1_Ar += F_FLUX(f,f_thread_artery)/den;
		}
		end_f_loop(f,f_thread_artery)
		
		count_i++;
		/*Message("Q_np1_Ar:		  P_inst_Ar:		  h_np1_Ar:		    P_np1_Ar:\n");*/
		P_inst_Ar = Rp_Ar*(Q_np1_Ar);
		h_np1_Ar = (exp(-dt/tau_Ar))*(h_n_Ar+.5*Q_n_Ar*dt)+.5*Q_np1_Ar*dt;
		P_np1_Ar = P_inst_Ar + 1/C_Ar*h_np1_Ar;

		/*Message("%16.16e    %16.16e    %16.16e    %16.16e\n",Q_np1_Ar, P_inst_Ar, h_np1_Ar, P_np1_Ar);*/
	}

	/* Adjusts the vein outlet pressure when initialising */
	if (myid == node_id_vein && stage > -0.5 && stage < 0.5)
	{
		Thread *f_thread_vein = Lookup_Thread(d,thread_out_vein_ID); 

		face_t f;
		Q_n_Ve = 0;

		begin_f_loop(f,f_thread_vein)
		{
			Q_n_Ve += F_FLUX(f,f_thread_vein)/den;
		}
		end_f_loop(f,f_thread_vein)
		
		P_np1_Ve = Rp_Ve*(Q_n_Ve);
	}
	/* Adjusts the vein outlet pressure for 3 element windkessel */
	else if (myid == node_id_vein && stage > 0.5 && stage < 1.5)
	{
		Thread *f_thread_vein = Lookup_Thread(d,thread_out_vein_ID); 

		face_t f;
		Q_np1_Ve = 0;

		begin_f_loop(f,f_thread_vein)
		{
			Q_np1_Ve += F_FLUX(f,f_thread_vein)/den;
		}
		end_f_loop(f,f_thread_vein)
		
		/*Message("Q_np1_Ve:		  P_inst_Ve:  		 h_np1_Ve:		    P_np1_Ve:\n");*/
		P_inst_Ve = Rp_Ve*(Q_np1_Ve);
		h_np1_Ve = (exp(-dt/tau_Ve))*(h_n_Ve+.5*Q_n_Ve*dt)+ .5*Q_np1_Ve*dt;
		P_np1_Ve = P_inst_Ve + 1/C_Ve*h_np1_Ve;

		/*Message("%16.16e    %16.16e    %16.16e    %16.16e\n",Q_np1_Ve, P_inst_Ve, h_np1_Ve, P_np1_Ve);*/
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
		p_var = i;

		begin_f_loop(f,t)
			{
				F_PROFILE(f,t,i) = P_np1_Ar; /* Applies pressure to face */
			}
		end_f_loop(f,t)
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
		p_var = i;

		begin_f_loop(f,t)
			{
				F_PROFILE(f,t,i) = P_np1_Ve; /* Applies pressure to face */
			}
		end_f_loop(f,t)
	}
}

/********************************************************************
   Macro for counting stuff
 *********************************************************************/
DEFINE_EXECUTE_AT_END(count_coup_iter)
{
if (myid == node_id_artery)
	{
		fp = fopen ("count_iter.txt", "a+");
		fprintf(fp, "%d\n", count_i);
		fclose(fp);
		count_i = 0;
	}
if (myid == node_id_artery && stage > 0.5 && stage < 1.5)
	{
		h_n_Ar = h_np1_Ar;
		Q_n_Ar = Q_np1_Ar;

		Message("Q_n_Ar:			h_n_Ar:\n");
		Message("%16.16e		%16.16e\n", Q_n_Ar, h_n_Ar);
	}
if (myid == node_id_vein && stage > 0.5 && stage < 1.5)
	{
		h_n_Ve = h_np1_Ve;
		Q_n_Ve = Q_np1_Ve;

		Message("Q_n_Ve:			h_n_Ve:\n");
		Message("%16.16e		%16.16e	\n", Q_n_Ve, h_n_Ve);
	}
}