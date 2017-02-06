#include "udf.h"
#include "unsteady.h"

//define IDOUT1 15

int p_var;								/* Initialize profile macro variable */
int no_face;							/* Initialize no. of faces integer variable*/
int thread_out_artery_ID = 6;			/* Initialize current time step artery flow rate */
int thread_out_vein_ID = 8;				/* Initialize current time step artery flow rate */

int WindInit = 0;
int WindAdjust = 0;

real den = 1060;						/* Blood density */
real dt = 0.01;							/* time step size */
real NV_VEC(A);							/* Initialize vector A for area calculations */
real Area = 0.;							/* Initialize Area for area calculations */

/* Artery outlet parameters and variables*/
int node_id_artery;						/* Initialize node id for artery */

real Rp_Ar = 1.729e+09;					/* Distal Resistance (Far)*/
real Rd_Ar = 8.648e+08;					/* Proximal Resistance (Near)*/
real C_Ar  = 3.336e-10;

real Q_n_Ar = 0.;						/* Initialize current time step artery flow rate */
real Q_np1_Ar = 0.;						/* Initialize previous time step artery flow rate */
real P_np1_Ar = 4.5447451;				/* Initialize previous time step artery flow rate */

real tau_Ar = 0;
real P_inst_Ar = 4.5447451;
real h_np1_Ar = 0;
real h_n_Ar = 0;

/* Vein outlet parameters and variables*/
int node_id_vein;						/* Initialize node id for vein */

real Rp_Ve = 6.874e+08;					/* Distal Resistance (Far)*/
real Rd_Ve = 4.124e+08;					/* Proximal Resistance (Near)*/
real C_Ve  = 1.399e-09;

real Q_n_Ve = 0.;						/* Initialize current time step vein flow rate */
real Q_np1_Ve = 0.;						/* Initialize previous time step vein flow rate */
real P_np1_Ve = 4.5102368;				/* Initialize previous time step vein flow rate */

real tau_Ve = 0;
real P_inst_Ve = 4.5102368;
real h_np1_Ve = 0;
real h_n_Ve = 0;

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

	tau_Ar = Rd_Ar*C_Ar;
	tau_Ve = Rd_Ve*C_Ve;

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
		Message("ARTERY \n - Initial Pressure (%e)Pa \n", P_np1_Ar);
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
		Message("VEIN \n - Initial Pressure (%e)Pa \n", P_np1_Ve);
		Message(" - Outlet is on node (%d)\n", myid);
		Message(" - Outlet is made up of (%e) faces with area (%8.8e)\n", no_face , Area);
	}
} 

 /********************************************************************
   Macro that calculates the initialized flow rates of the artery and vein outlets and sets the time to the current time so that other ADJUST and EXECUTE_AT_END Macros are used
 *********************************************************************/
DEFINE_ON_DEMAND(pre_windinit)
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

		WindInit = 1;

		Message("Node %d = myid = %d\n", node_id_artery, myid);
		Message("CURRENT TIME:			Initial Massflow 1:				P_np1_Ar:\n");
		Message("%16.16e    %16.16e \n", Q_n_Ar, P_np1_Ar);
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

		WindInit = 1;

		Message("Node %d = myid = %d\n", node_id_vein, myid);
		Message("CURRENT TIME:			Initial Massflow 1:				P_np1_Ve:\n");
		Message("%16.16e    %16.16e \n", Q_n_Ve, P_np1_Ve);
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

		Message("Node %d = myid = %d\n", node_id_artery, myid);
		Message("CURRENT TIME:			Initial Massflow 1:				P_np1_Ar:\n");
		Message("%16.16e    %16.16e \n", Q_n_Ar, P_np1_Ar);

		WindInit = 0;
		WindAdjust = 1;
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

		WindInit = 0;
		WindAdjust = 1;

		Message("Node %d = myid = %d\n", node_id_vein, myid);
		Message("CURRENT TIME:			Initial Massflow 1:				P_np1_Ve:\n");
		Message("%16.16e    %16.16e \n", Q_n_Ve, P_np1_Ve);
	}
}
 /********************************************************************
   Macro for adjusting the pressure at the artery and vein outlets based on the Windkessel algorithm using the Crank–Nicolson method
 *********************************************************************/

DEFINE_ADJUST(my_adjust,d)
{
	if (WindInit == 1)
	{
		/* Artery outlet pressure adjust */
		if (myid == node_id_artery)
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
		if (myid == node_id_vein)
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
	
	if (WindAdjust == 1)
	{
		/* Artery outlet pressure adjust */
		if (myid == node_id_artery)
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
			P_inst_Ar = Rp_Ar*(Q_np1_Ar);
			h_np1_Ar = (exp(-dt/tau_Ar))*(h_n_Ar+.5*Q_n_Ar*dt)+.5*Q_np1_Ar*dt;
			P_np1_Ar = P_inst_Ar + 1/C_Ar*h_np1_Ar;
			//Message("%16.16e    %16.16e\n",Q_np1_Ar, P_np1_Ar);
		}

		/* Vein outlet pressure adjust */
		if (myid == node_id_vein)
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
			P_inst_Ve = Rp_Ve*(Q_np1_Ve);
			h_np1_Ve = (exp(-dt/tau_Ve))*(h_n_Ve+.5*Q_n_Ve*dt)+ .5*Q_np1_Ve*dt;
			P_np1_Ve = P_inst_Ve + 1/C_Ve*h_np1_Ve;
			//Message("%16.16e    %16.16e\n",Q_np1_Ve, P_np1_Ve);
		}
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
	if (WindAdjust == 1)
	{
		if (myid == node_id_artery)
		{
			h_n_Ar = h_np1_Ar;
			Q_n_Ar = Q_np1_Ar;

			Message("Q_n_Ar:			P_np1_Ar:\n");
			Message("%8.8e		%8.8e\n", Q_n_Ar, P_np1_Ar);
		}

		if (myid == node_id_vein)
		{
			h_n_Ve = h_np1_Ve;
			Q_n_Ve = Q_np1_Ve;

			Message("Q_n_Ve:				P_np1_Ve:\n");
			Message("%8.8e		%8.8e\n", Q_n_Ve, P_np1_Ve);
		}
	}
}

 /********************************************************************
   Macro for resetting all the variables and parameters to start again
 *********************************************************************/
DEFINE_ON_DEMAND(restart)
{
WindInit = 0;
WindAdjust = 0;

/* Artery parameters*/

Rp_Ar = 1.729e+09;		/* Distal Resistance (Far)*/
Rd_Ar = 8.648e+08;		/* Proximal Resistance (Near)*/
C_Ar  = 3.336e-10;

Q_n_Ar = 0.;
Q_np1_Ar = 0.;
P_np1_Ar = 4.5447451;

tau_Ar = 0;
P_inst_Ar = 4.5447451;
h_np1_Ar = 0;
h_n_Ar = 0;

/* Vein parameters*/

Rp_Ve = 6.874e+08;		/* Distal Resistance (Far)*/
Rd_Ve = 4.124e+08;		/* Proximal Resistance (Near)*/
C_Ve  = 1.399e-09;

Q_n_Ve = 0.;
Q_np1_Ve = 0.;
P_np1_Ve = 4.5102368;

tau_Ve = 0;
P_inst_Ve = 4.5102368;
h_np1_Ve = 0;
h_n_Ve = 0;
}