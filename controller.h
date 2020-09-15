#ifndef _CONTROLLER_H_
#define _CONTROLLER_H_

#include <stdio.h>
#include <stdlib.h>


////////////// VERY IMPORTANT: TO CALIBRATE THIS VALUE EACH TIME //////////////

#define RHO 1.0224				//Air density [kg/m^3]

/////////////////////////////////////////////////////////////////////////////////


//Define
#define CLOCK_XMOS 100000000	//Xmos clock frequency
#define DT 0.01					//Time step [s]
#define DT_2 (DT/2)				//Half time step [s]
#define RAD 0.01745329	    	//Conversion to radians
#define	PI	3.141592			//180 deg
#define	TWO_PI 6.283185			//360 deg
#define G 9.81					//Gravity acceleration [m/s^2]


#define PWM_MIN 1000			//Minimum pwm radio signal [microsec]
#define PWM_MEAN 1500			//Mean pwm radio signal [microsec]
#define PWM_MAX 2000			//Max pwm radio signal [microsec]


#define ND 100					//Filter coefficient of the derivative action (as in Simulink)

#define ERROR_VALUES(reference,actual) ((reference)-(actual))
#define DERIVATIVE(err,prev_err) (((err)-(prev_err))/DT)
#define INTEGRAL_ERROR(err,prev_err,prev_int_err) ((prev_int_err)+((err)+(prev_err))*DT/2)
#define PID(err,int_err,P,I,der_action) (((P)*(err))+((I)*(int_err))+(der_action))


typedef struct
{
	int da;
	int de;
	int dr;
	int th;
	int err_V;
	int err_h;
	int err_psi;
	int int_err_V;
	int int_err_h;
	int int_err_psi;
	int der_action_V;
	int der_action_h;
	int der_action_psi;
	int ps;
	int qs;
	int rs;
	int alpha;
	int beta;
	int theta;
	int phi;
	int g2;
	int g3;
	int p;
	int	q;
	int r;
	int t;
}variables;


void control(int i, int dt, int V_in, int h_in, int psi_in, variables prev[]);


#endif
