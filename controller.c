
//BACKSTEPPING CONTROLLER FOR UAV


//Library definition
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "controller.h"


//Functions declaration
float derivative_action(float D, float der_err, float prev_der_action);
void backstepping_controller(int i, float alpha_ref, float p_ref, float V, float th, float err_psi, variables prev[]);
float calculation_ps(float ps_ref, float prev_ps);
float calculation_qs(float prev_qs, float alpha, float alpha_ref, float ps, float beta, float V, float T, float g2);
float calculation_rs(float prev_rs, float theta, float phi, float V, float beta);
float calculation_g2(float alpha, float theta, float phi);
float calculation_g3(float alpha, float theta, float phi, float beta);
float calculation_theta(float prev_theta, float q, float r, float phi);
float calculation_phi(float prev_phi, float theta, float p, float q, float r);
float calculation_alpha(float prev_alpha, float ps, float qs, float beta, float V, float g2, float T);
float calculation_beta(float prev_beta, float rs, float alpha, float V, float g3, float T);
float calculation_L(float p, float q, float r, float der_p, float der_r);
float calculation_M(float p, float r, float der_q);
float calculation_N(float p, float q, float r, float der_p, float der_r);
int calculation_da(float L, float V, float p, float r, float beta);
int calculation_de(float M, float V, float q, float alpha);
int calculation_dr(float N, float V, float p, float r, float beta);
float command_saturation(float d, float max);
int command_pwm(float d, float c1, float c2, float c3, float c4);


//Union prototype for conversion from int to float
typedef union
{
	int ii;
	float ff;

}f2i;

//Temporary store f2i variables
f2i temporary;


//Aircraft characteristics
const float m=1.959;				//Aircraft mass [kg]
const float S=0.3097;				//Wing surface [m^2]
const float b=1.27;					//Wingspan [m]
const float c=0.25;					//Aerodynamic mean chord [m]

const float CL_0=0.1086;			//Lift coefficient for zero incidence
const float CL_alpha=4.58;			//Lift coefficient derivative [rad^-1]
const float CY_beta=-0.83;			//Lateral force coefficient derivative [rad^-1]
const float Cl_beta=-0.0545;		//Roll moment coefficient derivative [rad^-1]
const float Cl_p=-0.4496;			//Roll moment coefficient derivative [rad^-1]
const float Cl_r=0.1086;			//Roll moment coefficient derivative [rad^-1]
const float Cl_da=-0.1646;			//Aileron command derivative [rad^-1]
const float CM_alpha=-0.723;		//Pitch moment coefficient derivative [rad^-1]
const float CM_q=-13.5664;			//Pitch moment coefficient derivative [rad^-1]
const float CM0=-0.0278;			//Pitch moment coefficient derivative [rad^-1]
const float CM_de=-0.8488;			//Elevator command derivative [rad^-1]
const float CN_beta=0.0723;			//Yaw moment coefficient derivative [rad^-1]
const float CN_p=0.118;				//Yaw moment coefficient derivative [rad^-1]
const float CN_r=-0.1833;			//Yaw moment coefficient derivative [rad^-1]
const float CN_dr=-0.1811;			//Yaw command derivative [rad^-1]

const float Ix=0.0715;				//Inertia [kg*m^2]
const float Iy=0.0864;				//Inertia [kg*m^2]
const float Iz=0.15364;				//Inertia [kg*m^2]
const float Ixz=0.014;				//Inertia [kg*m^2]

const float Pmax=600*0.8;			//Maximum engine power (multiplied by typical prop efficiency) [W]
const float th_max=1;				//Maximum engine throttle
const float th_min=0.1;				//Minimum engine throttle
float da_max=28*RAD;				//Maximum/minimum aileron deflection [rad]
float de_max=26*RAD;				//Maximum/minimum elevator deflection [rad]
float dr_max=26*RAD;				//Maximum/minimum rudder deflection [rad]
const float pitch_sf=0.6;			//Pitch angle safety factory (reduces maximum elevator deflection)
const float roll_sf=0.4;			//Bank angle safety factory (reduces maximum aileron deflection)
const float yaw_sf=0.8;				//Yaw angle safety factory (reduces maximum rudder deflection)


//Initial flight conditions
const float V0=17;
const float h0=1638;
const float	psi0=0*RAD;

//Trim conditions
const float de_trim=-5.5056*RAD;	//Trim elevator deflection [rad] (for these initial conditions)
const float da_trim=0.575*RAD;		//Trim aileron deflection [rad] (for these initial conditions)
const float dr_trim=0.1816*RAD;	    //Trim rudder deflection [rad] (for these initial conditions)
const float th_trim=0.5703;		    //Trim throttle (for these initial conditions)




//Backstepping controller gains
const float Kps=1, K_alpha_1=10, K_alpha_2=30, K_beta_1=5, K_beta_2=15;


//CUTOFF FREQUENCY 0.1
const float P_V=-0.0026, I_V=-0.004, D_V=0.028;
const float P_h=0.03, I_h=0.0012, D_h=0.027;
const float P_psi=0.25, I_psi=0.08, D_psi=2.7;




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////MAIN FUNCTION///////////////////////////////////////////////////////////////////////////////////////////////////////
//What main function does:
// - defines simulation parameters,
// - defines the reference signals,
// - gets feedbacks from aircraft,
// - gives the signals to the controller,
// - gives commands to the aircraft.



void control(int i, int dt, int V_in, int h_in, int psi_in, variables prev[])

{

	//Variables and constants definition//


	//Throttle
	float th;

	//Reference values for the controller
	float alpha_ref;		//alpha[rad],
	float p_ref;			// p[rad/s]

	//Reference values for aircraft states
	float V_ref;			//[m/s]
	float h_ref;			//[m]
	float psi_ref;			//[rad]


	//Errors
	float err_V, err_h, err_psi;
	float prev_err_V, prev_err_h, prev_err_psi;

	//Integral and derivative of errors
	float int_err_psi, int_err_V, int_err_h;
	float prev_int_err_psi, prev_int_err_V, prev_int_err_h;
	float der_err_psi, der_err_V, der_err_h;

	//Derivative action
	float der_action_V, der_action_h, der_action_psi;
	float prev_der_action_V, prev_der_action_h, prev_der_action_psi;

	//Feedback variables
	float V, h, psi;

	//Simulation time
	float t;



	//Computation starts//


	//Calculate time

	//Value from incoming struct
	temporary.ii=prev[0].t;
	t=temporary.ff;

	//Initialization
	if (i==0)
	{
		t=0;
	}

	//Update time
	t=t+(float)dt/CLOCK_XMOS;

	//Structure for output
	temporary.ff=t;
	prev[0].t=temporary.ii;



	//Reference values 
	V_ref=V0;
	h_ref=h0;
	psi_ref=psi0;


	//Assign incoming data
	temporary.ii=V_in;
	V=temporary.ff;

	//Control on V=0 (to avoid NaN generation in case of pitot failure or activation on the ground)
	if (V<1)
	{
		V=V_ref;
	}


	temporary.ii=h_in;
	h=temporary.ff;

	psi=(float)psi_in*RAD;


	//Calculate error
	err_psi=ERROR_VALUES(psi_ref,psi);

	//Define minimal error direction
	if (err_psi<-PI)
	{
		err_psi=err_psi+TWO_PI;
	}
	else if (err_psi>PI)
	{
		err_psi=err_psi-TWO_PI;
	}

	//Calculate errors
	err_V=ERROR_VALUES(V_ref,V);
	err_h=ERROR_VALUES(h_ref,h);



		//Calculate integral and derivative of error

		if (i==0)
		{
			prev_err_V=0;
			prev_err_h=0;
			prev_err_psi=0;

			prev_der_action_V=0;
			prev_der_action_h=0;
			prev_der_action_psi=0;

			int_err_psi=0;
			int_err_V=0;
			int_err_h=0;

			der_err_psi=0;
			der_err_V=0;
			der_err_h=0;

			der_action_psi=ND*D_psi*err_psi;
			der_action_V=ND*D_V*err_V;
			der_action_h=ND*D_h*err_h;
		}

		else

		{
			//Values from incoming struct
			temporary.ii=prev[0].err_V;
			prev_err_V=temporary.ff;
			temporary.ii=prev[0].err_h;
			prev_err_h=temporary.ff;
			temporary.ii=prev[0].err_psi;
			prev_err_psi=temporary.ff;

			temporary.ii=prev[0].der_action_V;
			prev_der_action_V=temporary.ff;
			temporary.ii=prev[0].der_action_h;
			prev_der_action_h=temporary.ff;
			temporary.ii=prev[0].der_action_psi;
			prev_der_action_psi=temporary.ff;

			temporary.ii=prev[0].int_err_V;
			prev_int_err_V=temporary.ff;
			temporary.ii=prev[0].int_err_h;
			prev_int_err_h=temporary.ff;
			temporary.ii=prev[0].int_err_psi;
			prev_int_err_psi=temporary.ff;

			int_err_psi=INTEGRAL_ERROR(err_psi,prev_err_psi,prev_int_err_psi);
			int_err_V=INTEGRAL_ERROR(err_V,prev_err_V,prev_int_err_V);
			int_err_h=INTEGRAL_ERROR(err_h,prev_err_h,prev_int_err_h);

			der_err_psi=DERIVATIVE(err_psi,prev_err_psi);
			der_err_V=DERIVATIVE(err_V,prev_err_V);
			der_err_h=DERIVATIVE(err_h,prev_err_h);

			der_action_psi=derivative_action(D_psi,der_err_psi,prev_der_action_psi);
			der_action_V=derivative_action(D_V,der_err_V,prev_der_action_V);
			der_action_h=derivative_action(D_h,der_err_h,prev_der_action_h);
		}



		//Structure for output
		temporary.ff=err_V;
		prev[0].err_V=temporary.ii;
		temporary.ff=err_h;
		prev[0].err_h=temporary.ii;
		temporary.ff=err_psi;
		prev[0].err_psi=temporary.ii;

		temporary.ff=int_err_V;
		prev[0].int_err_V=temporary.ii;
		temporary.ff=int_err_h;
		prev[0].int_err_h=temporary.ii;
		temporary.ff=int_err_psi;
		prev[0].int_err_psi=temporary.ii;

		temporary.ff=der_action_V;
		prev[0].der_action_V=temporary.ii;
		temporary.ff=der_action_h;
		prev[0].der_action_h=temporary.ii;
		temporary.ff=der_action_psi;
		prev[0].der_action_psi=temporary.ii;



		//Call PIDs
		p_ref=PID(err_psi,int_err_psi,P_psi,I_psi,der_action_psi);
		alpha_ref=PID(err_V,int_err_V,P_V,I_V,der_action_V);
		th=th_trim + PID(err_h,int_err_h,P_h,I_h,der_action_h);

		//Limit on throttle
		if (th>th_max)
		{
			th=th_max;
		}
		else if (th<th_min)
		{
			th=th_min;
		}

		//Throttle output in pwm
		prev[0].th=PWM_MIN+(PWM_MAX-PWM_MIN)*th;


		//Call backstepping controller
		backstepping_controller(i,alpha_ref,p_ref,V,th,err_psi,prev);


}





///////////////DERIVATIVE_ACTION/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//What derivative_action function does:
// - receives the derivative of the error, the derivative gain and the filter coefficient,
// - calculates the integral of the derivative action.

float derivative_action(float D, float der_err, float prev_der_action)
{
	float f1, D_der;		//Auxiliary expression
	float yt;				//Intermediate parameter of integration
	float der_action;


	D_der=D*der_err;
	f1=ND*(D_der-prev_der_action);

	yt=prev_der_action+DT*f1;
	der_action=prev_der_action+DT_2*(f1+ND*(D_der-yt));


	return der_action;
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////BACKSTEPPING_CONTROLLER////////////////////////////////////////////////////////////////////////////////////////
//What backstepping_controller function does:
// - receives the reference vectors of alpha and p,
// - calculates the controls de and da.



void backstepping_controller(int i, float alpha_ref, float p_ref, float V, float th, float err_psi, variables prev[])
{

	float alpha;		//Angle of attack [rad]
	float beta; 		//Sideslip attack [rad]

	float phi; 			//Roll angle [rad]
	float theta; 		//Pitch angle [rad]
	float ps;			//Roll rate in stability axes [rad/s]
	float qs;			//Pitch rate in stability axes [rad/s]
	float rs;			//Yaw rate in stability axes [rad/s]
	float p;			//Roll rate in body axes [rad/s]
	float q;			//Pitch rate in body axes [rad/s]
	float r;			//Yaw rate in body axes [rad/s]
	float g2;			//Gravity component [m/s^2]
	float g3;			//Gravity component [m/s^2]

	int da;				//Aileron command [rad]
	int de;				//Elevator command [rad]
	int dr;				//Rudder command [rad]

	float L;			//Roll moment [N*m]
	float M;			//Pitch moment [N*m]
	float N;			//Yaw moment [N*m]
	float ps_ref;		//Reference roll rate in stability axes [rad/s]
	float T;			//Engine thrust [N]

	float prev_alpha;
	float prev_beta;
	float prev_ps;
	float prev_qs;
	float prev_rs;
	float prev_g2;
	float prev_g3;
	float prev_theta;
	float prev_phi;
	float prev_p;
	float prev_q;
	float prev_r;
	float der_p;
	float der_q;
	float der_r;

	float cos_alpha;
	float sin_alpha;


	//Values from incoming struct
	temporary.ii=prev[0].ps;
	prev_ps=temporary.ff;
	temporary.ii=prev[0].qs;
	prev_qs=temporary.ff;
	temporary.ii=prev[0].rs;
	prev_rs=temporary.ff;

	temporary.ii=prev[0].alpha;
	prev_alpha=temporary.ff;
	temporary.ii=prev[0].beta;
	prev_beta=temporary.ff;

	temporary.ii=prev[0].g2;
	prev_g2=temporary.ff;
	temporary.ii=prev[0].g3;
	prev_g3=temporary.ff;;

	temporary.ii=prev[0].theta;
	prev_theta=temporary.ff;
	temporary.ii=prev[0].phi;
	prev_phi=temporary.ff;

	temporary.ii=prev[0].p;
	prev_p=temporary.ff;
	temporary.ii=prev[0].q;
	prev_q=temporary.ff;
	temporary.ii=prev[0].r;
	prev_r=temporary.ff;


	ps_ref=p_ref*cosf(alpha_ref);		//Transformation of roll rate to stability axes
	T=Pmax*th/V;						//Engine thrust [N]


	if (i==0)
	{
		//Initialization
		alpha=0;
		beta=0;
		phi=0;
		theta=0;
		ps=0;
		qs=0;
		rs=0;
		der_p=0;
		der_q=0;
		der_r=0;

		cos_alpha=cosf(alpha);
		sin_alpha=sinf(alpha);

		p=ps*cos_alpha-rs*sin_alpha;
		q=qs;
		r=ps*sin_alpha+rs*cos_alpha;

	}



	//Calculations
	if (i>0)
	{

		ps=calculation_ps(ps_ref,prev_ps);
		qs=calculation_qs(prev_qs,prev_alpha,alpha_ref,ps,prev_beta,V,T,prev_g2);
		rs=calculation_rs(prev_rs,prev_theta,prev_phi,V,prev_beta);

		cos_alpha=cosf(prev_alpha);
		sin_alpha=sinf(prev_alpha);

		p=ps*cos_alpha-rs*sin_alpha;
		q=qs;
		r=ps*sin_alpha+rs*cos_alpha;

		theta=calculation_theta(prev_theta,q,r,prev_phi);
		phi=calculation_phi(prev_phi,theta,p,q,r);

		alpha=calculation_alpha(prev_alpha,ps,qs,prev_beta,V,prev_g2,T);
		beta=calculation_beta(prev_beta,rs,alpha,V,prev_g3,T);

		der_p=DERIVATIVE(p,prev_p);
		der_q=DERIVATIVE(q,prev_q);
		der_r=DERIVATIVE(r,prev_r);

	}


	g2=calculation_g2(alpha,theta,phi);
	g3=calculation_g3(alpha,theta,phi,beta);

	L=calculation_L(p,q,r,der_p,der_r);
	M=calculation_M(p,r,der_q);
	N=calculation_N(p,q,r,der_p,der_r);


	da=calculation_da(L,V,p,r,beta);
	de=calculation_de(M,V,q,alpha);
	dr=calculation_dr(N,V,p,r,beta);



	//Define structure for output
	prev[0].da=da;
	prev[0].de=de;
	prev[0].dr=dr;

	temporary.ff=ps;
	prev[0].ps=temporary.ii;
	temporary.ff=qs;
	prev[0].qs=temporary.ii;
	temporary.ff=rs;
	prev[0].rs=temporary.ii;

	temporary.ff=alpha;
	prev[0].alpha=temporary.ii;
	temporary.ff=beta;
	prev[0].beta=temporary.ii;

	temporary.ff=theta;
	prev[0].theta=temporary.ii;
	temporary.ff=phi;
	prev[0].phi=temporary.ii;

	temporary.ff=g2;
	prev[0].g2=temporary.ii;
	temporary.ff=g3;
	prev[0].g3=temporary.ii;

	temporary.ff=p;
	prev[0].p=temporary.ii;
	temporary.ff=q;
	prev[0].q=temporary.ii;
	temporary.ff=r;
	prev[0].r=temporary.ii;

}




///////////////CALCULATION_PS////////////////////////////////////////////////////////////////////////////////////////
//What calculation_ps function does:
// - receives the reference roll rate in stability axes,
// - calculates the actual roll rate in stability axes.

float calculation_ps(float ps_ref, float prev_ps)
{
	float f1;		//Auxiliary expression
	float yt;		//Intermediate parameter of integration
	float ps;

	f1=Kps*(ps_ref-prev_ps);

	yt=prev_ps+DT*f1;
	ps=prev_ps+DT_2*(f1+Kps*(ps_ref-yt));

	return ps;
}



///////////////CALCULATION_QS////////////////////////////////////////////////////////////////////////////////////////
//What calculation_qs function does:
// - receives values inside the controller,
// - calculates the pitch rate in stability axes.

float calculation_qs(float prev_qs, float alpha, float alpha_ref, float ps, float beta, float V, float T, float g2)
{
	float f1, f2, f_din, CL_ref, Num;		//Auxiliary expressions
	float yt;								//Intermediate parameter of integration
	float qs, f_alpha, err_alpha;


	f_din=0.5*RHO*V*V*S;
	CL_ref=CL_alpha*alpha_ref+CL_0;
	Num=-f_din*CL_ref-T*sinf(alpha_ref)+m*g2;

	err_alpha=alpha-alpha_ref;
	f_alpha=-ps*tanf(beta)+Num/(m*V*cosf(beta));

	f1=K_alpha_1*err_alpha+f_alpha;
	f2=-K_alpha_2*(prev_qs+f1);

	yt=prev_qs+DT*f2;
	qs=prev_qs+DT_2*(f2-K_alpha_2*(yt+f1));


	return qs;
}



///////////////CALCULATION_RS////////////////////////////////////////////////////////////////////////////////////////
//What calculation_rs function does:
// - receives values inside the controller,
// - calculates the yaw rate in stability axes.

float calculation_rs(float prev_rs, float theta, float phi, float V, float beta)

{
	float f0, f1, f2, f3;	//Auxiliary expressions
	float yt;				//Intermediate parameter of integration
	float rs;

	f0=G/V*cosf(theta)*sinf(phi);
	f1=K_beta_1*beta+f0;
	f2=-prev_rs+f1;
	f3=K_beta_2*f2;

	yt=prev_rs+DT*f3;
	rs=prev_rs+DT_2*(f3+K_beta_2*(-yt+f1));


	return rs;
}




///////////////CALCULATION_G2////////////////////////////////////////////////////////////////////////////////////////
//What calculation_g2 function does:
// - receives the angle of attack and the pitch and roll angle,
// - calculates the variable g2 useful for the controller definition.

float calculation_g2(float alpha, float theta, float phi)
{
	float g2;

	g2=G*(cosf(alpha)*cosf(theta)*cosf(phi)+sinf(alpha)*sinf(theta));

	return g2;
}



///////////////CALCULATION_G3////////////////////////////////////////////////////////////////////////////////////////
//What calculation_g3 function does:
// - receives the angle of attack and the pitch, roll and sideslip angle,
// - calculates the variable g3 useful for the controller definition.

float calculation_g3(float alpha, float theta, float phi, float beta)
{
	float cos_theta, sin_beta; 	//Auxiliary expressions
	float g3;

	cos_theta=cosf(theta);
	sin_beta=sinf(beta);

	g3=G*(cosf(beta)*cos_theta*sinf(phi)+sin_beta*cosf(alpha)*sinf(theta)-sinf(alpha)*sin_beta*cos_theta*cosf(phi));

	return g3;
}



///////////////CALCULATION_THETA////////////////////////////////////////////////////////////////////////////////////////
//What calculation_theta function does:
// - receives pitch and yaw rate and roll angle,
// - calculates the pitch angle.

float calculation_theta(float prev_theta, float q, float r, float phi)

{
	float theta;

	theta=prev_theta+DT*(q*cosf(phi)-r*sinf(phi));

	return theta;
}


///////////////CALCULATION_PHI////////////////////////////////////////////////////////////////////////////////////////
//What calculation_phi function does:
// - receives angular rates, pitch and roll angle,
// - calculates the roll angle.

float calculation_phi(float prev_phi, float theta, float p, float q, float r)

{
	float f1, tan_theta;	//Auxiliary expressions
	float yt;				//Intermediate parameter of integration
	float phi;

	tan_theta=tanf(theta);
	f1=p+q*sinf(prev_phi)*tan_theta+r*cosf(prev_phi)*tan_theta;

	yt=prev_phi+DT*f1;
	phi=prev_phi+DT_2*(f1+(p+q*sinf(yt)*tan_theta+r*cosf(yt)*tan_theta));

	return phi;
}



///////////////CALCULATION_ALPHA////////////////////////////////////////////////////////////////////////////////////////
//What calculation_alpha function does:
// - receives values inside the controller,
// - calculates angle of attack inside the controller.

float calculation_alpha(float prev_alpha, float ps, float qs, float beta, float V, float g2, float T)

{
	float f_din, f1, f2, mg2, CL, Num;		//Auxiliary expressions
	float yt;								//Intermediate parameter of integration
	float alpha;

	mg2=m*g2;
	f_din=0.5*RHO*V*V*S;
	CL=CL_alpha*prev_alpha+CL_0;
	Num=-f_din*CL-T*sinf(prev_alpha)+mg2;
	f1=m*V*cosf(beta);
	f2=qs-ps*tanf(beta)+Num/f1;

	yt=prev_alpha+DT*f2;
	alpha=prev_alpha+DT_2*(f2+ (qs-ps*tanf(beta)+ (-f_din*(CL_alpha*yt+CL_0)-T*sinf(yt)+mg2)/f1 ));

	return alpha;
}




///////////////CALCULATION_BETA////////////////////////////////////////////////////////////////////////////////////////
//What calculation_beta function does:
// - receives values inside the controller,
// - calculates the sideslip angle inside the controller.

float calculation_beta(float prev_beta, float rs, float alpha, float V, float g3, float T)

{
	float f_din, f1, f2, one_mV, mg3, T_cos_alpha;	//Auxiliary expressions
	float yt;										//Intermediate parameter of integration
	float beta;

	one_mV=1/(m*V);
	mg3=m*g3;
	f_din=0.5*RHO*V*V*S;
	f2=f_din*CY_beta;
	T_cos_alpha=T*cosf(alpha);
	f1=-rs+one_mV*(f2*prev_beta-T_cos_alpha*sinf(prev_beta)+mg3);

	yt=prev_beta+DT*f1;
	beta=prev_beta+DT_2*(f1+(-rs+one_mV*(f2*yt-T_cos_alpha*sinf(yt)+mg3)));

	return beta;
}



///////////////CALCULATION_L////////////////////////////////////////////////////////////////////////////////////////
//What calculation_L function does:
// - receives the angular rates at moment i and at i-1,
// - calculates the roll moment L.

float calculation_L(float p, float q, float r, float der_p, float der_r)

{
	float L;

	L=der_p*Ix-der_r*Ixz+q*(-p*Ixz+r*Iz)-r*q*Iy;

	return L;
}


///////////////CALCULATION_M////////////////////////////////////////////////////////////////////////////////////////
//What calculation_M function does:
// - receives the angular rates at moment i and at i-1,
// - calculates the pitch moment M.

float calculation_M(float p, float r, float der_q)

{
	float M;

	M=der_q*Iy+r*(p*Ix-r*Ixz)-p*(-p*Ixz+r*Iz);

	return M;
}




///////////////CALCULATION_N////////////////////////////////////////////////////////////////////////////////////////
//What calculation_N function does:
// - receives the angular rates at moment i and at i-1,
// - calculates the yaw moment N.

float calculation_N(float p, float q, float r, float der_p, float der_r)

{
	float N;

	N=Iz*der_r+(Iy-Ix)*p*q-Ixz*(q*r-der_p);

	return N;
}



///////////////CALCULATION_DA////////////////////////////////////////////////////////////////////////////////////////
//What calculation_da function does:
// - receives the roll moment and lateral variables,
// - calculates the aileron deflection.

int calculation_da(float L, float V, float p, float r, float beta)

{
	float da, max, m_din;
	int da_pwm;
	float p_adim, r_adim;

	p_adim=p*b/2/V;
	r_adim=r*b/2/V;
	m_din=0.5*RHO*V*V*S*b;

	//Aileron deflection
	da=(L/m_din-(Cl_beta*beta+Cl_p*p_adim+Cl_r*r_adim))/Cl_da;
	da=da+da_trim;

	//Maximum deflection with safety factor
	max=da_max*roll_sf;


	//Output

	//Saturation
	da=command_saturation(da,max);

	//Transformation in pwm
	da_pwm=command_pwm(da,-291.8,804,-170,1090.2);

	return da_pwm;
}



///////////////CALCULATION_DE////////////////////////////////////////////////////////////////////////////////////////
//What calculation_de function does:
// - receives the pitch moment and longitudinal variables,
// - calculates the elevator deflection.

int calculation_de(float M, float V, float q, float alpha)

{
	float de, max, m_din;
	int de_pwm;
	float q_adim;


	q_adim=q*c/2/V;
	m_din=0.5*RHO*V*V*S*c;

	//Elevator deflection
	de=(M/m_din-(CM0+CM_alpha*alpha+CM_q*q_adim))/CM_de;
	de=de+de_trim;

	//Maximum deflection with safety factor
	max=de_max*pitch_sf;


	//Output

	//Saturation
	de=command_saturation(de, max);

	//Transformation in pwm (from experimental calibration data)
	de_pwm=command_pwm(de,64.3,1092,256.2,994.7);


	return de_pwm;
}




///////////////CALCULATION_DR////////////////////////////////////////////////////////////////////////////////////////
//What calculation_dr function does:
// - receives the heading error and the damper gain,
// - calculates the rudder deflection.

int calculation_dr(float N, float V, float p, float r, float beta)

{
	float dr, max, m_din;
	int dr_pwm;
	float p_adim, r_adim;


	p_adim=p*b/2/V;
	r_adim=r*b/2/V;
	m_din=0.5*RHO*V*V*S*b;

	//Rudder deflection
	dr=(N/m_din-(CN_beta*beta+CN_p*p_adim+CN_r*r_adim))/CN_dr;
	dr=dr+dr_trim;

	//Maximum deflection with safety factor
	max=dr_max*yaw_sf;

	//Saturation
	dr=command_saturation(dr, max);

	//Transformation in pwm (from experimental calibration data)
	dr_pwm=command_pwm(dr,740.8,-760.5,-250.4,-828.2);


	return dr_pwm;

}




///////////////COMMAND_SATURATION////////////////////////////////////////////////////////////////////////////////////////
//What command_saturation function does:
// - receives the command deflection and the maximum/minimum command deflection,
// - calculates the allowed command deflection.

float command_saturation(float d, float max)
{

	if (d>max)
	{
		d=max;
	}
	else if (d<-max)
	{
		d=-max;
	}

	return d;
}



///////////////COMMAND_PWM////////////////////////////////////////////////////////////////////////////////////////
//What command_pwm function does:
// - receives the command deflection and curve coefficients,
// - calculates the corresponding pwm deflection.

int command_pwm(float d, float c1, float c2, float c3, float c4)
{
	int d_pwm;

	if (d<0)
	{
		d_pwm=c1*d*d+c2*d+PWM_MEAN;
	}
	else
	{
		d_pwm=c3*d*d+c4*d+PWM_MEAN;
	}

	return d_pwm;
}
