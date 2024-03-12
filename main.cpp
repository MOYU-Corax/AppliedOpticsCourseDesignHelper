#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>

#define _PI_math 3.1415926
#define k9_v_D 64.06
#define k9_n 1.51630
#define lamda 0.00055

using namespace std;

typedef struct{
	double r;
	double d;
	double n_D;
	double v;
	void rayTraceInterface(double r,double d,double n_D,double v){
		this->r = r;
		this->d = d;
		this->n_D = n_D;
		this->v = v;
	};
}rayTraceInterface;
 
double AngleToRadian(double angle)
{
	return angle/180.0 * _PI_math;
}

double RadianToAngle(double radian)
{
	return radian/_PI_math * 180;
}

double RayTrace(int num,double D_obj,rayTraceInterface interface_rc[],double l1,bool type,double d_equal,double a,double U_1_alter)
{
	printf("\n\n\nRayTrace Result:\n");
	double l = l1,u = U_1_alter,n = 1,lu = D_obj/2;
	for(int i = 0;i < num;i++){
		printf("\nIndex %d\n\n",i+1);
		printf("l = %.9lf\n",l);
		printf("r = %.9lf\n",interface_rc[i].r);
		printf("l-r = %.9lf\n",l-interface_rc[i].r);
		printf("u = %.9lf\n",u);
		double i_in;
		if(interface_rc[i].r == 0){
			i_in = -u;
		}else{
			if(l == -1){
				i_in = D_obj/2 / interface_rc[i].r;
			}else{
				i_in = (l-interface_rc[i].r) / interface_rc[i].r * u;
			}
		}
		printf("i = %.9lf\n",i_in);
		
		printf("n/n` = %.9lf\n",n/interface_rc[i].n_D);
		double i_in_alter = i_in * n/interface_rc[i].n_D;
		if(i == num-1) i_in_alter = i_in * n;
		printf("i` = %.9lf\n",i_in_alter);

		double u_alter = u + i_in - i_in_alter;
		printf("u` = %.9lf\n",u_alter);
		
		double l_alter;
		if(interface_rc[i].r != 0){
			l_alter = i_in_alter * interface_rc[i].r / u_alter + interface_rc[i].r;
		}else{
			l_alter = -1;
		}
		printf("l`-r = %.9lf\n",l_alter-interface_rc[i].r);
		printf("l` = %.9lf\n",l_alter);
		
		lu = l * u;
		if(i==0) lu = D_obj/2;
		printf("lu = %.9lf\n",lu);

		double l_alter_check = lu / u_alter;
		printf("l`_check = %.9lf\n",l_alter_check);
		if(l_alter == -1) l_alter = l_alter_check;
		
		if(i==2 && type){
			l = d_equal + a;
			printf("d = %.9lf\n",l_alter - d_equal - a);
		}else{
			l = l_alter - interface_rc[i].d;
			printf("d = %.9lf\n",interface_rc[i].d);
		}

		u = u_alter;
		n = interface_rc[i].n_D;
	}
	return l;
}

double RayTraceType2(int num,double D_obj,double U_1_alter,rayTraceInterface interface_rc[],bool type,double d_equal,double a,double O_L,bool type_707)
{
	printf("\n\n\nRayTrace Result:\n");
	double L = O_L,U = 0,U_alter = U_1_alter,n = 1,lu = D_obj/2;
	for(int i=0;i<num;i++){
		printf("\nIndex %d\n\n",i+1);
		printf("L = %.9lf\n",L);
		printf("r = %.9lf\n",interface_rc[i].r);
		printf("L-r = %.9lf\n",L-interface_rc[i].r);
		double sin_i;
		printf("sinU = %.9lf\n",sin(U));
		if(interface_rc[i].r == 0){
			sin_i = -U;
		}else{
			if(i == 0){
				sin_i = D_obj/2 / interface_rc[i].r;
				if(type_707) sin_i *= 0.707;
			}else{
				sin_i = (L-interface_rc[i].r) / interface_rc[i].r * sin(U);
			}
		}
		printf("sinI = %.9lf\n",sin_i);
		
		printf("n/n` = %.9lf\n",n/interface_rc[i].n_D);
		double sin_i_alter = sin_i * n/interface_rc[i].n_D;
		printf("sinI` = %.9lf\n",sin_i_alter);
		
		double L_alter;
		printf("sinU` = %.9lf\n",sin(U_alter));
		if(interface_rc[i].r != 0){
			L_alter =  sin_i_alter * interface_rc[i].r / sin(U_alter) + interface_rc[i].r;
		}else{
			L_alter = L * sin(U);
		}
		printf("L`-r = %.9lf\n",L_alter-interface_rc[i].r);
		printf("L` = %.9lf\n",L_alter);
		printf("U = %.9lf\n",U);
		printf("I = %.9lf\n",atan(sin_i));
		printf("U+I = %.9lf\n",U+atan(sin_i));
		printf("I` = %.9lf\n",atan(sin_i_alter));
		printf("U` = %.9lf\n",U_alter);
		
		double IU_half = 0.5 * (asin(sin_i) - U);
		double IU_half_alter = 0.5 * (asin(sin_i_alter) - U_alter);
		printf("(I-U)/2 = %.9lf \n(I`-U`)/2 = %.9lf\n",IU_half,IU_half_alter);
		
		printf("LsinU = %.9lf\n",L*sin(U));
		printf("cos((I-U)/2) = %.9lf\n",cos(IU_half));
		U_alter = U + asin(sin_i) - asin(sin_i_alter);
		double PA = L * sin(U) / cos(IU_half);
		if(i == 0) PA = D_obj/2 / cos(IU_half);
		printf("PA = %.9lf\n",PA);
		
		printf("cos((I`-U`)/2) = %.9lf\n",cos(IU_half_alter));
		L_alter = PA * cos(IU_half_alter) / sin(U_alter);
		printf("sinU` = %.9lf\n",sin(U_alter));
		printf("L` = %.9lf\n",L_alter);

		if(i==2 && type){ 
			L = d_equal + a;
			printf("d = %.9lf\n", L_alter - d_equal - a);
			printf("L = %.9lf\n",L);
		}else{
			L = L_alter - interface_rc[i].d;
			printf("d = %.9lf\n",interface_rc[i].d);
			printf("L = %.9lf\n",L);
		}

		n = interface_rc[i].n_D;
		U = U_alter;
	}
	return L;
}

int Interpolate(char name[],double C_line,double P0_line,double C1,double P1,double C2,double P2)
{
	double k = (P2 - P1)/(C2 - C1);
	double b = P1 - k * C1;
	double P_0_fix = C_line * k + b;
	printf("No. %s	P_0_fix = %.9lf	P0_line = %.9lf\n",name,P_0_fix,P0_line);
	return 0;
}

double InterpolateValue(double X_line,double X1,double Y1,double X2,double Y2)
{
	double k = (Y2 - Y1)/(X2 - X1);
	double b = Y1 - k * X1;
	double Value_fix = X_line * k + b;
	return Value_fix;
}


int main()
{
	double ge,two_w;
	printf("Input ge & two_w:\n");
	scanf("%lf %lf",&ge,&two_w);
	double two_w_alter = 2.0 * atan(ge * tan(AngleToRadian(two_w)));
	double two_w_alter_use = two_w_alter * 1.05;
	printf("2w` = %.9lf	2w`_use = %.9lf\n",RadianToAngle(two_w_alter),RadianToAngle(two_w_alter_use));
	
	int lensNum;
	printf("\n\n\nInput Eyepiece Lens Num:\n");
	scanf("%d",&lensNum);
	rayTraceInterface* Eyepiece;
	Eyepiece = (rayTraceInterface*)malloc(sizeof(rayTraceInterface)*lensNum);
	memset(Eyepiece,0,sizeof(Eyepiece));
	printf("Input Each Lens Param:\n");
	printf("Format : \n r d n v\n");
	for(int i = 0;i < lensNum;i++)
	{
		printf("No.%d:\n",i);
		scanf("%lf %lf %lf %lf",&(Eyepiece[i].r),&(Eyepiece[i].d),&(Eyepiece[i].n_D),&(Eyepiece[i].v));
	}

	double Sf = RayTrace(lensNum,20,Eyepiece,-1,0,0,0,0);
	printf("\n\nSf = %.9lf\n",Sf);
	
	printf("\n\n\nInput f2:\n");
	double f2;
	scanf("%lf",&f2);
	double p_alter = Sf + f2/ge;
	printf("p` = %lf\n",p_alter);
	
	double f1 = ge * f2;
	printf("f1 = %lf\n",f1);
	
	double w = two_w/2.0;
	
	printf("\n\n\nInput D_eye:\n");
	double D_eye;
	scanf("%lf",&D_eye);
	double D_obj = D_eye * ge;
	printf("D_obj = %lf\n",D_obj);
	
	double relative_Aperture = D_obj/f1;
	printf("Relative_Aperture: %.9lf\n",relative_Aperture);

	double Df = 2.0 * tan(AngleToRadian(w)) * f1;
	printf("Df = %.9lf\n",Df);
	
	double alpha_f = 140.0/D_obj;
	printf("alpha_f = %.9lf\n",alpha_f);
	
	printf("\n\n\nInput K:\n");
	double K,n_prism = k9_n;
	scanf("%lf",&K);
	
	printf("\n\n\nInput k & a:\n");
	double k,a;
	scanf("%lf %lf",&k,&a);

	double tanAlpha = (k * D_obj - Df) / (2 * f1);
	printf("tanAlpha = %.9lf\n",tanAlpha);
	double D_prism = (n_prism * (Df + 2 * a * tanAlpha))/(n_prism - 2 * K * tanAlpha);
	printf("D_prism = %.9lf\n",D_prism);
	double D_prism_use = D_prism + 2;
	printf("D_prism_use = %.9lf\n",D_prism_use);
	
	double L = D_prism_use * K;
	double d_equal_prism = L / n_prism;
	printf("L = %.9lf\n",L);
	printf("d_equal = %.9lf\n",d_equal_prism);

	double b = f1 - a - d_equal_prism;
	printf("b = %lf\n",b);
	
	double u_m_alter = atan2(D_obj/2,f1);
	printf("u_m_alter = %.9lf\n",u_m_alter);
	double u_z = -w;
	
	double LA_prism_alter = ((n_prism * n_prism -1)/(2 * n_prism * n_prism * n_prism) * L * u_m_alter * u_m_alter);
	printf("LA_prism_alter = %.9lf\n",LA_prism_alter);
	double y_alter = Df / 2;
	double SC_prism_alter = LA_prism_alter * (AngleToRadian(u_z) / y_alter);
	printf("SC_prism_alter = %.9lf\n",SC_prism_alter);
	double delta_L_FC_prism_alter = (L/k9_v_D) * ((n_prism-1)/(n_prism*n_prism));
	printf("delta_L_FC_prism_alter = %.9lf\n",delta_L_FC_prism_alter);
	
	double LA_obj_alter = -LA_prism_alter;
	double SC_obj_alter = -SC_prism_alter;
	printf("SC_obj_alter = %.9lf\n",SC_obj_alter);
	double delta_L_FC_obj_alter = -delta_L_FC_prism_alter;
	printf("delta_L_FC_obj_alter = %.9lf\n",delta_L_FC_obj_alter);

	double J = 1 * u_m_alter * y_alter;
	printf("J = %.9lf\n",J);

	double sigma_SI = -2 * 1 * u_m_alter * u_m_alter * LA_obj_alter;
	printf("sigma_SI = %.9lf\n",sigma_SI);
	double sigma_SII = -2 * J * SC_obj_alter;
	printf("sigma_SII = %.9lf\n",sigma_SII);
	double sigma_CI = -1 * u_m_alter * u_m_alter * delta_L_FC_obj_alter;
	printf("sigma_CI = %.9lf\n",sigma_CI);
	
	double h = D_obj/2;
	double phi_obj=1/f1;
	double P_line = sigma_SI / (pow(h,4) * pow(phi_obj,3));
	printf("P_line = %.9lf\n",P_line);
	double W_line = sigma_SII / (h*h*phi_obj*phi_obj*J);
	printf("W_line = %.9lf\n",W_line);
	double C_line = sigma_CI / (h*h*phi_obj);
	printf("C_line = %.9lf\n",C_line);
	
	double P_0_line = P_line - 0.85 * pow(W_line +0.1,2);
	printf("P_0_line = %lf\nC_line = %lf\n",P_0_line,C_line);
	
	int lensGourpNum;
	printf("\n\n\nInput lensGourp Num:\n");
	scanf("%d",&lensGourpNum);
	double C1,C2;
	for(int i_1=0;i_1<lensGourpNum;i_1++){
		printf("Input Lens Gourp Para:\nFormat is\nname C1 P1 C2 P2\n");
		char name[20];
		double P1,P2;
		scanf("%s %lf %lf %lf %lf",&name,&C1,&P1,&C2,&P2);
		Interpolate(name,C_line,P_0_line,C1,P1,C2,P2);
	}

	printf("\n\n\nInput the Best Value of P_line_inter\n");
	double P_line_inter;
	scanf("%lf",&P_line_inter);
	
	double phi1_sjh,A_sjh,B_sjh,C_sjh,K_sjh,L_sjh,Q_0_sjh,P_0_sjh,W_sjh,P_sjh,n1_sjh,n2_sjh;

	printf("Input phi1 & phi2\n");
	double phi1_sjh_1,phi1_sjh_2;
	scanf("%lf %lf",&phi1_sjh_1,&phi1_sjh_2);
	phi1_sjh = InterpolateValue(C_line,C1,phi1_sjh_1,C2,phi1_sjh_2);
	printf("phi1_sjh = %.9lf\n",phi1_sjh);
	
	printf("Input A1 & A2\n");
	double A_1,A_2;
	scanf("%lf %lf",&A_1,&A_2);
	A_sjh = InterpolateValue(C_line,C1,A_1,C2,A_2);
	printf("A_sjh = %.9lf\n",A_sjh);
	
	printf("Input B1 & B2\n");
	double B_1,B_2;
	scanf("%lf %lf",&B_1,&B_2);
	B_sjh = InterpolateValue(C_line,C1,B_1,C2,B_2);
	printf("B_sjh = %.9lf\n",B_sjh);

	printf("Input C1 & C2\n");
	double C_1,C_2;
	scanf("%lf %lf",&C_1,&C_2);
	C_sjh = InterpolateValue(C_line,C1,C_1,C2,C_2);
	printf("C_sjh = %.9lf\n",C_sjh);

	printf("Input K1 & K2\n");
	double K_1,K_2;
	scanf("%lf %lf",&K_1,&K_2);
	K_sjh = InterpolateValue(C_line,C1,K_1,C2,K_2);
	printf("K_sjh = %.9lf\n",K_sjh);

	printf("Input L1 & L2\n");
	double L_1,L_2;
	scanf("%lf %lf",&L_1,&L_2);
	L_sjh = InterpolateValue(C_line,C1,L_1,C2,L_2);
	printf("L_sjh = %.9lf\n",L_sjh);

	printf("Input Q_0_1 & Q_0_2\n");
	double Q_0_1,Q_0_2;
	scanf("%lf %lf",&Q_0_1,&Q_0_2);
	Q_0_sjh = InterpolateValue(C_line,C1,Q_0_1,C2,Q_0_2);
	printf("Q_0_sjh = %.9lf\n",Q_0_sjh);

	printf("Input P_01 & P_0_2\n");
	double P_0_1,P_0_2;
	scanf("%lf %lf",&P_0_1,&P_0_2);
	P_0_sjh = InterpolateValue(C_line,C1,P_0_1,C2,P_0_2);
	printf("P_0_sjh = %.9lf\n",P_0_sjh);

	printf("Input W1 & W2\n");
	double W_1,W_2;
	scanf("%lf %lf",&W_1,&W_2);
	W_sjh = InterpolateValue(C_line,C1,W_1,C2,W_2);
	printf("W_sjh = %.9lf\n",W_sjh);

	printf("Input P1 & P2\n");
	double P_1,P_2;
	scanf("%lf %lf",&P_1,&P_2);
	P_sjh = InterpolateValue(C_line,C1,P_1,C2,P_2);
	printf("P_sjh = %.9lf\n",P_sjh);

	printf("\n\n\nInput n1 & n2:\n");
	scanf("%lf %lf",&n1_sjh,&n2_sjh);

	double Q_1_sub = Q_0_sjh - sqrt((P_line - P_line_inter)/A_sjh);
	printf("Q1_sub = %.9lf\n",Q_1_sub);
	double Q_1_add = Q_0_sjh + sqrt(((P_line - P_line_inter)/A_sjh));
	printf("Q1_add = %.9lf\n",Q_1_add);
	double Q_alter = Q_0_sjh +  (W_line - W_sjh)/K_sjh;
	printf("Q_alter = %.9lf\n",Q_alter);
	double Q_1_fin = abs(Q_1_add - Q_alter) < abs(Q_1_sub - Q_alter) ? Q_1_add:Q_1_sub;
	double Q_fin = (Q_alter+Q_1_fin)/2;
	printf("Q_fin = %.9lf\n",Q_fin);

	double rho2 = phi1_sjh + Q_fin;
	printf("rho2 = %.9lf\n",rho2);
	double rho1 = (n1_sjh/(n1_sjh-1)) * phi1_sjh + Q_fin;
	printf("rho1 = %.9lf\n",rho1);
	double rho3 = ((n2_sjh)/(n2_sjh - 1)) *phi1_sjh + Q_fin - 1/(n2_sjh - 1);
	printf("rho3 = %.9lf\n",rho3);

 	double r_1 = f1 / rho1;
	printf("r_1 = %lf\n",r_1);
 	double r_2 = f1 / rho2;
	printf("r_2 = %lf\n",r_2);
 	double r_3 = f1 / rho3;
	printf("r_3 = %lf\n",r_3);
	
	rayTraceInterface sjh[3];
	sjh[0].r = r_1;
	sjh[0].d = 0;
	sjh[0].n_D = n1_sjh;
	sjh[1].r = r_2;
	sjh[1].d = 0;
	sjh[1].n_D = n2_sjh;
	sjh[2].r = r_3;
	sjh[2].d = 0;
	sjh[2].n_D = 0;


	double l_3_alter = RayTrace(3,D_obj,sjh,-1,0,d_equal_prism,a,0);
	
	printf("\n\n\nInput u1` & u2` above this\n");
	double u_1_alter,u_2_alter;
	scanf("%lf %lf",&u_1_alter,&u_2_alter);
	
	double D_obj_use = D_obj + 2.5;
	double x_1 = r_1 > 0 ? r_1 - sqrt(pow(r_1,2)-pow(D_obj_use/2,2)):r_1 + sqrt(pow(r_1,2)-pow(D_obj_use/2,2));
	printf("x_1 = %.9lf\n",x_1);
	double x_2 = r_2 > 0 ? r_2 - sqrt(pow(r_2,2)-pow(D_obj_use/2,2)):r_2 + sqrt(pow(r_2,2)-pow(D_obj_use/2,2));
	printf("x_2 = %.9lf\n",x_2);
	double x_3 = r_3 > 0 ? r_3 - sqrt(pow(r_3,2)-pow(D_obj_use/2,2)):r_3 + sqrt(pow(r_3,2)-pow(D_obj_use/2,2));
	printf("x_3 = %.9lf\n",x_3);

	double t_1 = (D_obj_use + 3 * x_2 - 3 * x_1)/10;
	printf("t_1 = %.9lf\n",t_1);
	double t_2 = (D_obj_use + 8 * x_3 - 8 * x_2)/10;
	printf("t_2 = %.9lf\n",t_2);
	
	double d_1 = t_1 - x_2 + x_1;
	printf("d_1 = %.9lf\n",d_1);
	double d_2 = t_2 - x_3 + x_2;
	printf("d_2 = %.9lf\n",d_2);
	double d_2_use = d_2 + 1.5;
	printf("d_2_use = %.9lf\n",d_2_use);

	double h_1 = D_obj/2;
	printf("h_1 = %.9lf\n",h_1);
	double h_2 = h_1 - d_1 * u_1_alter;
	printf("h_2 = %.9lf\n",h_2);
	double h_3 = h_2 - d_2_use * u_2_alter;
	printf("h_3 = %.9lf\n",h_3);

	double r_1_fix = r_1 * h_1/h_1;
	printf("r_1_fix = %.9lf\n",r_1_fix);
	double r_2_fix = r_2 * h_2/h_1;
	printf("r_2_fix = %.9lf\n",r_2_fix);
	double r_3_fix = r_3 * h_3/h_1;
	printf("r_3_fix = %.9lf\n",r_3_fix);

	rayTraceInterface fin[5];
	memset(fin,0,sizeof(fin));
	fin[0].r = r_1_fix;
	fin[0].d = d_1;
	fin[0].n_D = n1_sjh;
	fin[1].r = r_2_fix;
	fin[1].d = d_2;
	fin[1].n_D = n2_sjh;
	fin[2].r = r_3_fix;
	fin[2].n_D = 1;
	fin[3].n_D = k9_n;
	fin[3].d = L;
	fin[4].n_D = 1;

	printf("\n\n\njzg\n");
	RayTrace(5,D_obj,fin,-1,1,d_equal_prism,a,0);

	printf("\n\n\nbyg\n");

	RayTraceType2(5,D_obj,u_1_alter,fin,1,d_equal_prism,a,-1,0);


	printf("Input PA_1:\n");
	double PA_1;
	scanf("%lf",&PA_1);
	double l_z_1 = pow(PA_1,2) / 2 / r_1;

	printf("\n\n\nzzg\n");
	RayTrace(5,D_obj,fin,l_z_1,1,d_equal_prism,a,u_1_alter);

	printf("\n\n\n0.707\n");
	RayTraceType2(5,D_obj,u_1_alter,fin,1,d_equal_prism,a,l_z_1,1);

	
	printf("\n\n\n");
	double delta_L_07_alter = 6 * lamda / pow(sin(u_m_alter),2);
	double delta_L_m_alter =  lamda / pow(sin(u_m_alter),2);
	double delta_l_FC_alter = lamda / pow(sin(u_m_alter),2);
	printf("delta_L_07_alter = %.9lf\ndelta_L_m_alter = %.9lf\ndelta_l_FC_alter = %.9lf\n",delta_L_07_alter,delta_L_m_alter,delta_l_FC_alter);
	

	system("pause");
	return 0;
}
