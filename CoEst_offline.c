#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <malloc.h>

#define tforder 	2	// order of transfer function
#define coeffnum 	5	// =tforder *2 +1
#define win_size 	20	// window size for parameter identification
#define time_data	1	// (=1/freq_data) sampling time for the data receiving (second)
#define time_ID 	5	// sampling time for the parameter id (second)
#define time_OB 	1	// sampling time for the observer (second)
#define N 			25	// max size of the matrix
#define N_1 		8000	// max size of the matrix
#define ratedcurr	1.36	// 1 C-rate

// claim the function used in the program
void print_PI_result();	// Print results of parameter identification
void print_OB_result();	// Print results of observer
void getParams();		// Parameter identification
void printdata();		// print the data in J matrix
void dataprocess();		// build the J matrix
void estimateSOC(int i);	// Observer
void lookup();			// SOC-VOC lookup table
void readdata();		// read data from Arduino
void readdata_1(int i);	// read data from pre-defined matrix
void readfile();		// read data - from a file
void savedata(int time, FILE * file2);
						// save data into a file
void coulombcount();	// Coulomb counting
float determinant(float [][25], float);
void cofactor(float [][25], float);
void transpose(float [][25], float [][25], float);

struct timespec beg,end;	// count the time based on system timer

// define global variables
// --- temp flag ---
int tmp0 = 0;
int tmp1 = 0;
int flag = 0;
// --- crucial ---
float curr_threshold  = 0.5;	// for regional-awareness parameter id
float volt_cutoff_up  = 4.2;	// upper cut-off voltage
int   freq_data = 10;	// frequency of the data receiving (Hz)
// --- not crucial ---
long counter = 0;
int status = -1;	// -1: discharging, 1: charging
float volt = 0.0;	// terminal votlage
float curr = 0.0;	// load current
float curr_old = 0.0;	// load current (past)
float curr_win1[N_1];
float volt_win1[N_1];
float out_Y[win_size];	// array for parameter identification
float mat_J[win_size][coeffnum];	// array for parameter identification
float mat_JT[coeffnum][win_size];	// array for parameter identification
float mat_JTJ[N][N];
					// Temp matrix of parameter identificaiton
float mat_JTJI[coeffnum][coeffnum];
					// Inverse matrix
float mat_JTJIJT[coeffnum][win_size];
					// JTJIJT
int det_a;			// det A index (0: cannot inverse)
float R0   = 0.09;	// R_0
float R1   = 0.005;	// R_1
float C1   = 400;	// C_1
float b1Qr = 0.0002;	// b_1 / Q_R
float para[coeffnum];		// wrap the identified parameters
float coeff_tfd[coeffnum];
					// coefficient of the identified discrete-time transfer function 
float Ad[tforder][tforder];	
					// A matrix in observer
float Bd[tforder];	// B matrix in observer
float Cd[tforder];	// C matrix in observer
float pole_c[tforder];	// polse assignments for the observer gain
float L_d[tforder];	// observer gain
float err = 0.0;	// observer error
float estVoc = 3.0;	// estimated V_OC  (x_1 in observer)
float estVrc1= 0.0;	// estimated V_RC1 (x_2 in observer)
float estVt  = 0.0;	// estimated V_T (y observer)
float estSOC = 0.0;	// estimated SOC (V_OC -> lookup table -> SOC)
float estSOC_nxt = 0.0;	// estimated SOC output
float SOC_ref = 0.0;	// SOC reference
float estSOC_max = 0.0;	// maximum SOC of the battery cell
float pole_d[tforder];
float coeff1, coeff2, detA;
static float table_SOC[100] = {\
	0.010900, 0.020570, 0.030250, 0.039920, 0.049590,\
	0.059270, 0.068940, 0.078620, 0.088290, 0.097960,\
	0.107640, 0.117310, 0.126980, 0.136660, 0.146330,\
	0.156010, 0.165680, 0.175350, 0.185030, 0.194700,\
	0.204370, 0.214050, 0.223720, 0.233400, 0.243070,\
	0.252740, 0.262420, 0.272090, 0.281760, 0.291440,\
	0.301110, 0.310790, 0.320460, 0.330130, 0.339810,\
	0.349480, 0.359150, 0.368830, 0.378500, 0.388180,\
	0.397850, 0.407520, 0.417200, 0.426870, 0.436540,\
	0.446220, 0.455890, 0.465570, 0.475240, 0.484910,\
	0.494590, 0.504260, 0.513930, 0.523610, 0.533280,\
	0.542960, 0.552630, 0.562300, 0.571980, 0.581650,\
	0.591320, 0.601000, 0.610670, 0.620350, 0.630020,\
	0.639690, 0.649370, 0.659040, 0.668710, 0.678390,\
	0.688060, 0.697740, 0.707410, 0.717080, 0.726760,\
	0.736430, 0.746100, 0.755780, 0.765450, 0.775130,\
	0.784800, 0.794470, 0.804150, 0.813820, 0.823490,\
	0.833170, 0.842840, 0.852520, 0.862190, 0.871860,\
	0.881540, 0.891210, 0.900880, 0.910560, 0.920230,\
	0.929910, 0.939580, 0.949250, 0.958930, 0.968600};
static float table_Voc[100] = {\
	3.325000, 3.344000, 3.364000, 3.383000, 3.403000,\
	3.422000, 3.441000, 3.461000, 3.480000, 3.500000,\
	3.518000, 3.530000, 3.542000, 3.555000, 3.567000,\
	3.579000, 3.591000, 3.603000, 3.615000, 3.627000,\
	3.639000, 3.649000, 3.660000, 3.670000, 3.681000,\
	3.691000, 3.702000, 3.713000, 3.723000, 3.734000,\
	3.742000, 3.746000, 3.749000, 3.752000, 3.756000,\
	3.759000, 3.763000, 3.766000, 3.770000, 3.773000,\
	3.776000, 3.777000, 3.779000, 3.781000, 3.782000,\
	3.784000, 3.785000, 3.787000, 3.788000, 3.790000,\
	3.793000, 3.796000, 3.800000, 3.803000, 3.807000,\
	3.810000, 3.814000, 3.817000, 3.821000, 3.824000,\
	3.829000, 3.836000, 3.842000, 3.849000, 3.855000,\
	3.862000, 3.868000, 3.875000, 3.881000, 3.888000,\
	3.895000, 3.903000, 3.911000, 3.919000, 3.926000,\
	3.934000, 3.942000, 3.950000, 3.958000, 3.965000,\
	3.975000, 3.984000, 3.994000, 4.003000, 4.013000,\
	4.022000, 4.032000, 4.041000, 4.051000, 4.060000,\
	4.071000, 4.081000, 4.091000, 4.101000, 4.111000,\
	4.121000, 4.132000, 4.142000, 4.152000, 4.162000};


// Main function
int main(void){
	// --- Read data from files and to array ---
	readfile();			// read data from file
	// --- Open files for data saving ---
	FILE * file2 = fopen("out_est_SOC.txt", "w");
	// --- Count the clock time ---
	clock_gettime(CLOCK_MONOTONIC_RAW,&beg);
	// Co-Estimation algorithm
	for (int Loop = 0; Loop <= N_1-1; Loop++){
		// --- Read data from the array ---
		readdata_1(Loop);		// read data - test
		// --- Coulomb counting: SOC reference ---
		coulombcount();			// Coulomb counting
		// --- Parameter Identification ---
		if (Loop % time_ID == 0){
			dataprocess();		// build the J matrix
			getParams();		// parameter identification
		}
		// --- Observer ---
		if (Loop % time_OB == 0){
			estimateSOC(Loop);	// observer
			print_OB_result();	// print result
		}
		// --- Save data into files ---
		savedata(Loop,file2);
								// save data into a file
	}
	// --- Close the saved data ---
	fclose(file2);
	return 0;
}

// read data - from a file
void readfile(){
	FILE * myfile_curr;
	FILE * myfile_volt;
	myfile_curr = fopen("current_pisa.txt","r");
	myfile_volt = fopen("voltage_pisa.txt","r");
	for(int i = 0; i <= N_1-1; i++){
		fscanf(myfile_curr, "%f", &curr_win1[i]);
		fscanf(myfile_volt, "%f", &volt_win1[i]);
	}
	fclose(myfile_curr);
	fclose(myfile_volt);
}

// read data - test - pre-defined matrix
void readdata_1(int i1){
	curr = curr_win1[i1];
	volt = volt_win1[i1];
}

// Coulomb counting
void coulombcount(){
	SOC_ref += curr * time_data / 5165.66;
}

// build the J matrix
void dataprocess(){
	// moves the data in the moving window
	for (int i = 0; i <= win_size-2; i++){
		out_Y[i] = out_Y[i+1];
		for (int j = 0; j <= coeffnum-1 ; j++)
			mat_J[i][j] = mat_J[i+1][j];
	}
	// assign the new data
	out_Y[win_size-1]    =  volt;
	mat_J[win_size-1][0] = -1 * out_Y[win_size-2];
	mat_J[win_size-1][1] = -1 * out_Y[win_size-3];
	mat_J[win_size-1][2] = curr;
	mat_J[win_size-1][3] = mat_J[win_size-2][2];
	mat_J[win_size-1][4] = mat_J[win_size-3][2];
}

// Parameter identification
void getParams(){
	// parameter identification
	// regional-awareness
	if ((mat_J[10][3]-mat_J[9][3]) >= curr_threshold){
		float detA = 0.0;
		printf("-- Performing the parameter identificaiton. --\n");

		// --- theta(k)=(J(k)^T .* J(k))^-1 .* J(k) * y(k) ---
		// --- Step 1: Transpose the matrix ---
		for (int i = 0; i <= win_size-1; i++){
			for (int j = 0; j <= coeffnum-1 ; j++)
				mat_JT[j][i] = mat_J[i][j];
		}
		// --- Step 2: mat_JTJ = J^T .* J ---
		for (int i = 0; i <= coeffnum-1 ; i++){
			for (int k = 0; k <= coeffnum-1 ; k++){
				for (int j = 0; j <= win_size-1; j++){
					mat_JTJ[i][k] += mat_JT[i][j] * mat_J[j][k];
				}
			}
		}
		// --- Step 3: mat_JTJI = (J^T .* J)^-1 ---
		detA = determinant(mat_JTJ, coeffnum);
		if (detA == 0)
			printf("\nInverse of Entered Matrix is not possible\n");
		else
			cofactor(mat_JTJ, coeffnum);
		// --- Step 4: mat_JTJIJ = (J^T.*J)^-1*J ---
		for (int i = 0; i <= coeffnum-1 ; i++){
			for (int k = 0; k <= win_size-1 ; k++){
				for (int j = 0; j <= coeffnum-1; j++){
					mat_JTJIJT[i][k] += mat_JTJI[i][j] * mat_JT[j][k];
				}
			}
		}
		// --- Step 5: para1 = mat_JTJIJY = (J^T.*J)^-1*J*Y ---
		for (int i = 0; i <= coeffnum-1; i++){
			for (int j = 0; j <= win_size-1 ; j++)
				coeff_tfd[i] = mat_JTJIJT[i][j] * out_Y[j];
		}
		// --- Step 6: Pull the coefficients out ---
		float a1 = coeff_tfd[0];
		float a2 = coeff_tfd[1];
		float c0 = coeff_tfd[2];
		float c1 = coeff_tfd[3];
		float c2 = coeff_tfd[4];
		// --- Step 7: convert from discrete tf to continuous tf ---
		float n0 = (c0-c1+c2)*time_ID*time_ID/4.0;
		float n1 = (c0-c2)*time_ID;
		float n2 = c0+c1+c2;
		float d0 = (1-a1+a2)*time_ID*time_ID/4.0;
		float d1 = (1-a2)*time_ID;
		float d2 = (1+a1+a2);
		// --- Step 8: Normalizing the first coefficient in the denominator to 1 ---
		d1 = d1/d0;
		d2 = d2/d0;
		n0 = n0/d0;
		n1 = n1/d0;
		n2 = n2/d0;
		//  Step 9 - Convert to the Pole-Zero form //
		float detden, pole_PI[tforder];
		// solving the a*x^2 + b*x + c
		// x = [ -b +- sqrt(b^2-4*a*c) ] / (2 * a)
		detden = d1 * d1 - 4 * d2;
		if(detden >= 0){
			pole_PI[0] = (-d1 + sqrt(detden))/2;
			pole_PI[1] = (-d1 - sqrt(detden))/2;
		}
		// --- Step 10: Pull the parameters out ---
		R0 = n0;
		n1 -= n0 * d1;
		n2 -= n0 * d2;
		printf("%f\t%f\n", pole_PI[0], pole_PI[1]);
		float RC = 1.0 / pole_PI[0];
		C1 = (pole_PI[0] - pole_PI[1]) / (n1 * pole_PI[0] + n2);
		R1 = RC / C1;
		b1Qr = (n1 * pole_PI[1] + n2) / (pole_PI[1] - pole_PI[0]);
		// --- Step 11: data processing ---
		R0 = fabs(R0);
		R1 = fabs(R1);
		C1 = fabs(C1);
		b1Qr = fabs(b1Qr);
		
		para[0] = R0;
		para[1] = R1;
		para[2] = C1;
		para[3] = b1Qr;

		// print the identified parameters
		printf("the identified para using ARX:\n");
		printf("R_0\t\tR_1\t\tC_1\t\tRC_1\t\tb1/QR\n");
		for (int i = 0; i<= coeffnum-1; i++){
			printf("%.10f\t", para[i]);
		}	printf("\n");
	}
}

// SOC-VOC lookup table
void lookup(void){
	int i = 0; // index of look-up table
	float slope = 0.0;
	// find the index
	while (estVoc > table_Voc[i] && i < 100){
		i++;
	}
	// find the slope
	slope  = (table_SOC[i] - table_SOC[i-1])/(table_Voc[i] - table_Voc[i-1]);
	estSOC = table_SOC[i-1] + slope * (estVoc - table_Voc[i-1]);

	// saturation
	if(estVoc >= table_Voc[99]){
		estSOC = table_SOC[99];
	}else if(estVoc <= table_Voc[0]){
		estSOC = table_SOC[0];
	}
}

// Observer
void estimateSOC(int i1){
	// Output the SOC_cell (x^(k))
	lookup();
	float Ad_alpha[tforder][tforder];
	float CA[tforder];
	float N_matI[tforder][tforder];
	float N_mat[tforder][tforder];
	float AA[tforder][tforder];
	float a1,a2,alpha1,alpha2;
	
	// update the SOC of next time stamp
	// update the state-space matrix
	Ad[0][0] = 1;
	Ad[0][1] = 0;
	Ad[1][0] = 0;
	Ad[1][1] = exp(-time_OB / R1 / C1);	// e^(-Ts/RC)
	Bd[0] = time_OB * b1Qr + time_OB / C1;
	Bd[1] = time_OB * b1Qr + (1-Ad[1][1]) * R1;
	Cd[0] = 1;
	Cd[1] = 1;

	// Pole placement using Ackermann's formula:
	// 		L_d =  phi(Ad) * N_matI * [0 0....1]'
	// Step 1: Assign the poles: pole_c
	pole_c[0] = -0.05;
	pole_c[1] = -3 / R1 / C1;
	// Step 2: Convert to Discrete time
	pole_d[0] = exp(pole_c[0] * time_OB);
	pole_d[1] = exp(pole_c[1] * time_OB);
	// Step 3: Calculate coefficients for equation phi(Ad)
	alpha1 = -(pole_d[0] + pole_d[1]);	
	alpha2 = (pole_d[0])*(pole_d[1]);
	// Step 4: Calculate Cd*Ad = CA
	CA[0] = Cd[0]*Ad[0][0] + Cd[1]*Ad[1][0];
	CA[1] = Cd[0]*Ad[0][1] + Cd[1]*Ad[1][1];
	// Step 5: N_mat = [Cd Cd*Ad ...]'
	N_mat[0][0] = Cd[0];
	N_mat[0][1] = Cd[1];
	N_mat[1][0] = CA[0];
	N_mat[1][1] = CA[1];  
	// Step 6: N_matI = N_mat^-1
	float den;
	float detN_mat;
	detN_mat = N_mat[0][0] * N_mat[1][1] - N_mat[0][1] * N_mat[1][0];
	if(detN_mat == 0)
		printf("Inverse not Possible in observer computation");
	else{
		den = 1/detN_mat;
		N_matI[0][0] =  N_mat[1][1] * den;
		N_matI[0][1] = -N_mat[0][1] * den;
		N_matI[1][0] = -N_mat[1][0] * den;
		N_matI[1][1] =  N_mat[0][0] * den;
	}
	// Step 7: Calculate Ad^2 =  A * A
	for (int i = 0; i < tforder; i++){
		for (int j = 0; j < tforder; j++){
			AA[i][j] = 0;
			for (int k = 0; k < tforder; k++)
				AA[i][j] += Ad[i][k] * Ad[k][j];
		}
	}
	// Step 8: alpha1 * Ad = Ad_alpha
	for (int i = 0; i < tforder; i++){
		for (int j = 0; j < tforder; j++){
			Ad_alpha[i][j] = alpha1 * Ad[i][j];
		}
	}
	// Step 9: AA = Ad^2 + alpha1 * Ad + alpha2 * I
	for (int i = 0; i < tforder; i++){
		for (int j = 0; j < tforder; j++){
			if(i == j)
				AA[i][j] = AA[i][j] + Ad_alpha[i][j] + alpha2;
			else
				AA[i][j] = AA[i][j] + Ad_alpha[i][j] + 0;			
		}
	}	
	// Step 10: L1 = phi(A) * [C CA]^-1  = AA * N_matI
	float L1[tforder][tforder];
	for (int i = 0; i < tforder; i++){
		for (int j = 0; j < tforder; j++){
			L1[i][j] = 0;
			for (int k = 0; k < tforder; k++)
				L1[i][j] += AA[i][k]*N_matI[k][j];
		}
	}
	// Step 11: observer gain L_d =  phi(Ad)*N_matI*[0 0 ... 1]'
	L_d[0] = L1[0][1];
	L_d[1] = L1[1][1];

	// --- y^(k)   = C * x^(k) + D * u(k) ---
	estVt = Cd[0]*estVoc + Cd[1]*estVrc1 + R0*curr;
	err = volt - estVt;
	// --- xhat = [SOC, V_RC1, V_RC2] ---
	// --- Update x^(k+1) ---
	// --- x^(k+1) = A * x^(k) + B * u(k) + L * (y - V_T) ---
	estVoc  = Ad[0][0] * estVoc  + Bd[0] * curr;
	estVrc1 = Ad[1][1] * estVrc1 + Bd[1] * curr;

	estVoc  += L_d[0] * err * fabs(curr / ratedcurr);
	estVrc1 += L_d[1] * err * fabs(curr / ratedcurr);
}

/*For calculating Determinant of the Matrix */
float determinant(float a[25][25], float k){
	float s = 1, det = 0, b[25][25];
	int i, j, m, n, c;
	if (k == 1){
    	return (a[0][0]);
	}else{
		det = 0;
		for (c = 0; c < k; c++){
			m = 0;
			n = 0;
			for (i = 0;i < k; i++){
				for (j = 0 ;j < k; j++){
					b[i][j] = 0;
					if (i != 0 && j != c){
						b[m][n] = a[i][j];
						if (n < (k - 2))
							n++;
						else{
							n = 0;
							m++;
						}
					}
				}
			}
		det = det + s * (a[0][c] * determinant(b, k - 1));
		s = -1 * s;
		}
	}
	return (det);
}
void cofactor(float num[25][25], float f){
	float b[25][25], fac[25][25];
	int p, q, m, n, i, j;
	for (q = 0;q < f; q++){
		for (p = 0;p < f; p++){
			m = 0;
			n = 0;
			for (i = 0;i < f; i++){
				for (j = 0;j < f; j++){
					if (i != q && j != p){
						b[m][n] = num[i][j];
						if (n < (f - 2))
							n++;
						else{
							n = 0;
							m++;
						}
					}
				}
			} fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
		}
	} transpose(num, fac, f);
}

/*Finding transpose of matrix*/
void transpose(float num[25][25], float fac[25][25], float r){
	int i, j;
	float b[25][25], d;
	for (i = 0;i < r; i++){
		for (j = 0;j < r; j++){
			b[i][j] = fac[j][i];
		}
	} d = determinant(num, r);
	for (i = 0;i < r; i++){
		for (j = 0;j < r; j++){
			mat_JTJI[i][j] = b[i][j] / d;
		}
	}
}

// print the data in J matrix
void printdata(){
	// check the Y matrix
	printf("the data in the Y matrix is:\n");
	for (int i = 0; i <= win_size-1; i++){
		printf("\t%.2f", out_Y[i]);
	}	printf("\n\n");

	// check the J matrix
	printf("the data in the J matrix is:\n");
	for (int i = 0; i <= win_size-1; i++){
		for (int j = 0; j <= coeffnum-1 ; j++){
			printf("\t%.2f", mat_J[i][j]);
		}	printf("\n");
	}	printf("\n");
}

// Print results of parameter identification
void print_PI_result(){
	// check the JT matrix
	printf("the data in the JT matrix is:\n");
	for (int i = 0; i <= coeffnum-1; i++){
		for (int j = 0; j <= win_size-1 ; j++){
			printf("\t%.2f", mat_JT[i][j]);
		}	printf("\n");
	}	printf("\n");

	// check the JTJ matrix
	printf("the JTJ matrix is:\n");
	for (int i = 0; i <= coeffnum-1; i++){
		for (int j = 0; j <= coeffnum-1 ; j++)
			printf("\t%.2f", mat_JTJ[i][j]);
		printf("\n");
	}	printf("\n");

	// check the JTJI matrix
	printf("The JTJI matrix is:\n");
	for (int i = 0; i <= coeffnum-1; i++){
		for (int j = 0; j <= coeffnum-1 ; j++)
			printf("\t%.2f", mat_JTJI[i][j]);
		printf("\n");
	}	printf("\n");
	// check the JTJI matrix
	printf("The JTJIJT matrix is:\n");
	for (int i = 0; i <= coeffnum-1; i++){
		for (int j = 0; j <= win_size-1 ; j++)
			printf("\t%.2f", mat_JTJIJT[i][j]);
		printf("\n");
	}	printf("\n");
	// check the identifed coefficietns of tf_d
	printf("the coefficient of the discrete time tf:\n");
	for (int i = 0; i <= coeffnum-1; i++){
		printf("%.2f\t", coeff_tfd[i]);
	}	printf("\n");
}

// Print results of observer
void print_OB_result(){
	// Print the estimated battery cell VOC
	printf("Estd cell VOC is %.2f Volt.\n", estVoc);
	// Print the estimated battery cell SOC
	printf("Estd cell SOC is %.2f percent.\n", estSOC*100);
}

// save data into a file
void savedata(int time, FILE * file2){
	// count the time of each cycles
	clock_gettime(CLOCK_MONOTONIC_RAW,&end);
	float time_interval = (end.tv_sec - beg.tv_sec)*1000 + (end.tv_nsec - beg.tv_nsec)/1000000;	// nanosecond
	time_interval /= 1000;	// milisecond
	char array1[50] ="time (ms) \testSoC (%) \tSoC_ref (%) \tSoC_err (%)";
	// Save data
	if(time == 0){
		fprintf(file2, "%s",array1);
	}
	else{
	fprintf(file2, "\n%f\t%f\t%f\t%f", time_interval, estSOC * 100, SOC_ref * 100, (estSOC - SOC_ref) * 100);
	}
}









