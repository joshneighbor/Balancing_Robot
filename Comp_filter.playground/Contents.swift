//Josh Neighbor -  Complimentary Filter
//Library
#include <robotics_cape.h> // From Strawson Designs
//Declarations
#define DEFAULT_SAMPLE_RATE 20
unsigned short gyro_fsr;
int sample_rate, print_yet, first_run;
float last_gyro, thetaYAccel, thetaZAccel, alpha_HP, alpha_LP, gyro_X, theta_LP, theta_HP, theta, h_inv, tau, thetaA_Ynaught, thetaA_Znaught, sample_constanter;
mpudata_t mpu;
float thetaG_Xnaught;
float accel_current, last_accel, alpha;
//The comp filter is called by IMU interrupt function, a thread spawned in main which polls
//the IMU for updated values of theta dot from the gyro and the theta x from accelerometer
//this function also prints the sampled values at the rate prescribed in the sample_constanter
int comp_filter(){
    memset(&mpu, 0, sizeof(mpudata_t)); //makes sure it's clean before starting
    if (mpu9150_read(&mpu) == 0) {
        //on first run, the theta values are initialized along with filter constants
        if(first_run==0){
            thetaA_Ynaught = mpu.calibratedAccel[VEC3_Y];
            thetaA_Znaught = mpu.calibratedAccel[VEC3_Z];
            thetaG_Xnaught = mpu.rawGyro[VEC3_X];
            //initialize theta_LP at current theta
            theta = theta_LP;
            //computes the constants so they are not caluclated repeatedly
            alpha_HP = (2*tau*h_inv/(2*tau+h_inv));
            alpha_LP = (h_inv/(2*tau+h_inv));
            alpha = ((2*tau-h_inv)/(2*tau+h_inv));
            //initialize the lowpass and highpass filters
            theta_LP = alpha_LP*atan2(thetaA_Znaught,thetaA_Ynaught);
            theta_HP = theta_LP;
            first_run = 1;
        }
        //collect gyro and accelerometer readings
        thetaYAccel = mpu.calibratedAccel[VEC3_Y];
        thetaZAccel = mpu.calibratedAccel[VEC3_Z];
        //converts the gyro reading to radians per second
        gyro_X = (mpu.rawGyro[VEC3_X]-thetaG_Xnaught)*gyro_fsr*DEG_TO_RAD/(32768*PI);
    }
    //calculates the arctan once so it may be stored later without recalculating
    accel_current = (atan2(thetaZAccel,thetaYAccel));
    
    //this is where the filter values are calculated
    theta_LP = alpha_LP*(atan2(thetaZAccel,thetaYAccel) + last_accel)+alpha*theta_LP;
    theta_HP = alpha_HP*(gyro_X + last_gyro)+alpha*theta_HP;
    //sums the two filters to make a complementary filter
    theta = theta_HP + theta_LP;
    
    //store current gyro/accel readings to be the previous time step inputs
    //in the next iteration
    last_gyro = gyro_X*h_inv;
    last_accel = accel_current;
    
    //Counter for print loop rate
    print_yet +=1;
    if (print_yet >= sample_constanter)
    {
        printf("\n");
        printf("%0.2f, ",	theta_HP);
        printf("%0.2f, ",	theta_LP);
        printf("%0.2f, ", 	theta);
        //this compares the comp filter estimates with the IMU's built-in
        //filter estimates
        printf("%0.2f", -mpu.fusedEuler[VEC3_X]+atan2(thetaA_Znaught,thetaA_Ynaught));
        print_yet = 0;
    }
    return 0;
}
//First initializes the cape, then reads-in the sample rate argument, if present
//if no specified sampling rate, default is used.
int main(int argc, char *argv[]){
    initialize_cape();
    if (argc==1){
        sample_rate = DEFAULT_SAMPLE_RATE;
        tau = .5;
    }
    else{
        sample_rate = atoi(argv[1]);
        tau = .2;
        //imposes max and min samle rate requirements
        if((sample_rate>MAX_SAMPLE_RATE)||(sample_rate<MIN_SAMPLE_RATE)){
            printf("sample rate should be between %d and %d\n", MIN_SAMPLE_RATE,MAX_SAMPLE_RATE);
            return -1;
        }
    }
    //initializes boolean values
    print_yet = 0;
    first_run = 0;
    h_inv = 1.0/sample_rate;
    //this value ensures the printer statement in the complementary filter
    //prints out at a constant rate, regardless of the sampling rate
    sample_constanter = sample_rate/10;
    fflush(stdout);
    
    //For initializing the IMU
    signed char orientation[9] = ORIENTATION_UPRIGHT;
    //signed char orientation[9] = ORIENTATION_FLAT;
    initialize_imu(sample_rate, orientation);
    last_gyro = 0;
    // read the gyro full-scale range
    mpu_get_gyro_fsr(&gyro_fsr);
    // start the filter
    set_imu_interrupt_func(&comp_filter); 
    
    //Program idles until terminated by user
    while(get_state()!=EXITING){
        sleep(1);
    }
    cleanup_cape();
    return 0;}
