
/*
NathanielBot-4.0
Microbots  from the island by the sea
Modified code by : Christopher M . Coballes
Phillippines 
*/

#include <wiringPi.h>
#include <ros/ros.h>
#include <ros/time.h>
#include <geometry_msgs/Vector3Stamped.h>
#include <tf/transform_broadcaster.h>
#include <geometry_msgs/Vector3.h>
#include <geometry_msgs/Twist.h>
#include <nav_msgs/Odometry.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <robot_specs.h>

//Motor Shield headers
#include <Wire.h>
#include <Adafruit_MotorShield.h>
#include "utility/Adafruit_PWMServoDriver.h"

#define encodPinA1      3     // encoder A pin
#define encodPinB1      8     // encoder B pin
#define encodPinA2      2
#define encodPinB2      7
#define LOOPTIME        100   // PID loop time(ms)
#define SMOOTH      10

#define sign(x) (x > 0) - (x < 0)

// Create the motor shield object with the default I2C address
Adafruit_MotorShield AFMS = Adafruit_MotorShield(); 
// Select which 'port' M1, M2, M3 or M4. 
Adafruit_DCMotor *motor1 = AFMS.getMotor(1);
Adafruit_DCMotor *motor2 = AFMS.getMotor(2);

// Base Controller

unsigned long lastMilli = 0;       // loop timing 
unsigned long lastMilliPub = 0;
double rpm_req1 = 0;
double rpm_req2 = 0;
double rpm_act1 = 0;
double rpm_act2 = 0;
double rpm_req1_smoothed = 0;
double rpm_req2_smoothed = 0;
int direction1 = FORWARD;
int direction2 = FORWARD;
int prev_direction1 = RELEASE;
int prev_direction2 = RELEASE;
int PWM_val1 = 0;
int PWM_val2 = 0;
volatile long count1 = 0;          // rev counter
volatile long count2 = 0;
long countAnt1 = 0;
long countAnt2 = 0;
float Kp =   0.5;
float Kd =   0;
float Ki =   0;
//ros::NodeHandle nh;

void handle_cmd( const geometry_msgs::Twist& cmd_msg) {
  double x = cmd_msg.linear.x;
  double z = cmd_msg.angular.z;
  if (z == 0) {     // go straight
    // convert m/s to rpm
    rpm_req1 = x*60/(pi*wheel_diameter);
    rpm_req2 = rpm_req1;
  }
  else if (x == 0) {
    // convert rad/s to rpm
    rpm_req2 = z*track_width*60/(wheel_diameter*pi*2);
    rpm_req1 = -rpm_req2;
  }
  else {
    rpm_req1 = x*60/(pi*wheel_diameter)-z*track_width*60/(wheel_diameter*pi*2);
    rpm_req2 = x*60/(pi*wheel_diameter)+z*track_width*60/(wheel_diameter*pi*2);
  }
}




double rpm_act1 = 0.0;
double rpm_act2 = 0.0;
double rpm_req1 = 0.0;
double rpm_req2 = 0.0;
double gyro_x = 0.0;
double gyro_y = 0.0;
double gyro_z = 0.0;
double rpm_dt = 0.0;
double x_pos = 0.0;
double y_pos = 0.0;
double theta = 0.0;
ros::Time current_time;
ros::Time rpm_time(0.0);
ros::Time last_time(0.0);

void handle_rpm( const geometry_msgs::Vector3Stamped& rpm) {
  rpm_act1 = rpm.vector.x;
  rpm_act2 = rpm.vector.y;
  rpm_dt = rpm.vector.z;
  rpm_time = rpm.header.stamp;
}

void handle_gyro( const geometry_msgs::Vector3& gyro) {
  gyro_x = gyro.x;
  gyro_y = gyro.y;
  gyro_z = gyro.z;
}

int main(int argc, char** argv){
ros::init(argc, argv, "base_controller");
  ros::NodeHandle n;
  ros::NodeHandle nh_private_("~");
  ros::Subscriber<geometry_msgs::Twist> sub("cmd_vel", handle_cmd);
  
  ros::Subscriber sub = n.subscribe("rpm", 50, handle_rpm);
  
  ros::Subscriber gyro_sub = n.subscribe("gyro", 50, handle_gyro);
  
  geometry_msgs::Vector3Stamped rpm_msg;
  ros::Publisher rpm_pub("rpm", &rpm_msg);
  ros::Publisher odom_pub = n.advertise<nav_msgs::Odometry>("odom", 50);
  
  tf::TransformBroadcaster broadcaster;
  ros::Time current_time;
  ros::Time last_time;

  double rate = 10.0;
  double linear_scale_positive = 1.0;
  double linear_scale_negative = 1.0;
  double angular_scale_positive = 1.0;
  double angular_scale_negative = 1.0;
  double angular_scale_accel = 1.0;
  double acc_theta = 0.0;
  double acc_x = 0.0;
  double acc_max_theta = 0.0;
  double acc_max_x = 0.0;
  double alpha = 0.0;
  bool publish_tf = true;
  bool use_imu = false;
  double dt = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dth_odom = 0.0;
  double dth_gyro = 0.0;
  double dth = 0.0;
  double dth_prev = 0.0;
  double dth_curr = 0.0;
  double dxy_prev = 0.0;
  double dxy_ave = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vth = 0.0;
  char base_link[] = "/base_link";
  char odom[] = "/odom";
  ros::Duration d(1.0);
  // this is the launh file
  nh_private_.getParam("publish_rate", rate);
  nh_private_.getParam("publish_tf", publish_tf);
  nh_private_.getParam("linear_scale_positive", linear_scale_positive);
  nh_private_.getParam("linear_scale_negative", linear_scale_negative);
  nh_private_.getParam("angular_scale_positive", angular_scale_positive);
  nh_private_.getParam("angular_scale_negative", angular_scale_negative);
  nh_private_.getParam("angular_scale_accel", angular_scale_accel);
  nh_private_.getParam("alpha", alpha);
  nh_private_.getParam("use_imu", use_imu);

  ros::Rate r(rate);
  initNathan();

 while(n.ok()){
  nh.spinOnce();
  unsigned long time = millis();
  if(time-lastMilli>= LOOPTIME)   {      // enter tmed loop
    getMotorData(time-lastMilli);
    PWM_val1 = updatePid(1, PWM_val1, rpm_req1, rpm_act1);
    PWM_val2 = updatePid(2, PWM_val2, rpm_req2, rpm_act2);

    if(PWM_val1 > 0) direction1 = FORWARD;
    else if(PWM_val1 < 0) direction1 = BACKWARD;
    if (rpm_req1 == 0) direction1 = RELEASE;
    if(PWM_val2 > 0) direction2 = FORWARD;
    else if(PWM_val2 < 0) direction2 = BACKWARD;
    if (rpm_req2 == 0) direction2 = RELEASE;
    motor1->run(direction1);
    motor2->run(direction2);

    motor1->setSpeed(abs(PWM_val1));
    motor2->setSpeed(abs(PWM_val2));
    
    publishRPM(time-lastMilli);
    lastMilli = time;
  }
  if(time-lastMilliPub >= LOOPTIME) {
  //  publishRPM(time-lastMilliPub);
    lastMilliPub = time;
  }
}
} //main  ROS

void getMotorData(unsigned long time)  {
 rpm_act1 = double((count1-countAnt1)*60*1000)/double(time*encoder_pulse*gear_ratio);
 rpm_act2 = double((count2-countAnt2)*60*1000)/double(time*encoder_pulse*gear_ratio);
 countAnt1 = count1;
 countAnt2 = count2;
}

int updatePid(int id, int command, double targetValue, double currentValue) {
  double pidTerm = 0;                            // PID correction
  double error = 0;
  double new_pwm = 0;
  double new_cmd = 0;
  static double last_error1 = 0;
  static double last_error2 = 0;
  static double int_error1 = 0;
  static double int_error2 = 0;
  
  error = targetValue-currentValue;
  if (id == 1) {
    int_error1 += error;
    pidTerm = Kp*error + Kd*(error-last_error1) + Ki*int_error1;
    last_error1 = error;
  }
  else {
    int_error2 += error;
    pidTerm = Kp*error + Kd*(error-last_error2) + Ki*int_error2;
    last_error2 = error;
  }
  new_pwm = constrain(double(command)*MAX_RPM/4095.0 + pidTerm, -MAX_RPM, MAX_RPM);
  new_cmd = 4095.0*new_pwm/MAX_RPM;
  return int(new_cmd);
}

void publishRPM(unsigned long time) {
  rpm_msg.header.stamp = nh.now();
  rpm_msg.vector.x = rpm_act1;
  rpm_msg.vector.y = rpm_act2;
  rpm_msg.vector.z = double(time)/1000;
  rpm_pub.publish(&rpm_msg);
  nh.spinOnce();
}

void encoder1() {
  if (digitalRead(encodPinA1) == digitalRead(encodPinB1)) count1++;
  else count1--;
}
void encoder2() {
  if (digitalRead(encodPinA2) == digitalRead(encodPinB2)) count2--;
  else count2++;
}

void publishOdom()
{

    // ros::topic::waitForMessage<geometry_msgs::Vector3Stamped>("rpm", n, d);
    current_time = ros::Time::now();
    dt = rpm_dt;
    dxy_ave = (rpm_act1+rpm_act2)*dt*wheel_diameter*pi/(60*2);
    dth_odom = (rpm_act2-rpm_act1)*dt*wheel_diameter*pi/(60*track_width);

    if (use_imu) dth_gyro = dt*gyro_z;
    dth = alpha*dth_odom + (1-alpha)*dth_gyro;

    if (dth > 0) dth *= angular_scale_positive;
    if (dth < 0) dth *= angular_scale_negative;
    if (dxy_ave > 0) dxy_ave *= linear_scale_positive;
    if (dxy_ave > 0) dxy_ave *= linear_scale_negative;

    dx = cos(dth) * dxy_ave;
    dy = -sin(dth) * dxy_ave;

    x_pos += (cos(theta) * dx - sin(theta) * dy);
    y_pos += (sin(theta) * dx + cos(theta) * dy);
    theta += dth;

    if(theta >= two_pi) theta -= two_pi;
    if(theta <= -two_pi) theta += two_pi;

    geometry_msgs::Quaternion odom_quat = tf::createQuaternionMsgFromYaw(theta);

    if(publish_tf) {
      geometry_msgs::TransformStamped t;
      t.header.frame_id = odom;
      t.child_frame_id = base_link;
      t.transform.translation.x = x_pos;
      t.transform.translation.y = y_pos;
      t.transform.translation.z = 0.0;
      t.transform.rotation = odom_quat;
      t.header.stamp = current_time;

      broadcaster.sendTransform(t);
    }

    nav_msgs::Odometry odom_msg;
    odom_msg.header.stamp = current_time;
    odom_msg.header.frame_id = odom;
    odom_msg.pose.pose.position.x = x_pos;
    odom_msg.pose.pose.position.y = y_pos;
    odom_msg.pose.pose.position.z = 0.0;
    odom_msg.pose.pose.orientation = odom_quat;
    if (rpm_act1 == 0 && rpm_act2 == 0){
      odom_msg.pose.covariance[0] = 1e-9;
      odom_msg.pose.covariance[7] = 1e-3;
      odom_msg.pose.covariance[8] = 1e-9;
      odom_msg.pose.covariance[14] = 1e6;
      odom_msg.pose.covariance[21] = 1e6;
      odom_msg.pose.covariance[28] = 1e6;
      odom_msg.pose.covariance[35] = 1e-9;
      odom_msg.twist.covariance[0] = 1e-9;
      odom_msg.twist.covariance[7] = 1e-3;
      odom_msg.twist.covariance[8] = 1e-9;
      odom_msg.twist.covariance[14] = 1e6;
      odom_msg.twist.covariance[21] = 1e6;
      odom_msg.twist.covariance[28] = 1e6;
      odom_msg.twist.covariance[35] = 1e-9;
    }
    else{
      odom_msg.pose.covariance[0] = 1e-3;
      odom_msg.pose.covariance[7] = 1e-3;
      odom_msg.pose.covariance[8] = 0.0;
      odom_msg.pose.covariance[14] = 1e6;
      odom_msg.pose.covariance[21] = 1e6;
      odom_msg.pose.covariance[28] = 1e6;
      odom_msg.pose.covariance[35] = 1e3;
      odom_msg.twist.covariance[0] = 1e-3;
      odom_msg.twist.covariance[7] = 1e-3;
      odom_msg.twist.covariance[8] = 0.0;
      odom_msg.twist.covariance[14] = 1e6;
      odom_msg.twist.covariance[21] = 1e6;
      odom_msg.twist.covariance[28] = 1e6;
      odom_msg.twist.covariance[35] = 1e3;
    }
    vx = (dt == 0)?  0 : dxy_ave/dt;
    vth = (dt == 0)? 0 : dth/dt;
    odom_msg.child_frame_id = base_link;
    odom_msg.twist.twist.linear.x = vx;
    odom_msg.twist.twist.linear.y = 0.0;
    odom_msg.twist.twist.angular.z = dth;

    odom_pub.publish(odom_msg);
    last_time = current_time;
     ros::spinOnce();
    r.sleep();
  }	
  
  void initNathan()
  {

 AFMS.begin();  // create with the default frequency 1.6KHz
 count1 = 0;
 count2 = 0;
 countAnt1 = 0;
 countAnt2 = 0;
 rpm_req1 = 0;
 rpm_req2 = 0;
 rpm_act1 = 0;
 rpm_act2 = 0;
 PWM_val1 = 0;
 PWM_val2 = 0;
 nh.initNode();
 nh.getHardware()->setBaud(57600);
 nh.subscribe(sub);
 nh.advertise(rpm_pub);
  
 pinMode(encodPinA1, INPUT); 
 pinMode(encodPinB1, INPUT); 
 digitalWrite(encodPinA1, HIGH);                // turn on pullup resistor
 digitalWrite(encodPinB1, HIGH);
 attachInterrupt(1, encoder1, RISING);

 pinMode(encodPinA2, INPUT); 
 pinMode(encodPinB2, INPUT); 
 digitalWrite(encodPinA2, HIGH);                // turn on pullup resistor
 digitalWrite(encodPinB2, HIGH);
 attachInterrupt(0, encoder2, RISING);
 motor1->setSpeed(0);
 motor2->setSpeed(0);
 motor1->run(FORWARD);
 motor1->run(RELEASE);
 motor2->run(FORWARD);
 motor2->run(RELEASE);
}
