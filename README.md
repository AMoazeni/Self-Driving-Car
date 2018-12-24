# Introduction

<br></br>
Take me to the [MATLAB Simulation Code](https://github.com/AMoazeni/Self-Driving-Car/blob/master/Code/Autonomous%20Drifting%20Simulation.m) for the Self Driving Car!


<br></br>
Development of autonomous vehicles has accelerated in the past decade due to advances in computing speed, sensor technology, and popular interest. This article explores the Software Architecture for the Self-Driving Car shown below. The controller uses a Model Predictive Control (MPC) algorithm to anticipate the car's future position, knowing the car's Vehicle Dynamics equations and measured position (current state).


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/01%20-%20Car%20Sensors.png" width=75% alt="Car-Sensor"></div>


<br></br>
The algorithm should work on any [Drive-By-Wire](https://en.wikipedia.org/wiki/Drive_by_wire) car which has electronically controlled steering wheel, gas, and brake pedals. You can read sensor values (Camera, Radar, Lidar, GPS) and control the car directly from the car's [CAN Bus](https://en.wikipedia.org/wiki/CAN_bus), check with the manufacturer to confirm.



<br></br>
<br></br>

# Model Predictive Control (MPC) Algorithm

<br></br>
This architecture lets you control the vehicle acceleration, brake, and steering using Model Predictive Control (MPC). The software architecture shown above has been tested on the Self-Driving Hyundai Sonata shown in the above picture. It has been successfully tested in highway speeds and city driving conditions.


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/02%20-%20Control%20System.png" alt="Control-System"></div>


<br></br>
At each sampling time step, beginning at the current state, an open loop optimal control problem is solved over a finite horizon. For each consecutive time step, time step a new optimal control problem based on new measurements is solved over a shifted horizon.


<br></br>
The optimal solution relies on a dynamic model of the process with respects to input constraints, output constraints, and minimizing a performance index (cost). The cost equation for this model is the simple distance formula, this keeps the car at the center of the lane by minimizing cost error which is the car's center position coordinate minus lane center coordinate.


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/03%20-%20MPC%20Algorithm.png" width=75% alt="MPC"></div>



<br></br>
<br></br>

# Vehicle Dynamics

<br></br>
Here is a description of a car's Vehicle Dynamics. You can find the constants for your car with some simple research or taking measurements. The simulated vehicle in this project tracks a circular trajectory, but the algorithm can follow any trajectory. Your real-world trajectory comes from lane detection which is generated by your camera sensor's Computer Vision (CV) system.


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/04%20-%20Vehicle%20Dynamics.png" width=75% alt="Vehicle-Dynamics"></div>


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/05%20-%20Vehicle%20Equations.png" width=60% alt="Vehicle-Equations"></div>


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/06%20-%20Vehicle%20Constants.png" width=30% alt="Vehicle-Constants"></div>


<br></br>
The only input considered is the vehicle steering angle (steering wheel) and acceleration (gas and brake pedal). Some constraints are implemented to smooth out driving commands like a limited steering angle range, and limited rate of change for steering angle which reduces jerk in steering. 


<br></br>
The Vehicle Model equations need to be [linearized](https://apmonitor.com/pdc/index.php/Main/ModelLinearization) which is like finding the next position of a point given its initial position and slope. Linearized equations are MUCH easier to calculate from a computer's perspective.


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/07%20-%20Vehicle%20Linear.png" width=75% alt="Vehicle-Linear"></div>



<br></br>
<br></br>

# Bicycle Model

<br></br>
Bicycle models are not as computationally expensive as car models, and also less prone to errors due to their simplicity. It's highly recommended to use a Bicycle Model for the real Self-Driving Car. Check out this research paper which compares [Vehicle versus Bicycle Models for Autonomous Vehicles](https://github.com/AMoazeni/Self-Driving-Car/blob/master/Research%20Papers/01%20-%20Vehicle%20Models%20for%20Autonomous%20Driving.pdf). Also read this research paper on [Autonomous Drifting](https://github.com/AMoazeni/Self-Driving-Car/blob/master/Research%20Papers/02%20-%20Autonomous%20Vehicle%20Drifting%20with%20MPC.pdf).


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/08%20-%20Bicycle%20Model.png" width=75% alt="Bicycle-Model"></div>


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/09%20-%20Bicycle%20Equations.png" width=50% alt="Bicycle-Equations"></div>


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/10%20-%20Bicycle%20Symbols.png" alt="Bicycle-Symbols"></div>



<br></br>
<br></br>

# Autonomous Drifting!

<br></br>
The tire model equation is not necessary for normal driving. It's only included for fun to simulate drifting. It forces tire friction to saturate (no traction) and get the car to drift! Remove the tire equations and friction saturation constraints to simulate regular driving.


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/11%20-%20Drifting.png" width=75% alt="Drifting"></div>


<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/12%20-%20Tire%20Model.png" width=75% alt="Tire-Model"></div>


<br></br>
<br></br>

# Results
<br></br>
The following images are made in MATLAB, they show the lane (Red lines), the car and steering wheel angle (Black box), and the covered trajectory (Blue line).


### Normal Driving

<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/13%20-%20Result%20Car.png" width=50% alt="Result-Car"></div>


### Drifting

<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/14%20-%20Result%20Drifting.gif" width=50% alt="Drifting-Result"></div>




<br></br>
<br></br>

# Robot Operating System (ROS)

<br></br>
The [Robot Operating System (ROS)](http://www.ros.org/) is an open source platform that's quickly becoming an industry standard in the Robotics field. ROS let's you quickly prototype and reuse code through their robust publisher-subscriber ([pub-sub](https://en.wikipedia.org/wiki/Publish%E2%80%93subscribe_pattern)) architecture. You can also use different languages to create your nodes which is a huge bonus for modular design.

<br></br>
ROS applies a Soft Real-Time system to your code which is fine for single projects but becomes problematic when implemented at scale. There are workarounds for getting Hard Real-Time performance, but that's a topic for another post. Follow along with these [ROS Tutorials](http://wiki.ros.org/ROS/Tutorials) to get you started, section 1.1 steps 1-13 are the most important tutorials. The following Self-Driving Car visualization was made in ROS using the Rviz package.

<br></br>
<div align="center"><img src="https://raw.githubusercontent.com/AMoazeni/Self-Driving-Car/master/Jupyter%20Notebook/Images/15%20-%20ROS.png" alt="ROS"></div>



<br></br>
<br></br>

# Code

<br></br>
This repository only contains the [MATLAB simulation](https://github.com/AMoazeni/Self-Driving-Car/blob/master/Code/Autonomous%20Drifting%20Simulation.m) for this project. However in the real world, this project was implemented on the Self-Driving Car using Robot Operating System (ROS). Find the MATLAB simulation in the 'Code' folder of this repository.


<br></br>
```shell
$ git clone https://github.com/AMoazeni/Self-Driving-Car.git
$ cd Self-Driving-Car
```

<br></br>
<br></br>
<br></br>
<br></br>
