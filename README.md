# SelfDrivingCar2-MPC
# A MPC controller for vehicle actuation and steering in a virtual track

## The Controller: Brief Overview

The cost function includes both orientation (OE) and cross-track (CTE) error terms, as well as regularizatyon on the steering and actuation inputs, as well as on their transition. Furthermore, deviation from a specific speed value is penalized to reduce gidder while the controller is in charge of actuation. The time horizon is 8 seconds and is divided in 33 steps, so that the interstep interval is greater than the system's latency time, which is roughly 110 ms (100 ms simulator latency, plus a little something to account for network lag).

Provisions are taken to account for the rather large latency time by forwarding the state 110 ms in the future and then using it as a starting state in the MPC solver. More details on the transition model, cost function formulation, choice of time horizon and step duration, as well as counter-latency time provisions are given [in this report]().  

### Performance and Tuning
I have to admit I am not very fond of the parameter tuning, but it does the job nonetheless. The controller can safely steer the vehicle throughout the track either with fixed throttle (in the range of 0.2), or with the actuation value provided by the solver. ** In any way, the vehicle stays on track, but there is clearly large margin for improvement to say the least**. I found that different weighting of the quadratic error terms in the cost function can make a huge difference in the outcome. All the above of course, provided that the solver does its job well. Further more, tweaking the time horizon in conjunction with the number if steps can also dramtically affect the behavior of the controller. Finally, another important is latency time. We know that it is above 100 ms, but as to how long it takes for the network to propage messages is unknown; I used 10 ms because anything over 150 ms (also below 100 ms) makes an apparent difference to the worse.    

## Compiling and Executing
Provided that ipopt is alrteady installed, from the root directory of the project,
```
mkdir build
cd build
cmake ..
make
./mpc
```
Then run the simulator in auto-driving mode. The simulator can be dowloaded from [here](https://github.com/udacity/self-driving-car-sim/releases).

