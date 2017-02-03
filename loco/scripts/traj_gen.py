"""
Generate open-loop control input trajectories for system identification tests.

Parameters for all tests: dt [ms], total time of test
Two control inputs: velocity (m/s) and steering (radians)

Tests:
* constant forward velocity with ramping steering ("ramp")
    * params: forward velocity, steering angle velocity, max angle
* fixed zero steering angle, short step inputs of throttle ("step")
    * params: gap, step_length, velocity

Format of txt control input files:
dt [timestep of control inputs]
u [velocity(t=0)] [steering(t=0)]
u [velocity(t=dt)] [steering(t=dt)]
....
"""

dt = 0.05  # [s] Publish frequency is 1/dt. Keep it under 60Hz?
dt_int = int(dt*1000)
T = 5000   # [ms] Total length of trajectory

experiment = "step"  # "ramp" or "step"

# Ramp experiment parameters
vx_ramp = 1.0    # Forward velocity, [m/s]
max_steering = 0.75

# Step experiment parameters
vx_step = 1.0  # [m/s]
gap = 1000  # [ms]
step_length = 300  # [ms]

outfile = "generated.txt"
file = open(outfile, "w")

if (experiment == "ramp"):
    for t in range(0, T+dt_int, dt_int):
        steer = (max_steering)*(t)/T
        command = "u " + str(vx_ramp) + " " + str(steer) + " " + str(dt)
        print(command)
        file.write(command + '\n')

if (experiment == "step"):
    tot_period = gap + step_length
    for t in range(0, T+dt_int, dt_int):
        if ((t % tot_period) < step_length):
            command = "u " + str(vx_step) + " " + str(0.0) + " " + str(dt)
        else:
            command = "u " + str(0.0) + " " + str(0.0) + " " + str(dt)
        print(command)
        file.write(command + '\n')

file.write("end")

file.close()
