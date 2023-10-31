# Drone state received from ground control station. See also custom Mavlink message called "NOMINAL_STATES_VECTOR"
# The topic states_vector describe the position, velocity, acceleration and jerk of the current target setpoint

uint64 timestamp 	# time since system start (microseconds)

float32[3] pos_stp 	# Target position setpoint x,y,z [m]
float32[3] vel_stp 	# Target velocity setpoint x,y,z [m/s]
float32[3] acc_stp 	# Target acceleration setpoint x,y,z [m/s^2]
float32[3] jerk_stp 	# Target jerk setpoint x,y,z [m/s^3]
float32 psi	 	# Target orientation setpoint [rad]
float32 psi_dot 	# Target angular rate setpoint [rad/s]

# TOPICS nominal_states_vector measured_states_vector
