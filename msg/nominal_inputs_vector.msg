# input vector containing highest order derivatives of the differential output vector (sigma) sent to the HOSM controller

uint64 timestamp 	# time since system start (microseconds)

float32[4] values 	# vector containing the highest order derivatives of the differential output vector (sigma)
			# values = [x_4dot, y_4dot, z_4dot, psi_2dot]
