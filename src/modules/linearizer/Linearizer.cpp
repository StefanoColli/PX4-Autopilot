/****************************************************************************
 *
 *   Copyright (c) 2023 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

#include "Linearizer.hpp"

Linearizer::Linearizer() :
	ModuleParams(nullptr),
	ScheduledWorkItem(MODULE_NAME, px4::wq_configurations::nav_and_controllers)
{
}

Linearizer::~Linearizer()
{
	perf_free(_loop_perf);
	perf_free(_loop_interval_perf);
}

bool Linearizer::init()
{
	ScheduleOnInterval(5000_us); // run every 5000 us interval -> 200 Hz rate

	// set integrators initial states
	_u1_2dot_int_state = 0;
	_u1_dot_int_state = _mass * CONSTANT_G;

	return true;
}

void Linearizer::update_state()
{
	// update input vector v
	if (_input_vector_sub.updated())
	{
		linearizer_input_s input_vector_msg;

		if (_input_vector_sub.copy(&input_vector_msg))
		{
			_input_vector[0] = input_vector_msg.values[0];
			_input_vector[1] = input_vector_msg.values[1];
			_input_vector[2] = input_vector_msg.values[2];
			_input_vector[3] = input_vector_msg.values[3];
		}
	}

	// update state vector z
	if (_feed_type == FF)	//Feedforward linearization -> state vector z from trajectory generator
	{
		if (_trajectory_vector_sub.updated())
		{
			trajectory_vector_s traj_vector;

			if (_trajectory_vector_sub.copy(&traj_vector))
			{
				_pos(0) = traj_vector.pos_stp[0];
				_pos(1) = traj_vector.pos_stp[1];
				_pos(2) = traj_vector.pos_stp[2];

				_vel(0) = traj_vector.vel_stp[0];
				_vel(1) = traj_vector.vel_stp[1];
				_vel(2) = traj_vector.vel_stp[2];

				_acc(0) = traj_vector.acc_stp[0];
				_acc(1) = traj_vector.acc_stp[1];
				_acc(2) = traj_vector.acc_stp[2];

				_jerk(0) = traj_vector.jerk_stp[0];
				_jerk(1) = traj_vector.jerk_stp[1];
				_jerk(2) = traj_vector.jerk_stp[2];

				_yaw = traj_vector.psi;
				_yaw_rate = traj_vector.psi_dot;
			}
		}

	}
	else	//Feedback linearization -> state vector z from quadcopter sensors
	{
		// update position and velocity
		if (_local_pos_sub.updated())
		{
			vehicle_local_position_s local_pos;

			if (_local_pos_sub.copy(&local_pos))
			{
				_pos(0) = local_pos.x;
				_pos(1) = local_pos.y;
				_pos(2) = local_pos.z;

				_vel(0) = local_pos.vx;
				_vel(1) = local_pos.vy;
				_vel(2) = local_pos.vz;

				// TODO this or use accelerometer data?
				//_acc(0) = local_pos.ax;
				//_acc(1) = local_pos.ay;
				//_acc(2) = local_pos.az;
			}
		}
		// update acceleration from accelerometer data
		if (_sensor_accel_sub.updated())
		{
			sensor_accel_s accel;

			if (_sensor_accel_sub.copy(&accel))
			{
				_acc(0) = accel.x;
				_acc(1) = accel.y;
				_acc(2) = accel.z;
			}
		}
		// update jerk //TODO use differentiator
		_jerk(0) = 0.001;
		_jerk(1) = 0.001;
		_jerk(2) = 0.001;
		// update yaw
		if (_vehicle_attitude_sub.updated())
		{
			vehicle_attitude_s attitude;

			if (_vehicle_attitude_sub.copy(&attitude))
			{
				matrix::Quatf quat = matrix::Quatf(attitude.q);
				double num = 2 * (quat(3) * quat(2) + quat(0) * quat(1));
    				double den = 1 - 2 * (quat(1) * quat(1) + quat(2) * quat(2));
    				_yaw = std::atan2(num, den);
			}
		}

		//update yaw rate
		if (_angular_vel_sub.updated())
		{
			vehicle_angular_velocity_s angular_velocity;

			if (_angular_vel_sub.copy(&angular_velocity))
			{
				_yaw_rate = angular_velocity.xyz[2];
			}
		}
	}
}

float Linearizer::compute_u1(float z[14], float v[4])
{
	// intemediate vectors
	matrix::Vector3f alpha = {z[2], z[6], z[10] + CONSTANT_G};
	matrix::Vector3f beta = {z[3], z[7], z[11]};
	matrix::Vector3f gamma = {v[0], v[1], -v[2]}; // - for NED frame

	float u1_2dot = _mass/alpha.norm() * (alpha.dot(gamma) + beta.norm_squared() -
		   		powf(alpha.dot(beta) / alpha.norm(), 2));

	//double integration to get u1
	float u1_dot = discrete_integration(u1_2dot, _u1_2dot_int_state, 1.0f);
	return discrete_integration(u1_dot, _u1_dot_int_state, 1.0f);
}

matrix::Vector3f Linearizer::compute_u2(float z[14], float v[4], float u1)
{
	// intermediate vectors
	matrix::Vector3f alpha = {z[2], z[6], z[10] + CONSTANT_G};
	matrix::Vector3f beta = {z[3], z[7], z[11]};

	// angles in world frame
	float eta   = sqrt(powf(alpha(2), 2) + powf((sin(z[12])*z[2] - cos(z[12])*z[6]), 2));
	float psi   = z[12];
	float phi   = atan2(sin(psi)*z[2] - cos(psi)*z[6], alpha(2));
	float theta = atan2(cos(psi)*z[2] + sin(psi)*z[6], eta);

	// body to World rotation matrix
	matrix::Vector3f R_BW_x = { cos(psi)*cos(theta)-sin(psi)*sin(phi)*sin(theta),
           		 	    sin(psi)*cos(theta)+cos(psi)*sin(phi)*sin(theta),
                                    -cos(phi)*sin(theta) };
	matrix::Vector3f R_BW_y = { -sin(psi)*cos(phi),
            		 	    cos(psi)*cos(phi),
                     	 	    sin(phi) };
	matrix::Vector3f R_BW_z = { cos(psi)*sin(theta)+sin(psi)*sin(phi)*cos(theta),
           		 	    sin(psi)*sin(theta)-cos(psi)*sin(phi)*cos(theta),
                         	    cos(phi)*cos(theta) };

	float M_rot_data[9] = { R_BW_x(0), R_BW_y(0), R_BW_z(0),
			  	R_BW_x(1), R_BW_y(1), R_BW_z(1),
			   	R_BW_x(2), R_BW_y(2), R_BW_z(2) };

	matrix::Matrix3f M_rot{M_rot_data};

	// rotational velocities in body frame
	float p = -_mass/u1 * R_BW_y.dot(beta);
	float q =  _mass/u1 * R_BW_x.dot(beta);
	float r = tan(theta) * (p + cos(phi)*sin(theta)*z[13]) + cos(phi)*cos(theta)*z[13];

	matrix::Vector3f pqr(p, q, r);

	// rotational velocities in world frame
	float dpsi   = z[13];
	float dphi   = (p + cos(phi)*sin(theta)*dpsi)/cos(theta);
	float dtheta = q - sin(phi)*dpsi;

	// Rotational accelerations in body frame (dp, dq)
	matrix::Vector3f v_sliced = {v[0], v[1], -v[2]}; //- for NED frame
	float dp = -_mass/u1 * (R_BW_y.dot(v_sliced) + 2*p*alpha.dot(beta)/alpha.norm()) + r*q;
	float dq =  _mass/u1 * (R_BW_x.dot(v_sliced) - 2*q*alpha.dot(beta)/alpha.norm()) - r*p;

	// Rotational accelerations in world frame
	float M_rot_acc_data[9] = {cos(psi), R_BW_y(0), 0,
				   sin(psi), R_BW_y(1), 0,
				      0,     R_BW_y(2), 1};
	matrix::Matrix3f M_rot_acc{M_rot_acc_data};
	matrix::Matrix3f As = M_rot.transpose() * M_rot_acc;

	matrix::Vector3f vect_a(-sin(psi)*dpsi, cos(psi)*dpsi, 0);
	matrix::Vector3f vect_b(-r, 0, p);

	matrix::Vector3f Bs = M_rot.transpose() * vect_a * dphi + vect_b*dtheta;

	float d2psi   = v[3];
	float d2phi   =  (As(1,1)*(dp-Bs(0)) + As(0,1)*(-dq+Bs(1)) + d2psi*(As(0,1)*As(1,2) - As(0,2)*As(1,1))) /
				(As(1,1)*As(0,0) - As(0,1)*As(1,0));
	float d2theta = -(As(1,0)*(dp-Bs(0)) + As(0,0)*(-dq+Bs(1)) + d2psi*(As(0,0)*As(1,2) - As(0,2)*As(1,0))) /
				(As(1,1)*As(0,0) - As(0,1)*As(1,0));

	// Rotational accelerations in body frame (dr)
	float dr = Bs(2) + As(2,0)*d2phi + As(2,1)*d2theta + As(2,2)*d2psi;

	// u_2 equation
	return _inertia*matrix::Vector3f{dp, dq, dr} + pqr.cross(_inertia * pqr);
}

float Linearizer::discrete_integration(float input, float& state, float gain)
{
	// Integration using backward Euler method

	// Compute output
	float output = state + gain * _Ts * input;
	// Update state
	state = output;

	return output;
}

void Linearizer::linearize()
{
	float z[14] = {_pos(0), _vel(0), _acc(0), _jerk(0),
		       _pos(1), _vel(1), _acc(1), _jerk(1),
		       _pos(2), _vel(2), _acc(2), _jerk(2),
		       _yaw, _yaw_rate};

	// compute thrust (u1) and roll, pitch and yaw moments (u2.xyz)
	float u1 = compute_u1(z, _input_vector);
	matrix::Vector3f u2 = compute_u2(z, _input_vector, u1);

	//PX4_WARN("u1 = %f, u2 = %f, u3 = %f, u4 = %f", (double) u1, (double) u2(0), (double) u2(1), (double) u2(2));

	// U -> OMEGA: OMEGA^2 = _thrust_to_speed * U
	matrix::Vector<float, 4> U, W;
	U(0) = u1;
	U.slice<3,1>(1,0) = u2;

	W = _thrust_to_speed * U; // apply equations described in the header to obtain the square of each propeller speed
	W = W.sqrt();

	// remap and publish propeller speeds
	propeller_speeds_s prop_speeds_msg;
	for (size_t i = 0; i < 4; i++)
	{
		//PX4_WARN("Desired propeller %d speed = %f", (int) i+1, (double) W(i));

		// remap propeller speeds from [0; max_speed] to [-1; 1]
		float midpoint = (_max_propeller_speed - _min_propeller_speed) / 2;
		W(i) = math::constrain( (W(i) - midpoint) / midpoint, -1.f, 1.f);

		//PX4_WARN("Converted propeller %d speed = %f", (int) i, (double) W(i));

		prop_speeds_msg.propeller_speed[i] = PX4_ISFINITE(W(i))? W(i) : -1.f;
	}
	prop_speeds_msg.timestamp = hrt_absolute_time();

	_propeller_speeds_pub.publish(prop_speeds_msg);
}


void Linearizer::Run()
{
	if (should_exit()) {
		ScheduleClear();
		exit_and_cleanup();
		return;
	}

	perf_begin(_loop_perf);
	perf_count(_loop_interval_perf);

	// Check if parameters have changed
	if (_parameter_update_sub.updated()) {
		// clear update
		parameter_update_s param_update;
		_parameter_update_sub.copy(&param_update);
		updateParams(); // update module parameters (in DEFINE_PARAMETERS)
	}

	update_state();
	linearize();

	perf_end(_loop_perf);
}

int Linearizer::task_spawn(int argc, char *argv[])
{
	Linearizer *instance = new Linearizer();

	if (instance) {
		_object.store(instance);
		_task_id = task_id_is_work_queue;

		if (instance->init()) {
			return PX4_OK;
		}

	} else {
		PX4_ERR("alloc failed");
	}

	delete instance;
	_object.store(nullptr);
	_task_id = -1;

	return PX4_ERROR;
}

int Linearizer::print_status()
{
	perf_print_counter(_loop_perf);
	perf_print_counter(_loop_interval_perf);
	return 0;
}

int Linearizer::custom_command(int argc, char *argv[])
{
	return print_usage("unknown command");
}

int Linearizer::print_usage(const char *reason)
{
	if (reason) {
		PX4_WARN("%s\n", reason);
	}

	PRINT_MODULE_DESCRIPTION(
		R"DESCR_STR(
### Description
Module in charge of the quadcopter linearization

)DESCR_STR");

	PRINT_MODULE_USAGE_NAME("linearizer", "controller");
	PRINT_MODULE_USAGE_COMMAND("start");
	PRINT_MODULE_USAGE_DEFAULT_COMMANDS();

	return 0;
}

extern "C" __EXPORT int linearizer_main(int argc, char *argv[])
{
	return Linearizer::main(argc, argv);
}
