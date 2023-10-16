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

#pragma once

#include <parameters/param.h>

#include <px4_platform_common/defines.h>
#include <px4_platform_common/module.h>
#include <px4_platform_common/module_params.h>
#include <px4_platform_common/posix.h>
#include <px4_platform_common/px4_work_queue/ScheduledWorkItem.hpp>

#include <drivers/drv_hrt.h>
#include <lib/perf/perf_counter.h>
#include <mathlib/mathlib.h>

#include <uORB/Publication.hpp>
#include <uORB/Subscription.hpp>
#include <uORB/SubscriptionCallback.hpp>
#include <uORB/topics/parameter_update.h>
#include <uORB/topics/sensor_accel.h>
#include <uORB/topics/vehicle_local_position.h>
#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/vehicle_angular_velocity.h>
#include <uORB/topics/trajectory_vector.h>
#include <uORB/topics/control_type.h>
#include <uORB/topics/propeller_speeds.h>
#include <uORB/topics/linearizer_input.h>

using namespace time_literals;

class Linearizer : public ModuleBase<Linearizer>, public ModuleParams, public px4::ScheduledWorkItem
{
public:
	Linearizer();
	~Linearizer() override;

	/** @see ModuleBase */
	static int task_spawn(int argc, char *argv[]);

	/** @see ModuleBase */
	static int custom_command(int argc, char *argv[]);

	/** @see ModuleBase */
	static int print_usage(const char *reason = nullptr);

	bool init();

	int print_status() override;

private:
	void Run() override;

	/**
	 * @brief Update linearizer state (inputs) checking if new messages have been published
	 *
	 */
	void update_state();

	/**
	 * @brief Compute vector U (linearization), then translate it to desired propeller speeds.
	 *  Finally publish the speeds to the "propeller_speeds" topic.
	 *
	 */
	void linearize();

	/**
	 * @brief Compute u1 = total force in the direction of the rotors axes
	 *
	 * @param z state vector [x, x_dot, x_2dot, x_3dot, y, y_dot, y_2dot, y_3dot, z, z_dot, z_2dot, z_3dot, psi, psi_dot]
	 * @param v input vector [x_4dot, y_4dot, z_4dot, psi_2dot]
	 * @return float
	 */
	float compute_u1(float z[14], float v[4]);

	/**
	 * @brief Compute vector[u2, u3, u4] where u2 = roll moment, u3 = pitch moment, u4 = yaw moment
	 *
	 * @param z state vector [x, x_dot, x_2dot, x_3dot, y, y_dot, y_2dot, y_3dot, z, z_dot, z_2dot, z_3dot, psi, psi_dot]
	 * @param v input vector [x_4dot, y_4dot, z_4dot, psi_2dot]
	 * @param u1 total force in the direction of the rotors axes (computed using "compute_u1")
	 * @return matrix::Vector3f
	 */
	matrix::Vector3f compute_u2(float z[14], float v[4], float u1);

	/**
	 * @brief Discrete integration using backward Euler method [s = (z-1)/(Ts*z) ]
	 *
	 * @param input signal to integrate
	 * @param state	integrator state
	 * @param gain integrator gain
	 * @return float
	 */
	float discrete_integration(float input, float& state, float gain);

	// Subscriptions
	uORB::Subscription _local_pos_sub{ORB_ID(vehicle_local_position)};
	uORB::Subscription _sensor_accel_sub{ORB_ID(sensor_accel)};
	uORB::Subscription _vehicle_attitude_sub{ORB_ID(vehicle_attitude)};
	uORB::Subscription _angular_vel_sub{ORB_ID(vehicle_angular_velocity)};
	uORB::Subscription _trajectory_vector_sub{ORB_ID(trajectory_vector)};
	uORB::Subscription _parameter_update_sub{ORB_ID(parameter_update)};
	uORB::Subscription _input_vector_sub{ORB_ID(linearizer_input)};

	// Publications
	uORB::Publication<propeller_speeds_s> _propeller_speeds_pub{ORB_ID(propeller_speeds)};

	// Performance (perf) counters
	perf_counter_t	_loop_perf{perf_alloc(PC_ELAPSED, MODULE_NAME": cycle")};
	perf_counter_t	_loop_interval_perf{perf_alloc(PC_INTERVAL, MODULE_NAME": interval")};

	matrix::Vector3f _pos; // position in NED frame
	matrix::Vector3f _vel; // velocity in NED frame
	matrix::Vector3f _acc; // acceleration in NED frame
	matrix::Vector3f _jerk; // jerk in NED frame
	float _input_vector[4] = {0.0, 0.0, 0.0, 0.0}; // [x_4dot, y_4dot, z_4dot, psi_2dot]
	float _yaw;
	float _yaw_rate;

	enum linearization_type {FF, FB}; 	// Feedback linearization or feedforward linearization
	linearization_type _feed_type = FF;

	static constexpr float CONSTANT_G = 9.80665f;

	const float _Ts = 0.005; 	// sampling time
	float _u1_2dot_int_state; 	// integrator state (u1_2dot -> u1_dot)
	float _u1_dot_int_state;	// integrator state (u1_dot -> u1)

	float _mass = 1.714 ;		 // gazebo solo model mass
	float inertia_data[9] = {0.011,     0,       0,
				   0,  	  0.015,     0,
				   0,       0,     0.021};
	matrix::Matrix3f _inertia = matrix::Matrix3f(inertia_data);	// inertia matrix

	float _b = 8.54858e-06; 	// lift coefficient
	float _d = 1.140e-7; 		// drag coefficient
	float _l_1 = 0.14525; 		// distance from the rotor to the vertical symmetry axis
	float _l_2 = 0.14745; 		// distance from the rotor to the horizontal symmetry axis

	/* Thrust to square rotor speed matrix implementing relations: (quadcopter with X configuration)
		W1^2 = 1/(4b)*U1 - 1/(4bl_1)*U2 - 1/(4bl_2)*U3 - 1/(4d)*U4
		W2^2 = 1/(4b)*U1 + 1/(4bl_1)*U2 + 1/(4bl_2)*U3 - 1/(4d)*U4
		W3^2 = 1/(4b)*U1 + 1/(4bl_1)*U2 - 1/(4bl_2)*U3 + 1/(4d)*U4
		W4^2 = 1/(4b)*U1 - 1/(4bl_1)*U2 + 1/(4bl_2)*U3 + 1/(4d)*U4
	   where:
	     -- W1..W4 are the rotational speeds of the 4 propellers,
	     -- l_1 is the distance from the rotor to the vertical symmetry axis
	     -- l_2 is the distance from the rotor to the horizontal symmetry axis
	     -- b is the rotor lift coefficient
	     -- d is the rotor drag coefficient
	*/

	float _thrust_to_speed_data[16] = { 1/(4*_b), -1/(4*_b*_l_1), -1/(4*_b*_l_2), -1/(4*_d),
			  		    1/(4*_b),  1/(4*_b*_l_1),  1/(4*_b*_l_2), -1/(4*_d),
			   		    1/(4*_b),  1/(4*_b*_l_1), -1/(4*_b*_l_2),  1/(4*_d),
					    1/(4*_b), -1/(4*_b*_l_1),  1/(4*_b*_l_2),  1/(4*_d) };

	matrix::SquareMatrix<float, 4> _thrust_to_speed{_thrust_to_speed_data};

	float _min_propeller_speed = 0;		// minimum propeller speed [rad/s]
	float _max_propeller_speed = 1500;	// maximum propeller speed [rad/s]
};
