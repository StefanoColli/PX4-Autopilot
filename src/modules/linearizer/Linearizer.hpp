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
#include <uORB/topics/propeller_speeds.h>
#include <uORB/topics/nominal_inputs_vector.h>
#include <uORB/topics/states_vector.h>
#include <uORB/topics/control_type.h>
#include <uORB/topics/actuator_controls.h>

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
	 * @brief High Order Sliding Mode controller implementation
	 *
	 * @param nominal_states_vector nominal trajectory state vector
	 * 		[x, x_dot, x_2dot, x_3dot, y, y_dot, y_2dot, y_3dot, z, z_dot, z_2dot, z_3dot, psi, psi_dot]
	 * @param measured_states_vector actual trajectory state vector
	 * 		[x, x_dot, x_2dot, x_3dot, y, y_dot, y_2dot, y_3dot, z, z_dot, z_2dot, z_3dot, psi, psi_dot]
	 * @param v_ref nominal trajectory input vector [x_4dot, y_4dot, z_4dot, psi_2dot]
	 * @param v_out corrected input vector for trajectory tracking
	 */
	void HOSM_control(const float nominal_states_vector[14], const float measured_states_vector[14], const float v_ref[4], float v_out[4]);

	/**
	 * @brief Pole Placement controller implementation
	 *
	 * @param nominal_states_vector nominal trajectory state vector
	 * 		[x, x_dot, x_2dot, x_3dot, y, y_dot, y_2dot, y_3dot, z, z_dot, z_2dot, z_3dot, psi, psi_dot]
	 * @param measured_states_vector actual trajectory state vector
	 * 		[x, x_dot, x_2dot, x_3dot, y, y_dot, y_2dot, y_3dot, z, z_dot, z_2dot, z_3dot, psi, psi_dot]
	 * @param v_out corrected input vector [x_4dot, y_4dot, z_4dot, psi_2dot] for trajectory tracking
	 * @param k_1 gain matrix for the 3 fourth-order systems (x, y and z)
	 * @param k_2 gain matrix for the second-order system (psi)
	 */
	void PP_control(const float nominal_states_vector[14], const float measured_states_vector[14], const float k_1[4],
				const float k_2[2], float v_out[4]);

	/**
	 * @brief Compute the quadrotor model inputs: u1 = thrust force, vector u2 = moments along the 3 body axes
	 *
	 * @param input_vector vector containing the fourth time derivative of the state variables
	 * @return matrix::Vector<float, 4> quadrotor model inputs (u1, u2_vector)
	 */
	matrix::Vector<float, 4> linearize(const float input_vector[4]);

	/**
	 * @brief Compute u1 = total force in the direction of the rotors axes
	 *
	 * @param z state vector [x, x_dot, x_2dot, x_3dot, y, y_dot, y_2dot, y_3dot, z, z_dot, z_2dot, z_3dot, psi, psi_dot]
	 * @param v input vector [x_4dot, y_4dot, z_4dot, psi_2dot]
	 * @return float
	 */
	float compute_u1(const float z[14], const float v[4]);

	/**
	 * @brief Compute vector[u2, u3, u4] where u2 = roll moment, u3 = pitch moment, u4 = yaw moment
	 *
	 * @param z state vector [x, x_dot, x_2dot, x_3dot, y, y_dot, y_2dot, y_3dot, z, z_dot, z_2dot, z_3dot, psi, psi_dot]
	 * @param v input vector [x_4dot, y_4dot, z_4dot, psi_2dot]
	 * @param u1 total force in the direction of the rotors axes (computed using "compute_u1")
	 * @return matrix::Vector3f
	 */
	matrix::Vector3f compute_u2(const float z[14], const float v[4], float u1);

	/**
	 * @brief Discrete integration using backward Euler method [s = (z-1)/(Ts*z) ]
	 *
	 * @param input signal to integrate
	 * @param state	integrator state
	 * @param gain integrator gain
	 * @return float
	 */
	float discrete_integration(float input, float& state, float gain);

	/**
	 * @brief Apply the mixing rules then publish the desired propeller speeds to the actuators
	 *
	 * @param U quadrotor model inputs (u1, u2_vector)
	 */
	void mix_and_publish(matrix::Vector<float, 4>& U);

	/**
	 * @brief Math sign function
	 *
	 * @param val
	 * @return float
	 */
	static inline float signf(float val)
	{
		if (val > 0.0f)
		{
			return 1.0f;
		}
		else if (val < 0.0f)
		{
			return -1.0f;
		}
		else
		{
			return 0.0f;
		}
	}

	uint8_t _control_type = 0; // 0 = PID, 1 = HOSM with linearization, 2 = Pole Placement with linearization

	// Subscriptions
	uORB::Subscription _nominal_inputs_vector_sub{ORB_ID(nominal_inputs_vector)};
	uORB::Subscription _nominal_states_sub{ORB_ID(nominal_states_vector)};
	uORB::Subscription _measured_states_sub{ORB_ID(measured_states_vector)};
	uORB::Subscription _control_type_sub{ORB_ID(control_type)};

	// Publications
	uORB::Publication<propeller_speeds_s> _propeller_speeds_pub{ORB_ID(propeller_speeds)};
	uORB::Publication<actuator_controls_s>	_actuators_pub{ORB_ID(actuator_controls_0)};

	// Performance (perf) counters
	perf_counter_t	_loop_perf{perf_alloc(PC_ELAPSED, MODULE_NAME": cycle")};
	perf_counter_t	_loop_interval_perf{perf_alloc(PC_INTERVAL, MODULE_NAME": interval")};

	/* State vector z = [x, x_dot, x_2dot, x_3dot, y, y_dot, y_2dot, y_3dot, z, z_dot, z_2dot, z_3dot, psi, psi_dot]:
	* 	- Nominal states if the states are taken from the trajectory generator
	*	- Measured states if the states are measures taken from the quadcopter sensors
	*/
	float _nominal_states_vector[14] = {0.f};
	float _measured_states_vector[14] = {0.f};

	float _nominal_inputs[4] = {0.f}; // [x_4dot, y_4dot, z_4dot, psi_2dot] coming from the trajectory generator

	enum linearization_type {FF, FB}; 	// Feedback linearization or feedforward linearization
	linearization_type _feed_type = FF;

	static constexpr float CONSTANT_G = 9.80665f;

	const float _Ts = 0.004f; 	// sampling time for the integrators
	float _u1_2dot_int_state; 	// integrator state (u1_2dot -> u1_dot)
	float _u1_dot_int_state;	// integrator state (u1_dot -> u1)

	const float _mass = 1.5f;	 	// gazebo solo model mass
	const float inertia_data[9] = {0.015,    0,       0,
				   	0,     0.015,     0,
				   	0,       0,     0.021};
	matrix::Matrix3f _inertia = matrix::Matrix3f(inertia_data);	// inertia matrix

	const float _b = 8.54858e-06f; 		// lift coefficient
	const float _d = _b * 0.06f; 		// drag coefficient
	const float _l_1 = 0.14525; 		// distance from the rotor to the vertical symmetry axis
	const float _l_2 = 0.14745; 		// distance from the rotor to the horizontal symmetry axis

	/* Thrust to square rotor speed matrix implementing relations: (quadcopter with X configuration)
		W1^2 = 1/(4b)*U1 - 1/(4bl_2)*U2 - 1/(4bl_1)*U3 - 1/(4d)*U4
		W2^2 = 1/(4b)*U1 + 1/(4bl_2)*U2 + 1/(4bl_1)*U3 - 1/(4d)*U4
		W3^2 = 1/(4b)*U1 + 1/(4bl_2)*U2 - 1/(4bl_1)*U3 + 1/(4d)*U4
		W4^2 = 1/(4b)*U1 - 1/(4bl_2)*U2 + 1/(4bl_1)*U3 + 1/(4d)*U4
	   where:
	     -- W1..W4 are the rotational speeds of the 4 propellers,
	     -- l_1 is the distance from the rotor to the vertical symmetry axis
	     -- l_2 is the distance from the rotor to the horizontal symmetry axis
	     -- b is the rotor lift coefficient
	     -- d is the rotor drag coefficient
	*/

	float _thrust_to_speed_data[16] = { 1/(4*_b), -1/(4*_b*_l_2), -1/(4*_b*_l_1), -1/(4*_d),
			  		    1/(4*_b),  1/(4*_b*_l_2),  1/(4*_b*_l_1), -1/(4*_d),
			   		    1/(4*_b),  1/(4*_b*_l_2), -1/(4*_b*_l_1),  1/(4*_d),
					    1/(4*_b), -1/(4*_b*_l_2),  1/(4*_b*_l_1),  1/(4*_d) };

	matrix::SquareMatrix<float, 4> _thrust_to_speed{_thrust_to_speed_data};

	float _min_propeller_speed = 0;		// minimum propeller speed [rad/s]
	float _max_propeller_speed = 1300;	// maximum propeller speed [rad/s]

	// HOSM controller gains
	float _gamma[4] = {200, 200, 200, 0}; //TODO tune gains

	// PP controller gains
	float _k1[4] = { -144, -168, -73, -14};
	float _k2[2] = { -12.f, -7.f};
};
