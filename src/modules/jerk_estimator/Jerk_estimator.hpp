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

#include <px4_platform_common/defines.h>
#include <px4_platform_common/module.h>
#include <px4_platform_common/module_params.h>
#include <px4_platform_common/posix.h>
#include <px4_platform_common/px4_work_queue/ScheduledWorkItem.hpp>

#include <drivers/drv_hrt.h>
#include <lib/perf/perf_counter.h>
#include <differentiator/Differentiator.hpp>
#include <mathlib/math/filter/LowPassFilter2p.hpp>

#include <uORB/Publication.hpp>
#include <uORB/Subscription.hpp>
#include <uORB/topics/vehicle_local_position.h>
#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/vehicle_angular_velocity.h>
#include <uORB/topics/states_vector.h>

using namespace time_literals;

class JerkEstimator : public ModuleBase<JerkEstimator>, public ModuleParams, public px4::ScheduledWorkItem
{
public:
	JerkEstimator();
	~JerkEstimator() override;

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

	void update_state();

	void publish_state();

	static float unwrap(const float prev_angle, const float prev_meas, const float curr_meas);

	// Publication
	uORB::Publication<states_vector_s> _measured_states_pub{ORB_ID(measured_states_vector)};

	// Subscriptions
	uORB::Subscription _local_pos_sub{ORB_ID(vehicle_local_position)};
	uORB::Subscription _vehicle_attitude_sub{ORB_ID(vehicle_attitude)};
	uORB::Subscription _vehicle_angular_velocity_sub{ORB_ID(vehicle_angular_velocity)};

	// Performance (perf) counters
	perf_counter_t	_loop_perf{perf_alloc(PC_ELAPSED, MODULE_NAME": cycle")};
	perf_counter_t	_loop_interval_perf{perf_alloc(PC_INTERVAL, MODULE_NAME": interval")};

	/*Measured states vector z = [x, x_dot, x_2dot, x_3dot, y, y_dot, y_2dot, y_3dot, z, z_dot, z_2dot, z_3dot, psi, psi_dot]:
	  where the measures are taken from the quadcopter sensors
	*/
	float _measured_states_vector[14];
	matrix::Quatf _B_to_W = matrix::Quatf(0.f,0.f,0.f,0.f); 	// Quaternion rotation from the FRD body frame to the NED earth frame

	const float _Ts = 0.004f; 	// sampling time
	float _psi_prev = 0.0f; 	// yaw angle in the previous iteration
	float _psi_meas_prev = -M_PI_2_F; 	// last sampled yaw angle

	// RED Jerk estimators
	Robust_exact_differentiator _jerk_x_diff = Robust_exact_differentiator(7.0, 1, _Ts);
	Robust_exact_differentiator _jerk_y_diff = Robust_exact_differentiator(7.0, 1, _Ts);
	Robust_exact_differentiator _jerk_z_diff = Robust_exact_differentiator(7.0, 1, _Ts);

};
