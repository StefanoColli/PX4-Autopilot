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

#include "Jerk_estimator.hpp"

JerkEstimator::JerkEstimator() :
	ModuleParams(nullptr),
	ScheduledWorkItem(MODULE_NAME, px4::wq_configurations::test1)
{
}

JerkEstimator::~JerkEstimator()
{
	perf_free(_loop_perf);
	perf_free(_loop_interval_perf);
}

bool JerkEstimator::init()
{
	// Run on fixed interval
	ScheduleOnInterval(4000_us); // 4000 us interval -> @250Hz rate

	return true;
}

void JerkEstimator::update_state()
{
	// update measured state vector (data coming from the quadcopter sensors)

	if (_local_pos_sub.updated())
	{
		vehicle_local_position_s local_pos;  //NED earth-fixed frame

		if (_local_pos_sub.copy(&local_pos))
		{
			// applied NED to ENU conversion
			if (PX4_ISFINITE(local_pos.x) && PX4_ISFINITE(local_pos.y) && PX4_ISFINITE(local_pos.z))
			{
				_measured_states_vector[0] = local_pos.y;
				_measured_states_vector[4] = local_pos.x;
				_measured_states_vector[8] = -local_pos.z;
			}

			if (PX4_ISFINITE(local_pos.vx) && PX4_ISFINITE(local_pos.vy) && PX4_ISFINITE(local_pos.vz))
			{
				_measured_states_vector[1] = local_pos.vy;
				_measured_states_vector[5] = local_pos.vx;
				_measured_states_vector[9] = -local_pos.vz;
			}
			if (PX4_ISFINITE(local_pos.ax) && PX4_ISFINITE(local_pos.ay) && PX4_ISFINITE(local_pos.az))
			{
				_measured_states_vector[2] = local_pos.ay;
				_measured_states_vector[6] = local_pos.ax;
				_measured_states_vector[10] = -local_pos.az;
			}

		}
	}

	// update yaw (psi)
	if (_vehicle_attitude_sub.updated())
	{
		vehicle_attitude_s attitude;

		if (_vehicle_attitude_sub.copy(&attitude))
		{
			if (PX4_ISFINITE(attitude.q[0]) && PX4_ISFINITE(attitude.q[1]) && PX4_ISFINITE(attitude.q[2]) && PX4_ISFINITE(attitude.q[3]))
			{
				_B_to_W = matrix::Quatf(attitude.q); //Quaternion rotation from the FRD body frame to the NED earth frame
				float psi_wrapped = -matrix::Eulerf(_B_to_W).psi(); // applied NED to ENU conversion but in [-pi;pi] range
				_measured_states_vector[12] = unwrap(_psi_prev, _psi_meas_prev, psi_wrapped); // unwrap from [-pi;pi] to [-inf;inf]

				_psi_meas_prev = psi_wrapped;
				_psi_prev = _measured_states_vector[12];
			}
		}
	}

	// update yaw rate (psi_dot)
	if (_vehicle_angular_velocity_sub.updated())
	{
		vehicle_angular_velocity_s angular_vel;	//FRD body frame

		//take yaw rate from the gyro
		if (_vehicle_angular_velocity_sub.copy(&angular_vel))
		{
			if (PX4_ISFINITE(angular_vel.xyz[0]) && PX4_ISFINITE(angular_vel.xyz[1]) && PX4_ISFINITE(angular_vel.xyz[2]))
			{
				matrix::Vector3f angular_vel_frd = matrix::Vector3f(angular_vel.xyz);
				_measured_states_vector[13] = -_B_to_W.rotateVector(angular_vel_frd)(2); // FRD to NED then NED to ENU
			}
		}
	}

	// estimate jerks from accelerations using RED differentiators
	matrix::Vector<double,1> res;
	_jerk_x_diff.evaluate(_measured_states_vector[2]);
	_jerk_x_diff.get_z(res);
	_measured_states_vector[3] = res(0);

	_jerk_y_diff.evaluate(_measured_states_vector[6]);
	_jerk_y_diff.get_z(res);
	_measured_states_vector[7] = res(0);

	_jerk_z_diff.evaluate(_measured_states_vector[10]);
	_jerk_z_diff.get_z(res);
	_measured_states_vector[11] = res(0);

}

void JerkEstimator::publish_state()
{
	states_vector_s measured_state_msg;
	// position
	measured_state_msg.pos_stp[0] = _measured_states_vector[0];
	measured_state_msg.pos_stp[1] = _measured_states_vector[4];
	measured_state_msg.pos_stp[2] = _measured_states_vector[8];
	// velocity
	measured_state_msg.vel_stp[0] = _measured_states_vector[1];
	measured_state_msg.vel_stp[1] = _measured_states_vector[5];
	measured_state_msg.vel_stp[2] = _measured_states_vector[9];
	// acceleration
	measured_state_msg.acc_stp[0] = _measured_states_vector[2];
	measured_state_msg.acc_stp[1] = _measured_states_vector[6];
	measured_state_msg.acc_stp[2] = _measured_states_vector[10];
	// jerk
	measured_state_msg.jerk_stp[0] = _measured_states_vector[3];
	measured_state_msg.jerk_stp[1] = _measured_states_vector[7];
	measured_state_msg.jerk_stp[2] = _measured_states_vector[11];
	// yaw and yaw rate
	measured_state_msg.psi = _measured_states_vector[12];
	measured_state_msg.psi_dot = _measured_states_vector[13];
	// timestamp
	measured_state_msg.timestamp = hrt_absolute_time();

	_measured_states_pub.publish(measured_state_msg);
}

void JerkEstimator::Run()
{
	if (should_exit()) {
		ScheduleClear();
		exit_and_cleanup();
		return;
	}

	perf_begin(_loop_perf);
	perf_count(_loop_interval_perf);

	update_state();
	publish_state();

	perf_end(_loop_perf);
}

int JerkEstimator::task_spawn(int argc, char *argv[])
{
	JerkEstimator *instance = new JerkEstimator();

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

int JerkEstimator::print_status()
{
	perf_print_counter(_loop_perf);
	perf_print_counter(_loop_interval_perf);
	return 0;
}

int JerkEstimator::custom_command(int argc, char *argv[])
{
	return print_usage("unknown command");
}

int JerkEstimator::print_usage(const char *reason)
{
	if (reason) {
		PX4_WARN("%s\n", reason);
	}

	PRINT_MODULE_DESCRIPTION(
		R"DESCR_STR(
### Description
The Jerk_estimator module takes the measurements from the quadcopter sensors, estimates the jerk starting from the
accelerometer data (using differentiator library) then publishes it all on the measured_states_vector topic for the
linearizer.

)DESCR_STR");

	PRINT_MODULE_USAGE_NAME("jerk_estimator", "dispatcher");
	PRINT_MODULE_USAGE_COMMAND("start");
	PRINT_MODULE_USAGE_DEFAULT_COMMANDS();

	return 0;
}

extern "C" __EXPORT int jerk_estimator_main(int argc, char *argv[])
{
	return JerkEstimator::main(argc, argv);
}

float JerkEstimator::unwrap(const float prev_angle, const float prev_meas, const float curr_meas)
{
	float diff = curr_meas - prev_meas;
	if (diff > M_PI_F)
	{
		return prev_angle + (curr_meas - M_PI_F);
	}
	else if (diff < -M_PI_F)
	{
		return prev_angle + (curr_meas + M_PI_F);
	}
	else return prev_angle + diff;
}
