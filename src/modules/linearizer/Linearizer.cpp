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
	ScheduledWorkItem(MODULE_NAME, px4::wq_configurations::test2)
{
}

Linearizer::~Linearizer()
{
	perf_free(_loop_perf);
	perf_free(_loop_interval_perf);
}

bool Linearizer::init()
{
	ScheduleOnInterval(4000_us); // run every 4000 us interval -> @250Hz rate

	// set integrators initial states
	_u1_2dot_int_state = 0;
	_u1_dot_int_state = _mass * CONSTANT_G;

	return true;
}

void Linearizer::update_state()
{
	// check control approach: PID or HOSM with linearization or PP with linearization
	control_type_s control_type_msg;
	if (_control_type_sub.update(&control_type_msg))
	{
		_control_type = control_type_msg.control_type;
	}

	// update nominal input vector (HOSM controller input)
	// NED -> ENU conversion
	if (_nominal_inputs_vector_sub.updated())
	{
		nominal_inputs_vector_s nominal_inputs_vector_msg;

		if (_nominal_inputs_vector_sub.copy(&nominal_inputs_vector_msg))
		{
			_nominal_inputs[0] = nominal_inputs_vector_msg.values[1];
			_nominal_inputs[1] = nominal_inputs_vector_msg.values[0];
			_nominal_inputs[2] =- nominal_inputs_vector_msg.values[2];
			_nominal_inputs[3] = nominal_inputs_vector_msg.values[3];
		}
	}

	// update nominal state vector (data coming from the trajectory generator)
	// NED -> ENU conversion
	if (_nominal_states_sub.updated())
	{
		states_vector_s nominal_states_msg;

		if (_nominal_states_sub.copy(&nominal_states_msg))
		{
			_nominal_states_vector[0] = nominal_states_msg.pos_stp[1];
			_nominal_states_vector[1] = nominal_states_msg.vel_stp[1];
			_nominal_states_vector[2] = nominal_states_msg.acc_stp[1];
			_nominal_states_vector[3] = nominal_states_msg.jerk_stp[1];

			_nominal_states_vector[4] = nominal_states_msg.pos_stp[0];
			_nominal_states_vector[5] = nominal_states_msg.vel_stp[0];
			_nominal_states_vector[6] = nominal_states_msg.acc_stp[0];
			_nominal_states_vector[7] = nominal_states_msg.jerk_stp[0];

			_nominal_states_vector[8] = - nominal_states_msg.pos_stp[2];
			_nominal_states_vector[9] = - nominal_states_msg.vel_stp[2];
			_nominal_states_vector[10] = - nominal_states_msg.acc_stp[2];
			_nominal_states_vector[11] = - nominal_states_msg.jerk_stp[2];

			_nominal_states_vector[12] = -nominal_states_msg.psi + M_PI_2_F;
			_nominal_states_vector[13] = -nominal_states_msg.psi_dot;
		}
	}

	// update measured state vector (data coming from the jerk estimator so it is already in NED frame)
	if (_measured_states_sub.updated())
	{
		states_vector_s measured_states_msg;

		if (_measured_states_sub.copy(&measured_states_msg))
		{
			_measured_states_vector[0] = measured_states_msg.pos_stp[0];
			_measured_states_vector[1] = measured_states_msg.vel_stp[0];
			_measured_states_vector[2] = measured_states_msg.acc_stp[0];
			_measured_states_vector[3] = measured_states_msg.jerk_stp[0];

			_measured_states_vector[4] = measured_states_msg.pos_stp[1];
			_measured_states_vector[5] = measured_states_msg.vel_stp[1];
			_measured_states_vector[6] = measured_states_msg.acc_stp[1];
			_measured_states_vector[7] = measured_states_msg.jerk_stp[1];

			_measured_states_vector[8] = measured_states_msg.pos_stp[2];
			_measured_states_vector[9] = measured_states_msg.vel_stp[2];
			_measured_states_vector[10] = measured_states_msg.acc_stp[2];
			_measured_states_vector[11] = measured_states_msg.jerk_stp[2];

			_measured_states_vector[12] = measured_states_msg.psi;
			_measured_states_vector[13] = measured_states_msg.psi_dot;
		}
	}
}

void Linearizer::HOSM_control(const float nominal_states_vector[14], const float measured_states_vector[14],
				const float v_ref[4], float v_out[4])
{
	float v_SM[4];
	for (size_t i = 0; i < 3; i++)
	{
		float si = measured_states_vector[4*i] - nominal_states_vector[4*i];
		float si_dot = measured_states_vector[4*i+1] - nominal_states_vector[4*i+1];
		float si_2dot = measured_states_vector[4*i+2] - nominal_states_vector[4*i+2];
		float si_3dot = measured_states_vector[4*i+3] - nominal_states_vector[4*i+3];

		v_SM[i] = - _gamma[i] * signf(si_3dot + 3 * powf(powf(si_2dot, 6) + powf(si_dot, 4) + powf(abs(si),3), 1/12.f) *
					signf(si_2dot + powf(powf(si_dot, 4) + powf(abs(si),3), 1/6.f)) *
					signf(si_dot + 1/2.f * powf(abs(si),3/4.f) * signf(si)));

	}

	float s4 = measured_states_vector[12] - nominal_states_vector[12];
	float s4_dot = measured_states_vector[13] - nominal_states_vector[13];
	v_SM[3] = - _gamma[3] * signf(s4_dot + powf(abs(s4), 1/2.f) * signf(s4));

	// control law
	for (size_t i = 0; i < 4; i++)
	{
		v_out[i] = v_ref[i] + v_SM[i];
	}
	//PX4_WARN("HOSM output: [%f, %f, %f, %f]", (double) v_out[0], (double) v_out[1], (double) v_out[2], (double) v_out[3]);

}

void Linearizer::PP_control(const float nominal_states_vector[14], const float measured_states_vector[14], const float k_1[4],
				const float k_2[2], float v_out[4])
{

	matrix::Vector<float,4> k1 = matrix::Vector<float,4>(k_1);
	matrix::Vector<float,2> k2 = matrix::Vector<float,2>(k_2);

	// for the 3 coordinates x,y and z we have 3 fourth-order state spaces
	for (int i = 0; i < 3; i++)
	{
		float diff_data[4] = {measured_states_vector[i*4] - nominal_states_vector[i*4],
				      measured_states_vector[i*4+1] - nominal_states_vector[i*4+1],
				      measured_states_vector[i*4+2] - nominal_states_vector[i*4+2],
				      measured_states_vector[i*4+3] - nominal_states_vector[i*4+3]};
		matrix::Vector<float,4> diff = matrix::Vector<float,4>(diff_data);
		v_out[i] = k1.dot(diff);
	}

	// for the psi we have one second-order state space

	float diff_data[2] = {measured_states_vector[12] - nominal_states_vector[12],
			      measured_states_vector[13] - nominal_states_vector[13]};
	matrix::Vector<float,2> diff = matrix::Vector<float,2>(diff_data);
	v_out[3] = k2.dot(diff);
}

float Linearizer::compute_u1(const float z[14], const float v[4])
{
	// intemediate vectors
	matrix::Vector3f alpha = {z[2], z[6], z[10] + CONSTANT_G};
	matrix::Vector3f beta = {z[3], z[7], z[11]};
	matrix::Vector3f gamma = {v[0], v[1], v[2]};

	// TODO handle alfa.norm() == 0
	float u1_2dot = _mass/alpha.norm() * (alpha.dot(gamma) + beta.norm_squared() -powf(alpha.dot(beta) / alpha.norm(), 2));
	//double integration to get u1
	float u1_dot = discrete_integration(u1_2dot, _u1_2dot_int_state, 1.0f);
	return discrete_integration(u1_dot, _u1_dot_int_state, 1.0f);
}

matrix::Vector3f Linearizer::compute_u2( const float z[14], const float v[4], float u1)
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
	matrix::Vector3f v_sliced = {v[0], v[1], v[2]};
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

matrix::Vector<float, 4> Linearizer::linearize(const float input_vector[4])
{
	// feedforward linearization -> z is the nominal state vectors
	// feedback linearization -> z is the measured state vectors
	float *z = _feed_type== FF ? _nominal_states_vector : _measured_states_vector;

	// compute thrust (u1) and roll, pitch and yaw moments (u2.xyz)
	float u1 = compute_u1(z, input_vector);
	matrix::Vector3f u2 = compute_u2(z, input_vector, u1);

	matrix::Vector<float, 4> U;
	U(0) = u1;
	U.slice<3,1>(1,0) = u2;
	return U;
}

void Linearizer::mix_and_publish(matrix::Vector<float, 4>& U)
{
	// U <-> OMEGA relation: OMEGA^2 = _thrust_to_speed * U
	matrix::Vector<float, 4> W;
	W = _thrust_to_speed * U; // apply equations described in the header to obtain the square of each propeller speed
	W = W.sqrt();

	// remap and publish propeller speeds
	propeller_speeds_s prop_speeds_msg;
	for (size_t i = 0; i < 4; i++)
	{
		if (W(i) < _min_propeller_speed || W(i) > _max_propeller_speed)
		{
			PX4_WARN("Requested speed %f for propeller %zu exceeds range [%g,%g]",
					(double) W(i), i, (double) _min_propeller_speed, (double) _max_propeller_speed );
		}

		//map from rad/s to normalized PWM (in [0;1] range), where rads = PWMnorm*1200 + 100
		W(i) = (W(i) - 100)/1200; //TODO correct with _min_propeller_speed and _max_propeller_speed

		// from [0;1] to [-1;1]
		W(i) = math::constrain((2.f * W(i) - 1.f), -1.f, 1.f);

		if (isnan(W(i)))
		{
			W(i) = -1;
		}

		prop_speeds_msg.propeller_speed[i] = W(i);
	}
	prop_speeds_msg.timestamp = hrt_absolute_time();
	_propeller_speeds_pub.publish(prop_speeds_msg);

	// publish empty actuator controls message to maintain the lockstep synch
	actuator_controls_s actuators{};
	actuators.timestamp = hrt_absolute_time();
	_actuators_pub.publish(actuators);
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

	update_state();

	if (_control_type != 0)
	{
		float v_out[4] = {0.f};

		// compute control action
		if (_control_type == 1) //HOSM
		{
			HOSM_control(_nominal_states_vector, _measured_states_vector, _nominal_inputs, v_out);
		}
		else if (_control_type == 2) //PP
		{
			PP_control(_nominal_states_vector, _measured_states_vector, _k1, _k2, v_out);
		}

		// apply control action to the linearized quadrotor model
		matrix::Vector<float, 4> U = linearize(v_out);

		mix_and_publish(U);
	}

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
