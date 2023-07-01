/****************************************************************************
 *
 *   Copyright (c) 2019 PX4 Development Team. All rights reserved.
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

/**
 * @file RateControl.cpp
 */

#include <RateControl.hpp>
#include <px4_platform_common/defines.h>

using namespace matrix;

void RateControl::setGains(const Vector3f &P, const Vector3f &I, const Vector3f &D)
{
	_rate_x_controller.updatePIDParameters(P(0), P(0)/I(0), D(0)/P(0), NAN, NAN, NAN, -_lim_int(0), _lim_int(0));
	_rate_y_controller.updatePIDParameters(P(1), P(1)/I(1), D(1)/P(1), NAN, NAN, NAN, -_lim_int(1), _lim_int(1));
	_rate_z_controller.updatePIDParameters(P(2), P(2)/I(2), D(2)/P(2), NAN, NAN, NAN, -_lim_int(2), _lim_int(2));
}

void RateControl::setSaturationStatus(const Vector<bool, 3> &saturation_positive,
				      const Vector<bool, 3> &saturation_negative)
{
	_control_allocator_saturation_positive = saturation_positive;
	_control_allocator_saturation_negative = saturation_negative;
}

Vector3f RateControl::update(const Vector3f &rate, const Vector3f &rate_sp, const Vector3f &angular_accel,
			     const float dt, const bool landed)
{
	// PID angular rate control
	if (!landed)
	{
		float u_x, u_y, u_z;
		if (PX4_ISFINITE(rate_sp(0)) && PX4_ISFINITE(rate(0)))
		{
			_rate_x_controller.evaluate(rate(0), rate_sp(0), dt, 0.0, u_x);
		}
		else
		{
			u_x = NAN;
		}

		if (PX4_ISFINITE(rate_sp(1)) && PX4_ISFINITE(rate(1)))
		{
			_rate_y_controller.evaluate(rate(1), rate_sp(1), dt, 0.0, u_y);
		}
		else
		{
			u_y = NAN;
		}

		if (PX4_ISFINITE(rate_sp(2)) && PX4_ISFINITE(rate(2)))
		{
			_rate_z_controller.evaluate(rate(2), rate_sp(2), dt, 0.0, u_z);
		}
		else
		{
			u_z = NAN;
		}

		return Vector3f(u_x, u_y, u_z) + _gain_ff.emult(rate_sp);
	}
	else
	{
		return Vector3f(0.0f, 0.0f, 0.0f);
	}
}

void RateControl::getRateControlStatus(rate_ctrl_status_s &rate_ctrl_status)
{
	rate_ctrl_status.rollspeed_integ = _rate_x_controller.getIntegral();
	rate_ctrl_status.pitchspeed_integ = _rate_y_controller.getIntegral();
	rate_ctrl_status.yawspeed_integ = _rate_z_controller.getIntegral();
}
