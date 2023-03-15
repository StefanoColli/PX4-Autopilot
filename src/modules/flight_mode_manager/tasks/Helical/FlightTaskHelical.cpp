#include "FlightTaskHelical.hpp"

/*
 * Called when switching to the task, and allows to initialize its state
 * and take over gently from the passed over setpoints (last_setpoint)
 * the previous task was just applying
 */
bool FlightTaskHelical::activate(const vehicle_local_position_setpoint_s &last_setpoint)
{
	bool ret = FlightTask::activate(last_setpoint);

	_starting_pos = _position; //the curent position is the starting point of the helix

	//calculate the finish time to perform the requested number of loops
	_finish_time = (float) _loops * (2 * (float) MATH_PI / _trajectory_speed); //each loop takes 2*pi/speed seconds

	PX4_INFO("Initiating Helical trajectory");
	return ret;
}

/*
 * Called on every loop iteration during the execution and contains the core behaviour
 * implementation producing setpoints.
 * The parametric equation of the helix is
 * x = r*cos(t)
 * y = r*sen(t)
 * z = c*t
 * Where r is the radius of the helix and 2*pi*c gives the vertical separation of the helix loops
 */
bool FlightTaskHelical::update()
{
	//compute time (in seconds) elapsed from the starting of the trajectory
	float time = (_time_stamp_current - _time_stamp_activate) / 1e6f; //timestamps are in ms

	if (time <= _finish_time)
	{
		_position_setpoint(0) = _starting_pos(0) + _helix_radius * cosf(time*_trajectory_speed);
		_position_setpoint(1) = _starting_pos(1) + _helix_radius * sinf(time*_trajectory_speed);
		_position_setpoint(2) = _starting_pos(2) - _vertical_separation * time*_trajectory_speed; //NED frame
	}
	else
	{
		//end holding reached position
		_position_setpoint = _position;
		PX4_INFO("Ended Helical trajectory of %d loops", _loops);
	}

	return true;
}
