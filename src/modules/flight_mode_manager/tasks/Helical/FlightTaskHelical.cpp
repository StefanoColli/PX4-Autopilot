#include "FlightTaskHelical.hpp"

/*
 * Called when switching to the task, and allows to initialize its state
 * and take over gently from the passed over setpoints (last_setpoint)
 * the previous task was just applying
 */
bool FlightTaskHelical::activate(const vehicle_local_position_setpoint_s &last_setpoint)
{
	_helical_trajectory_publisher.advertise();

	bool ret = FlightTask::activate(last_setpoint);

	//the starting point of the helix is at _helix_radius distance (along x) from the current position
	_starting_pos(0) = _position(0) + _helix_radius;
	_starting_pos(1) = _position(1);
	_starting_pos(2) = _position(2);

	//the helix origin is the current position when we call the task
	_helix_origin = _position;

	//calculate the finish time to perform the requested number of loops
	_trajectory_duration = (float) _loops * (2 * (float) MATH_PI / _trajectory_speed); //each loop takes 2*pi/speed seconds

	//PX4_INFO("Initiating Helical trajectory");
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
	helical_trajectory_s msg{};
	msg.timestamp = hrt_absolute_time();

	if (((_position - _starting_pos).longerThan(_tolerance_radius) ||
		_velocity.longerThan(_tolerance_vel)) &&
		_time_stamp_start_helix == 0)
	{
		//we still are too far from the starting position, or too fast
		_position_setpoint = _starting_pos;

		//PX4_INFO("Approaching starting position");
		msg.started = false;
		msg.ended = false;
	}
	else
	{
		//we are in the neighborhood of the starting position

		_velocity_setpoint.setNaN();

		//set the timestamp at which the actual trajectory starts
		if (_time_stamp_start_helix == 0)
		{
			_time_stamp_start_helix = hrt_absolute_time();
		}

		//compute time (in seconds) elapsed from the starting of the trajectory
		float trajectory_time = (_time_stamp_current - _time_stamp_start_helix) / 1e6f; //timestamps are in ms

		//helical trajectory
		if (trajectory_time <= _trajectory_duration)
		{
			_position_setpoint(0) = _helix_origin(0) + _helix_radius * cosf(trajectory_time*_trajectory_speed);
			_position_setpoint(1) = _helix_origin(1) + _helix_radius * sinf(trajectory_time*_trajectory_speed);
			_position_setpoint(2) = _helix_origin(2) - _vertical_separation * trajectory_time*_trajectory_speed; //NED frame

			//PX4_INFO("Following helical trajectory");
			msg.started = true;
			msg.ended = false;
		}
		else
		{
			//end holding reached position
			//_position_setpoint = _position;


			//PX4_INFO("Helical trajectory terminated");
			msg.started = true;
			msg.ended = true;
		}
	}

	_helical_trajectory_publisher.publish(msg);
	return true;
}
