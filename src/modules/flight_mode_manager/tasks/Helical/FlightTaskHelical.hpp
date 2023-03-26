#pragma once

#include "FlightTask.hpp"
#include <matrix/matrix/math.hpp>
#include <uORB/topics/helical_trajectory.h>
#include <uORB/Publication.hpp>

class FlightTaskHelical : public FlightTask
{

public:
	FlightTaskHelical() = default;
	virtual ~FlightTaskHelical() = default;

	bool update() override;
	bool activate(const vehicle_local_position_setpoint_s &last_setpoint) override;

private:
	matrix::Vector3f _starting_pos{ 0.f, 0.f, 0.f}; //starting position coord of the helix
	matrix::Vector3f _helix_origin{0.f, 0.f, 0.f}; //origin of the helical trajectory
	static constexpr float _helix_radius = 5.f; //radius of the helix (in meters)
	static constexpr float _vertical_separation = 1.f; //vertical separation between a coil and another (in meters)
	static constexpr float _trajectory_speed = 0.4f; //factor to scale the speed of the trajectory
	static constexpr uint _loops = 3; //number of loops in the helical trajectory
	static constexpr float _tolerance_radius = 0.1f; //tolerance radius for the starting point
	static constexpr float _tolerance_vel = 0.1f; //tolerance velocity for starting the helical trajectory
	hrt_abstime _time_stamp_start_helix = 0; //timestamp at which the actual helical trajectory started (not the flight task)
	float _trajectory_duration; //estimated time duration of the helical trajectory
	uORB::Publication<helical_trajectory_s> _helical_trajectory_publisher{ORB_ID(helical_trajectory)}; //handle to publish trajectory status
};
