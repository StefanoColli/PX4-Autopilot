#pragma once

#include "FlightTask.hpp"
#include <matrix/matrix/math.hpp>

class FlightTaskHelical : public FlightTask
{

public:
	FlightTaskHelical() = default;
	virtual ~FlightTaskHelical() = default;

	bool update() override;
	bool activate(const vehicle_local_position_setpoint_s &last_setpoint) override;

private:
	matrix::Vector3f _starting_pos{ 0.f, 0.f, 0.f}; //starting z coord of the helix
	static constexpr float _helix_radius = 5.f; //radius of the helix (in meters)
	static constexpr float _vertical_separation = 1.f; //vertical separation between a coil and another (in meters)
	static constexpr float _trajectory_speed = 0.5f; //factor to scale the speed of the trajectory
	static constexpr uint _loops = 3; //number of loops in the helical trajectory
	float _finish_time;
};
