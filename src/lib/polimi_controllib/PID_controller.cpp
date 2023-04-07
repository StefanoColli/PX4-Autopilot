#include "PID_controller.h"

PID_controller::PID_controller(float Kc, float Ti, float Td, float N, float b, float c, float uMin, float uMax) {

    // Check parameter consistency
    if ((Kc <= FLT_EPSILON) && (Ti <= FLT_EPSILON) && (Td <= FLT_EPSILON))
        PX4_ERR("PID parameters cannot be all zero");
    if (uMin > uMax)
        PX4_ERR("Lower control saturation cannot be greater then higher");

    // Initialise PID parameters
    this->_Kc = Kc;
    this->_Ti = Ti;
    this->_Td = Td;
    this->_N = N;
    this->_b = b;
    this->_c = c;
    this->_uMin = uMin;
    this->_uMax = uMax;


    // Initialise PID state variables
    _y_old = _ysp_old = _u_old = _uI = _uD = 0.0;

    // Initialize PID logic states
    _controller_state = PID_state::AUTO;
    _actuation_state = actuator_state::NO_SATURATION;
    _controlSignal_state = control_state::NO_FREEZE;
}

PID_controller::~PID_controller() {
    // Do nothing
}

void PID_controller::evaluate(float y, float ysp, float delta_T, float trk, float &u) {

    if (delta_T <= FLT_EPSILON)
        PX4_ERR("PID sampling time cannot be zero");

    // Compute PID equation internal coefficients a1,b1 and b2
    float a1, b1, b2;

    if (_Ti <= FLT_EPSILON)                        // No integral action
        a1 = 0.0;
    else {
        if (_Kc <= FLT_EPSILON)                    // Integral action without proportional action
            a1 = delta_T / _Ti;
        else                              // Integral action with proportional action
            a1 = _Kc * delta_T / _Ti;
    }

    if (_Td <= FLT_EPSILON)                        // No derivative action
        b1 = b2 = 0.0;
    else {
        if (_Kc <= FLT_EPSILON)                    // Derivative action without proportional action
        {
            b1 = _Td / (_N * delta_T + _Td);
            b2 = _N * b1;
        } else                            // Derivative action with proportional action
        {
            b1 = _Td / (_N * delta_T + _Td);
            b2 = _Kc * _N * b1;
        }
    }

    // Compute error
    float e = ysp - y;

    // Compute derivative action
    _uD = b1 * _uD + b2 * (_c *(ysp - _ysp_old) - (y - _y_old));

    // Manage AUTO/TRACKING behaviour
    if (_controller_state == PID_state::AUTO) {
        u = _Kc * (_b * ysp - y) + _uI + _uD;
    } else if (_controller_state == PID_state::TRACKING) {
        u = trk;
    } else {
        u = 0.0;
        PX4_ERR("Unknown PID state");
    }

    // Check windup and freeze conditions and update integral contribution
    if (u > _uMax) {
        u = _uMax;
        _actuation_state = actuator_state::SATURATION_UP;
    } else if (u < _uMin) {
        u = _uMin;
        _actuation_state = actuator_state::SATURATION_DOWN;
    } else if ((_controller_state == PID_state::AUTO) && (_controlSignal_state == control_state::NO_FREEZE)) {
        _uI += a1 * e;
        _actuation_state = actuator_state::NO_SATURATION;
    } else if ((_controller_state == PID_state::AUTO) && (_controlSignal_state == control_state::FREEZE_UP)) {
        if (u > _u_old) {
            u = _u_old;
            _uI = u - _Kc * (_b * ysp - y) - _uD;
        } else {
            _uI += a1 * e;
        }
        _actuation_state = actuator_state::NO_SATURATION;
    } else if ((_controller_state == PID_state::AUTO) && (_controlSignal_state == control_state::FREEZE_DOWN)) {
        if (u < _u_old) {
            u = _u_old;
            _uI = u - _Kc * (_b * ysp - y) - _uD;
        } else {
            _uI += a1 * e;
        }
        _actuation_state = actuator_state::NO_SATURATION;
    } else {
        //in tracking mode the integral action tracks the integral part of the trk signal given
        // so to have a bumpless transition between tracking and automatic mode
        _uI = u - _Kc * (_b * ysp - y) - _uD;
        _actuation_state = actuator_state::NO_SATURATION;
    }

    // Update PID states
    _y_old = y;
    _ysp_old = ysp;
    _u_old = u;
}

void PID_controller::resetState() {
    // Reset PID state
    _y_old = _ysp_old = _u_old = _uI = _uD = 0.0;
}
