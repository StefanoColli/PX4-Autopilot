#ifndef PID_CONTROLLER_H_
#define PID_CONTROLLER_H_

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PID controller class
//
// Implements the transfer function of a digital PID controller, that is obtained applying Backward-Euler transformation
// to the following continuous time transfer function
//
//           [       1          s Td     ]
// R(s) = Kc [ 1 + ------ + ------------ ]
//           [      s Ti     1 + s Td/N  ]
//
// An anti-windup, based on conditional integration, and a tracking action are also implemented.
// When PID is in TRACKING mode the output signal is equal to the tracking signal.
// No bumpless transition from AUTO to TRACKING is guaranteed as during tracking mode the PID is disconnected from the
// plant; bumpless transition from TRACKING to AUTO is guaranteed if the tracking signal is modulated accordingly.
//
// For cascaded control:
// - when the inner controller saturates up/down set the outer controller to freeze up/down
// - when the inner controller is set to TRACKING the outer controller should be set to TRACKING as well
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "float.h"

class PID_controller {
public:
    enum class PID_state {
        AUTO, TRACKING
    };
    enum class actuator_state {
        NO_SATURATION, SATURATION_UP, SATURATION_DOWN
    };
    enum class control_state {
        NO_FREEZE, FREEZE_UP, FREEZE_DOWN
    };

    /////// Named Constructors ///////

    //Proportional controller standard ISA form (P -> Td = N = Ti = 0)
    static inline PID_controller P_std(float Kc, float b, float c, float uMin, float uMax)
                    {return PID_controller(Kc, 0.0f, 0.0f, 0.0f, b, c, uMin, uMax); }
    //Proportional integral controller in standard ISA form (PI -> Td = N = 0)
    static inline PID_controller PI_std(float Kc, float Ti, float b, float c, float uMin, float uMax)
                    {return PID_controller(Kc, Ti, 0.0f, 0.0f, b, c, uMin, uMax); }
    //Proportional derivative controller in standard ISA form (PD -> Ti = 0)
    static inline PID_controller PD_std(float Kc, float Td, float N, float b, float c, float uMin, float uMax)
                    {return PID_controller(Kc, 0.0f, Td, N, b, c, uMin, uMax); }
    //PID controller in standard ISA form
    static inline PID_controller PID_std(float Kc, float Ti, float Td, float N, float b, float c, float uMin, float uMax)
                    {return PID_controller(Kc, Ti, Td, N, b, c, uMin, uMax); }

    //Proportional controller in academic form (P -> Td = N = Ti = 0, b = c = 1)
    static inline PID_controller P_academic(float Kc, float uMin, float uMax)
                    {return PID_controller(Kc, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, uMin, uMax); }
    //Proportional integral controller in academic form (PI -> Td = N = 0, b = c = 1)
    static inline PID_controller PI_academic(float Kc, float Ti, float uMin, float uMax)
                    {return PID_controller(Kc, Ti, 0.0f, 0.0f, 1.0f, 1.0f, uMin, uMax); }
    //Proportional derivative controller in academic form (PD -> Ti = 0, b = c = 1)
    static inline PID_controller PD_academic(float Kc, float Td, float N, float uMin, float uMax)
                    {return PID_controller(Kc, 0.0f, Td, N, 1.0f, 1.0f, uMin, uMax); }
    //PID controller in academic form (b = c = 1)
    static inline PID_controller PID_academic(float Kc, float Ti, float Td, float N, float uMin, float uMax)
                    {return PID_controller(Kc, Ti, Td, N, 1.0f, 1.0f, uMin, uMax); }

    //Proportional controller in parallel form (P -> Ki = Kd = 0, b = c = 1)
    static inline PID_controller P_parallel(float Kp, float uMin, float uMax)
                    {return PID_controller(Kp, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, uMin, uMax); }
    //Proportional integral controller in parallel form (PI -> Kd = 0, b = c = 1)
    static inline PID_controller PI_parallel(float Kp, float Ki, float uMin, float uMax)
                    {return PID_controller(Kp, Kp/Ki, 0.0f, 0.0f, 1.0f, 1.0f, uMin, uMax); }
    //Proportional derivative controller in parallel form (PD -> Ki = 0, b = c = 1)
    static inline PID_controller PD_parallel(float Kp, float Kd, float N, float uMin, float uMax)
                    {return PID_controller(Kp, 0.0f, Kd/Kp, N, 1.0f, 1.0f, uMin, uMax); }
    //PID controller in parallel form (b = c = 1)
    static inline PID_controller PID_parallel(float Kp, float Ki, float Kd, float N, float uMin, float uMax)
                    {return PID_controller(Kp, Kp/Ki, Kd/Kp, N, 1.0f, 1.0f, uMin, uMax); }

    ////////////////////////////

    ~PID_controller();

    /**
     * @brief Compute control action
     *
     * @param y current measure of the control variable
     * @param ysp target value of the control variable
     * @param delta_T time elapsed from last evaluation
     * @param trk reference to be followed in tracking mode
     * @param u output control action
     */
    void evaluate(float y, float ysp, float delta_T, float trk, float &u);

    void resetState();


    PID_state getControllerState() { return _controller_state; }

    void setControllerState(PID_state controller_state);

    /**
     * @brief Update Controller parameters, use NAN to leave the corresponding parameter as it is
     *
     * @param Kc
     * @param Ti
     * @param Td
     * @param N
     * @param b
     * @param c
     * @param uMin
     * @param uMax
     */
    void updatePIDParameters(float Kc, float Ti, float Td, float N, float b, float c, float uMin, float uMax);

    actuator_state getActuationState() { return _actuation_state; }

    control_state getControlSignalState() { return _controlSignal_state; }

    void setControlSignalState(control_state controlSignal_state) { this->_controlSignal_state = controlSignal_state; }

protected:
    /* Class constructor (used by the named constructors)*/
    PID_controller(float Kc, float Ti, float Td, float N, float b, float c, float uMin, float uMax);

private:
    /* PID parameters */
    float _Kc, _Ti, _Td, _N;            // PID parameters (proportional gain, integral/derivative time)
    float _uMin, _uMax;                 // Output saturation values
    float _b, _c;                       //

    /* PID state variables */
    float _u_old, _y_old, _ysp_old, _uI, _uD;     // PID internal states

    /* PID logic states */
    PID_state _controller_state;
    actuator_state _actuation_state;
    control_state _controlSignal_state;
};

#endif /* PID_CONTROLLER_H_ */
