//! SOGI-PLL implementation
//!
//! Based on: <https://www.researchgate.net/publication/345738110_Instantaneous_Symmetrical_Component_Estimator_Using_Second_Order_Generalized_Integrator_under_Distorted_Voltage_Conditions>
//!
//!
//! Usage:
//!
//! ```rust
//! let mut pll = SogiPll::new(config);
//!
//! // Call .update every sample_time s
//! let result = pll.update(measured_voltage);
//! ```
//!
//! Recommended parameters:
//!
//! ```rust
//! sogi_k: 1.0
//! pi_proportional_gain: 178.0,
//! pi_integral_gain: 0.0001
//! sample_rate > 1000Hz
//! ```
//! K, Kp, Ki for 50Hz taken from [here](https://ieeexplore.ieee.org/document/6636494)

use core::f32::consts::{PI, SQRT_2};

use pid::Pid;

const SQRT_1_2: f32 = 1.0 / SQRT_2;

struct ThirdOrderIntegrator {
    sample_time: f32,
    samples: [f32; 3],
}

impl ThirdOrderIntegrator {
    /// Create a new integrator with a given sample time
    fn new(sample_time: f32) -> Self {
        Self {
            sample_time,
            samples: [0.0; 3],
        }
    }

    /// Update the integrator with a new input
    fn update(&mut self, x: f32) {
        self.samples[2] = self.samples[1];
        self.samples[1] = self.samples[0];
        self.samples[0] += x * (self.sample_time / 12.0);
    }

    /// Get the current value of the integrator
    fn value(&self) -> f32 {
        self.samples[0] * 23.0 - self.samples[1] * 16.0 + self.samples[2] * 5.0
    }
}

/// SOGI implementation that uses third order integrators
pub struct Sogi<T> {
    k: T,
    integrator_1: ThirdOrderIntegrator,
    integrator_2: ThirdOrderIntegrator,
}

impl Sogi<f32> {
    /// Create a new SOGI with a given k and sample time(used for the integrators)
    pub fn new(k: f32, sample_time: f32) -> Self {
        Self {
            k,
            integrator_1: ThirdOrderIntegrator::new(sample_time),
            integrator_2: ThirdOrderIntegrator::new(sample_time),
        }
    }

    /// Update SOGI with a new voltage measurement
    /// Returns (v_alpha, v_beta)
    pub fn update(&mut self, v: f32, omega: f32) -> (f32, f32) {
        let v_alpha = self.integrator_1.value();
        let v_beta = self.integrator_2.value();

        let integrator_1_in = ((v - v_alpha) * self.k - v_beta) * omega;
        let integrator_2_in = v_alpha * omega;

        self.integrator_1.update(integrator_1_in);
        self.integrator_2.update(integrator_2_in);

        (v_alpha, v_beta)
    }
}

fn alpha_beta_to_q(alpha: f32, beta: f32, theta: f32) -> f32 {
    // let d = alpha * theta.cos() + beta * theta.sin();
    -alpha * theta.sin() + beta * theta.cos()
}

/// Configuration for the SOGI-PLL
pub struct PllConfig {
    pub sample_time: f32,
    pub sogi_k: f32,
    pub pi_proportional_gain: f32,
    pub pi_integral_gain: f32,
    pub omega_zero: f32,
    pub frequency_limit: f32,
}

/// Result returned by the SOGI-PLL
pub struct PllResult {
    pub v_beta: f32,
    pub theta: f32,
    pub v_rms: f32,
}

/// SOGI-PLL implementation
pub struct SogiPll {
    config: PllConfig,
    sogi: Sogi<f32>,
    pid: Pid<f32>,
    pi_value: f32,
    z1: f32,
}

impl SogiPll {
    /// Create a new SOGI-PLL with a given configuration
    pub fn new(config: PllConfig) -> SogiPll {
        let sogi = Sogi::new(config.sogi_k, config.sample_time);

        let mut pid = Pid::new(0.0, config.frequency_limit);
        pid.p(config.pi_proportional_gain, config.frequency_limit);
        pid.i(config.pi_integral_gain, config.frequency_limit);

        SogiPll {
            config,
            sogi,
            pid,
            pi_value: 0.0,
            z1: 0.0,
        }
    }

    /// Update the PLL with a new voltage measurement
    pub fn update(&mut self, v: f32) -> PllResult {
        let omega = self.pi_value + self.config.omega_zero;
        let (v_alpha, v_beta) = self.sogi.update(v, omega);

        let q = alpha_beta_to_q(v_alpha, v_beta, self.z1);

        self.z1 = (omega * self.config.sample_time + self.z1) % (2.0 * PI);

        // The controller internally does setpoint - measured to get the error
        // We are supplying the error directly so we need to negate it
        let pid_out = self.pid.next_control_output(-q);
        self.pi_value = pid_out.output;

        let v_rms = SQRT_1_2 * (v_alpha.powi(2) + v_beta.powi(2)).sqrt();

        PllResult {
            v_beta,
            theta: self.z1,
            v_rms,
        }
    }
}
