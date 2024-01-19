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

#![no_std]

use core::f32::consts::{FRAC_1_SQRT_2, PI};

use micromath::F32Ext;

const PI2: f32 = PI * 2.0;

struct ThirdOrderIntegrator {
    /// Gain is sample_time / 12
    integrator_gain: f32,
    z1: f32,
    z2: f32,
    z3: f32,
}

impl ThirdOrderIntegrator {
    /// Create a new integrator with a given sample time
    fn new(integrator_gain: f32) -> Self {
        Self {
            integrator_gain,
            z1: 0.0,
            z2: 0.0,
            z3: 0.0,
        }
    }

    /// Update the integrator with a new input
    fn update(&mut self, x: f32) {
        self.z3 = self.z2;
        self.z2 = self.z1;
        self.z1 += x * self.integrator_gain;
    }

    /// Get the current value of the integrator
    fn value(&self) -> f32 {
        self.z1 * 23.0 - self.z2 * 16.0 + self.z3 * 5.0
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
        let integrator_gain = sample_time / 12.0;
        Self {
            k,
            integrator_1: ThirdOrderIntegrator::new(integrator_gain),
            integrator_2: ThirdOrderIntegrator::new(integrator_gain),
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
    let (sin, cos) = theta.sin_cos();
    // let d = alpha * cos + beta * sin;
    -alpha * sin + beta * cos
}

/// Configuration for the SOGI-PLL
pub struct PllConfig {
    pub sample_time: f32,
    pub sogi_k: f32,
    pub pi_proportional_gain: f32,
    pub pi_integral_gain: f32,
    pub omega_zero: f32,
}

/// Result returned by the SOGI-PLL
pub struct PllResult {
    pub v_alpha: f32,
    pub v_beta: f32,
    pub omega: f32,
    pub theta: f32,
}

impl PllResult {
    pub fn v_rms(&self) -> f32 {
        FRAC_1_SQRT_2 * (self.v_alpha * self.v_alpha + self.v_beta * self.v_beta).sqrt()
    }
}

/// SOGI-PLL implementation
pub struct SogiPll {
    config: PllConfig,
    sogi: Sogi<f32>,
    pi_integral: f32,
    pi_value: f32,
    z1: f32,
}

impl SogiPll {
    /// Create a new SOGI-PLL with a given configuration
    pub fn new(config: PllConfig) -> SogiPll {
        let sogi = Sogi::new(config.sogi_k, config.sample_time);

        SogiPll {
            config,
            sogi,
            pi_integral: 0.0,
            pi_value: 0.0,
            z1: 0.0,
        }
    }

    /// Update the PLL with a new voltage measurement
    pub fn update(&mut self, v: f32) -> PllResult {
        let omega = self.pi_value + self.config.omega_zero;
        let (v_alpha, v_beta) = self.sogi.update(v, omega);

        let q = alpha_beta_to_q(v_alpha, v_beta, self.z1);

        self.z1 = (omega * self.config.sample_time + self.z1) % PI2;

        self.pi_integral += self.pi_value * self.config.sample_time;
        self.pi_value =
            q * self.config.pi_proportional_gain + self.pi_integral * self.config.pi_integral_gain;

        PllResult {
            v_alpha,
            v_beta,
            omega,
            theta: self.z1,
        }
    }
}
