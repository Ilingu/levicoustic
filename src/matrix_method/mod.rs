pub mod potential_field;
pub mod pressure_field;
pub mod radiation_force;
use ndarray::Array2;

use num_complex::Complex;
use std::{f64::consts::PI, fmt::Display};

pub type Field = Array2<Complex<f64>>;

#[allow(dead_code)]
#[derive(Clone, Copy, PartialEq)]
#[repr(u8)]
pub enum FieldType {
    Pressure,
    Velocity,
    RadiationPotential,
    RadiationForce,
    Intensity,
    LimitDensity,
}

impl FieldType {
    pub fn to_unit(self) -> String {
        match self {
            FieldType::Pressure => "Acoustic pressure (Pa)",
            FieldType::Velocity => "Velocity (m/s)",
            FieldType::RadiationPotential => "Relative acoutic potential (N/m²)",
            FieldType::RadiationForce => "Radiation force (N)",
            FieldType::Intensity => "Intensité sonore (W/m²)",
            FieldType::LimitDensity => "Masse volumique limite (kg/m³)",
        }
        .to_string()
    }
}

impl Display for FieldType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FieldType::Pressure => write!(f, "Pressure field"),
            FieldType::Velocity => write!(f, "Velocity field"),
            FieldType::RadiationPotential => write!(f, "Relative acoustic radiation potential"),
            FieldType::RadiationForce => write!(f, "Radiation force"),
            FieldType::Intensity => write!(f, "Intensité"),
            FieldType::LimitDensity => write!(f, "Masse volumique limite"),
        }
    }
}

/* CONSTANT DEFINITION */

pub const MM: f64 = 1e-3;
const I: Complex<f64> = Complex::I;

// Define air constants at 20°C (https://en.wikipedia.org/wiki/Density_of_air)
// pub const RHO: f64 = 1.2041; // Density of air in Raleigh, kg/m^3
// pub const C: f64 = 343.21; // Speed of sound in air, m/s

pub const RHO: f64 = 1.214; // Density of air in Raleigh, kg/m^3
pub const C: f64 = 340.1; // Speed of sound in air, m/s

/* Parameters */

#[derive(Debug, Clone, Copy)]
pub struct SimulationParametersArgs {
    // Grid dimensions
    pub x_min: f64,
    pub x_max: f64,
    /// Location of the reflector
    pub z_min: f64,
    /// Location of the transducer
    pub z_max: f64,

    // Simulation accuracy
    pub nb_reflection: u8,
    /// Discretization step size, m
    pub disc: f64,

    // Transducer parameters
    /// Transducer offset on x-axis, m
    pub offset: f64,
    /// Transducer radius, m
    pub radius: f64,
    /// Transducer hole, m
    pub hole_radius: f64,
    /// Resonante Frequency, Hz
    pub freq: f64,
    /// Transducer phase offset, degree
    pub phase: f64,
    /// Displacement amplitude, m          
    pub u_0: f64,
    /// Transducer tilt, deg
    pub inclination: f64,

    // Ball parameter
    /// Radius of the a small sphere in the field, m.
    ///
    /// **Should be much smaller than the wavelength**
    pub sphere_radius: f64,
}

impl From<SimulationParametersArgs> for SimulationParameters {
    fn from(
        SimulationParametersArgs {
            x_min,
            x_max,
            z_min,
            z_max,

            nb_reflection,
            disc,

            offset,
            radius,
            hole_radius,
            freq,
            phase,
            u_0,
            inclination,
            ..
        }: SimulationParametersArgs,
    ) -> Self {
        let mut sp = Self {
            x_min,
            x_max,
            z_min,
            z_max,

            nb_reflection,
            disc,

            offset,
            radius,
            transducer_area: PI * radius * radius,
            hole_radius,
            hole_area: PI * hole_radius * hole_radius,
            freq,
            phase,
            u_0,
            inclination,

            reflector_area: PI * ((x_max * x_max) / 4.0),

            wavelength: C / freq,
            omega: 2.0 * PI * freq,
            // to be computed
            wavenumber: 0.0,
            e: Complex::ZERO,
            d: 0.0,
        };
        sp.wavenumber = sp.omega / C;
        sp.e = Complex::new(0.0, 1.0 / sp.wavelength);
        sp.d = (sp.omega * RHO * C) / sp.wavelength;
        sp
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SimulationParameters {
    // Grid dimensions
    x_min: f64,
    x_max: f64,
    z_min: f64, // Location of the reflector
    z_max: f64, // Location of the transducer

    // Simulation accuracy
    nb_reflection: u8,
    disc: f64, // Discretization step size, m

    // Transducer parameters
    offset: f64,          // Transducer offset on x-axis, m
    radius: f64,          // Transducer radius, m
    transducer_area: f64, // Transducer area, m^2
    hole_radius: f64,     // Transducer hole, m
    hole_area: f64,       // Transducer hole area, m^2
    freq: f64,            // Resonante Frequency, Hz
    phase: f64,           // Transducer phase offset, degree
    u_0: f64,             // Displacement amplitude, m
    inclination: f64,     // Transducer tilt, degree

    // Reflector parameters
    reflector_area: f64,

    // Other computed constants
    wavelength: f64, // Wavelength, m
    omega: f64,      // Angular frequency, rad/s
    wavenumber: f64, // Wavenumber - k
    e: Complex<f64>, // Constant used for reflected waves
    d: f64,          // Constant used for transmitted wave
}

/* DEFAULT SIMULATION */

impl Default for SimulationParametersArgs {
    fn default() -> Self {
        Self {
            x_min: -50.0 * MM,
            x_max: 50.0 * MM,
            z_min: 0.0 * MM,
            z_max: 50.0 * MM,

            nb_reflection: 4,
            disc: 0.4 * MM,

            offset: 0.0 * MM,
            radius: 15.0 * MM,
            hole_radius: 2.0 * MM,
            freq: 56000.0,
            phase: 0.0,
            u_0: 0.0000060,
            inclination: 0.0,

            sphere_radius: 0.1 * MM, // wavelength=6.1mm
        }
    }
}

impl Default for SimulationParameters {
    fn default() -> Self {
        SimulationParametersArgs::default().into()
    }
}

/* HELPERS */
fn multiple_matrix_product(matrices: Vec<&Array2<Complex<f64>>>) -> Array2<Complex<f64>> {
    let mut prod = matrices[0].to_owned();
    for matrix in matrices.into_iter().skip(1) {
        prod = prod.dot(matrix);
    }
    prod
}

fn add_matrices(matrices: Vec<&Array2<Complex<f64>>>) -> Array2<Complex<f64>> {
    let mut sum = matrices[0].to_owned();
    for matrix in matrices.into_iter().skip(1) {
        sum += matrix;
    }
    sum
}
fn add_owned_matrices(matrices: Vec<Array2<Complex<f64>>>) -> Array2<Complex<f64>> {
    let mut sum = matrices[0].to_owned();
    for matrix in matrices.into_iter().skip(1) {
        sum += &matrix;
    }
    sum
}

#[allow(dead_code)]
pub fn u0_from_velocity_amp(velocity_amp: f64, freq: f64) -> f64 {
    velocity_amp / freq
}
