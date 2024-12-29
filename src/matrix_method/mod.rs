pub mod potential_field;
pub mod pressure_field;
use ndarray::Array2;

use num_complex::Complex;
use std::{f64::consts::PI, fmt::Display};

pub type Field = Array2<Complex<f64>>;

#[derive(Clone, Copy)]
#[repr(u8)]
pub enum FieldType {
    Pressure,
    Velocity,
    RadiationPotential,
}

impl FieldType {
    pub fn to_unit(self) -> String {
        match self {
            FieldType::Pressure => "Acoustic pressure (Pa)",
            FieldType::Velocity => "Velocity (m/s)",
            FieldType::RadiationPotential => "Relative acoutic potential (N/m²)",
        }
        .to_string()
    }
}

impl Display for FieldType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FieldType::Pressure => write!(f, "Pressure field"),
            FieldType::Velocity => write!(f, "Velocity field"),
            FieldType::RadiationPotential => write!(f, "Acoustic radiation potential"),
        }
    }
}

/* CONSTANT DEFINITION */

pub const MM: f64 = 1e-3;
const I: Complex<f64> = Complex::I;

// Define air constants at 20°C (https://en.wikipedia.org/wiki/Density_of_air)
// const RHO: Complex<f64> = Complex::new(1.2041, 0.0); // Density of air in Raleigh, kg/m^3
// const C: Complex<f64> = Complex::new(343.21, 0.0); // Speed of sound in air, m/s

const RHO: Complex<f64> = Complex::new(1.214, 0.0); // Density of air in Raleigh, kg/m^3
const C: Complex<f64> = Complex::new(340.1, 0.0); // Speed of sound in air, m/s

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
    pub nb_of_reflection: u8,
    /// Discretization step size, m
    pub disc: f64,

    // Transducer parameters
    /// Transducer radius, m
    pub radius: f64,
    /// Transducer hole, m
    pub hole_radius: f64,
    /// Resonante Frequency, Hz
    pub freq: Complex<f64>,
    /// Displacement amplitude, m          
    pub u_0: f64,
}

impl From<SimulationParametersArgs> for SimulationParameters {
    fn from(
        SimulationParametersArgs {
            x_min,
            x_max,
            z_min,
            z_max,
            nb_of_reflection,
            disc,
            radius,
            hole_radius,
            freq,
            u_0,
            ..
        }: SimulationParametersArgs,
    ) -> Self {
        let mut sp = Self {
            x_min,
            x_max,
            z_min,
            z_max,
            nb_of_reflection,
            disc,
            radius,
            transducer_area: Complex::new(PI * radius * radius, 0.0),
            hole_radius,
            hole_area: Complex::new(PI * hole_radius * hole_radius, 0.0),
            freq,
            u_0,
            reflector_area: Complex::new(PI * ((x_max * x_max) / 4.0), 0.0),
            wavelength: Complex::new(C.re / freq.re, 0.0),
            omega: Complex::new(2.0 * PI * freq.re, 0.0),
            // to be computed
            wavenumber: Complex::ZERO,
            e: Complex::ZERO,
            d: Complex::ZERO,
        };
        sp.wavenumber = Complex::new(sp.omega.re / C.re, 0.0);
        sp.e = Complex::new(0.0, 1.0 / sp.wavelength.re);
        sp.d = Complex::new((sp.omega.re * RHO.re * C.re) / sp.wavelength.re, 0.0);
        sp
    }
}

#[derive(Debug, Clone)]
pub struct SimulationParameters {
    // Grid dimensions
    x_min: f64,
    x_max: f64,
    z_min: f64, // Location of the reflector
    z_max: f64, // Location of the transducer

    // Simulation accuracy
    nb_of_reflection: u8,
    disc: f64, // Discretization step size, m

    // Transducer parameters
    radius: f64,                   // Transducer radius, m
    transducer_area: Complex<f64>, // Transducer area, m^2
    hole_radius: f64,              // Transducer hole, m
    hole_area: Complex<f64>,       // Transducer hole area, m^2
    freq: Complex<f64>,            // Resonante Frequency, Hz
    u_0: f64,                      // Displacement amplitude, m

    // Reflector parameters
    reflector_area: Complex<f64>,

    // Other computed constants
    wavelength: Complex<f64>, // Wavelength, m
    omega: Complex<f64>,      // Angular frequency, rad/s
    wavenumber: Complex<f64>, // Wavenumber - k
    e: Complex<f64>,          // Constant used for reflected waves
    d: Complex<f64>,          // Constant used for transmitted wave
}

/* DEFAULT SIMULATION */

impl Default for SimulationParametersArgs {
    fn default() -> Self {
        Self {
            x_min: -50.0 * MM,
            x_max: 50.0 * MM,
            z_min: 0.0 * MM,
            z_max: 50.0 * MM,
            nb_of_reflection: 4,
            disc: 0.25 * MM,
            radius: 15.0 * MM,
            hole_radius: 2.0 * MM,
            freq: Complex::new(56000.0, 0.0),
            u_0: 0.0000060,
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
