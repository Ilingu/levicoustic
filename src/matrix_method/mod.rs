pub mod pressure_field;
use static_assertions::const_assert;

use num_complex::Complex;
use std::f64::consts::PI;

/* CONSTANT DEFINITION */

const MM: f64 = 1e-3;

// Grid dimensions
pub const X_MIN: f64 = -50.0 * MM;
pub const X_MAX: f64 = 50.0 * MM;
pub const Z_MIN: f64 = 0.0; // Location of the reflector
pub const Z_MAX: f64 = 50.0 * MM; // Location of the transducer
const DISC: f64 = 0.50 * MM; // Discretization step size

const_assert!(X_MAX > X_MIN);
const_assert!(Z_MAX > Z_MIN);

// Define air constants
const RHO: Complex<f64> = Complex::new(1.214, 0.0); // Density of air in Raleigh, kg/m^3
const C: Complex<f64> = Complex::new(340.1, 0.0); // Speed of sound in air, m/s

// Define transducer parameters
const R: f64 = 15.0 * MM; // Transducer radius, m
const TRANSDUCER_AREA: Complex<f64> = Complex::new(PI * R * R, 0.0); // Transducer area, m^2
const R2: f64 = 2.0 * MM; // Transducer hole, m
const HOLE_AREA: Complex<f64> = Complex::new(PI * R2 * R2, 0.0); // Transducer hole area, m^2
const FREQ: Complex<f64> = Complex::new(56000.0, 0.0); // Frequency, Hz
const U_0: f64 = 0.0000060; // Displacement amplitude, m

const_assert!(R > R2);

// Define reflector parameters
const REFLECTOR_AREA: Complex<f64> = Complex::new(PI * ((X_MAX * X_MAX) / 4.0), 0.0); // Reflector area, m^2

// Constants
const WAVELENGTH: Complex<f64> = Complex::new(C.re / FREQ.re, 0.0); // Wavelength, m
const OMEGA: Complex<f64> = Complex::new(2.0 * PI * FREQ.re, 0.0); // Angular frequency, rad/s
const WAVENUMBER: Complex<f64> = Complex::new(OMEGA.re / C.re, 0.0); // Wavenumber - k
const E: Complex<f64> = Complex::new(0.0, 1.0 / WAVELENGTH.re); // Constant used for reflected waves
const D: Complex<f64> = Complex::new((OMEGA.re * RHO.re * C.re) / WAVELENGTH.re, 0.0); // Constant used for transmitted wave
