use std::f64::consts::PI;

use ndarray::Array2;
use num_complex::Complex;

use super::{Field, SimulationParametersArgs};

/// Fz_rad=-∂U/∂z with U = 2πR³Ũ
pub fn compute_radiation_force_field(
    relative_acoustic_radiation_field: &Field,
    sp: SimulationParametersArgs,
) -> Field {
    let SimulationParametersArgs {
        z_max,
        z_min,
        sphere_radius,
        ..
    } = sp;

    // compute acoustic radiation for a small sphere of radius r
    let acoustic_radiation_field =
        relative_acoustic_radiation_field.map(|ar| 2.0 * PI * sphere_radius.powi(3) * ar);

    let mut fz_rad: Array2<Complex<f64>> = Array2::zeros(acoustic_radiation_field.raw_dim());
    let (row, col) = match *acoustic_radiation_field.shape() {
        [r, c] => (r, c),
        _ => panic!(),
    };

    // compute ∂U/∂z
    for i in 0..row {
        for j in 0..col {
            if i == 0 {
                // forward
                let (z1, z2) = (
                    z_min + (i as f64) * (z_max - z_min) / (row - 1) as f64,
                    z_min + ((i + 1) as f64) * (z_max - z_min) / (row - 1) as f64,
                );
                fz_rad[[i, j]] = (acoustic_radiation_field[[i + 1, j]]
                    - acoustic_radiation_field[[i, j]])
                    / (z2 - z1)
            } else if i == row - 1 {
                // backward
                let (z1, z2) = (
                    z_min + ((i - 1) as f64) * (z_max - z_min) / (row - 1) as f64,
                    z_min + (i as f64) * (z_max - z_min) / (row - 1) as f64,
                );
                fz_rad[[i, j]] = (acoustic_radiation_field[[i, j]]
                    - acoustic_radiation_field[[i - 1, j]])
                    / (z2 - z1)
            } else {
                // Centered
                let (z1, z2) = (
                    z_min + ((i - 1) as f64) * (z_max - z_min) / (row - 1) as f64,
                    z_min + ((i + 1) as f64) * (z_max - z_min) / (row - 1) as f64,
                );
                fz_rad[[i, j]] = (acoustic_radiation_field[[i + 1, j]]
                    - acoustic_radiation_field[[i - 1, j]])
                    / (2.0 * (z2 - z1))
            }
        }
    }

    // transforms it to -∂U/∂z
    fz_rad.map(|du| -du)
}
