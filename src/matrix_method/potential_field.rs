use std::f64::consts::PI;

use ndarray::Array2;
use num_complex::Complex;

use crate::matrix_method::C;

use super::{SimulationParameters, I, RHO};

fn velocity_potential(
    pressure: &Array2<Complex<f64>>,
    sp: impl Into<SimulationParameters>,
) -> Array2<Complex<f64>> {
    let sp: SimulationParameters = sp.into();
    pressure.map(|p| -p / (I * sp.omega * RHO))
}

pub fn velocity_field(
    pressure: &Array2<Complex<f64>>,
    sp: impl Into<SimulationParameters>,
) -> (Array2<Complex<f64>>, Array2<Complex<f64>>) {
    let sp: SimulationParameters = sp.into();
    let velocity_potential = velocity_potential(pressure, sp.clone());
    // compute the gradient of the velocity potential
    let mut vf_x: Array2<Complex<f64>> = Array2::zeros(velocity_potential.raw_dim());
    let mut vf_z: Array2<Complex<f64>> = Array2::zeros(velocity_potential.raw_dim());
    let (row, col) = match *velocity_potential.shape() {
        [r, c] => (r, c),
        _ => panic!(),
    };

    let SimulationParameters {
        z_max,
        z_min,
        x_min,
        x_max,
        ..
    } = sp;
    for i in 0..row {
        for j in 0..col {
            // for z direction
            {
                if i == 0 {
                    // forward
                    let (z1, z2) = (
                        z_min + (i as f64) * (z_max - z_min) / row as f64,
                        z_min + ((i + 1) as f64) * (z_max - z_min) / row as f64,
                    );
                    vf_z[[i, j]] =
                        (velocity_potential[[i + 1, j]] - velocity_potential[[i, j]]) / (z2 - z1)
                } else if i == row - 1 {
                    // backward
                    let (z1, z2) = (
                        z_min + ((i - 1) as f64) * (z_max - z_min) / row as f64,
                        z_min + (i as f64) * (z_max - z_min) / row as f64,
                    );
                    vf_z[[i, j]] =
                        (velocity_potential[[i, j]] - velocity_potential[[i - 1, j]]) / (z2 - z1)
                } else {
                    // Centered
                    let (z1, z2) = (
                        z_min + ((i - 1) as f64) * (z_max - z_min) / row as f64,
                        z_min + ((i + 1) as f64) * (z_max - z_min) / row as f64,
                    );
                    vf_z[[i, j]] = (velocity_potential[[i + 1, j]] - velocity_potential[[i - 1, j]])
                        / (2.0 * (z2 - z1))
                }
            }

            // for x direction: not used here
            {
                if j == 0 {
                    // forward
                    let (x1, x2) = (
                        x_min + (j as f64) * (x_max - x_min) / col as f64,
                        x_min + ((j + 1) as f64) * (x_max - x_min) / col as f64,
                    );
                    vf_x[[i, j]] =
                        (velocity_potential[[i, j + 1]] - velocity_potential[[i, j]]) / (x2 - x1)
                } else if j == col - 1 {
                    // backward
                    let (x1, x2) = (
                        x_min + ((j - 1) as f64) * (x_max - x_min) / col as f64,
                        x_min + (j as f64) * (x_max - x_min) / col as f64,
                    );
                    vf_x[[i, j]] =
                        (velocity_potential[[i, j]] - velocity_potential[[i, j - 1]]) / (x2 - x1);
                } else {
                    // Centered
                    let (x1, x2) = (
                        x_min + ((j - 1) as f64) * (x_max - x_min) / col as f64,
                        x_min + ((j + 1) as f64) * (x_max - x_min) / col as f64,
                    );
                    vf_x[[i, j]] = (velocity_potential[[i, j + 1]] - velocity_potential[[i, j - 1]])
                        / (2.0 * (x2 - x1))
                }
            }
        }
    }

    (vf_x, vf_z)
}

pub fn acoustic_radiation_potential(
    pressure: &Array2<Complex<f64>>,
    sp: impl Into<SimulationParameters>,
    sphere_radius: Option<f64>,
) -> Array2<f64> {
    let vf = velocity_field(pressure, sp);
    assert_eq!(pressure.shape(), vf.0.shape());

    // compute the mean square amplitude of the pressure & velocity field
    let (row, col) = match *pressure.shape() {
        [r, c] => (r, c),
        _ => panic!(),
    };

    let mut relative_acoustic_potential: Array2<f64> = Array2::zeros(pressure.raw_dim());
    for i in 0..row {
        for j in 0..col {
            let (pmsa, vmsa) = (
                pressure[[i, j]].norm_sqr(),
                vf.0[[i, j]].norm_sqr() + vf.1[[i, j]].norm_sqr(),
            );
            relative_acoustic_potential[[i, j]] =
                (pmsa / (3.0 * RHO.re * C.re.powi(2))) - ((RHO.re * vmsa) / 2.0);
        }
    }

    match sphere_radius {
        Some(r) => relative_acoustic_potential.map(|ap| 2.0 * PI * r.powi(3) * ap),
        None => relative_acoustic_potential,
    }
}
