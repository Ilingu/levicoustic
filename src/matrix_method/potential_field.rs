use ndarray::Array2;
use num_complex::Complex;

use crate::matrix_method::C;

use super::{Field, SimulationParameters, I, RHO};

fn velocity_potential(pressure: &Field, sp: impl Into<SimulationParameters>) -> Field {
    let sp: SimulationParameters = sp.into();
    pressure.map(|p| -p / (I * sp.omega * RHO))
}

pub fn compute_velocity_field(
    pressure: &Field,
    sp: impl Into<SimulationParameters>,
) -> (Field, Field) {
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
                        z_min + (i as f64) * (z_max - z_min) / (row - 1) as f64,
                        z_min + ((i + 1) as f64) * (z_max - z_min) / (row - 1) as f64,
                    );
                    vf_z[[i, j]] =
                        (velocity_potential[[i + 1, j]] - velocity_potential[[i, j]]) / (z2 - z1)
                } else if i == row - 1 {
                    // backward
                    let (z1, z2) = (
                        z_min + ((i - 1) as f64) * (z_max - z_min) / (row - 1) as f64,
                        z_min + (i as f64) * (z_max - z_min) / (row - 1) as f64,
                    );
                    vf_z[[i, j]] =
                        (velocity_potential[[i, j]] - velocity_potential[[i - 1, j]]) / (z2 - z1)
                } else {
                    // Centered
                    let (z1, z2) = (
                        z_min + ((i - 1) as f64) * (z_max - z_min) / (row - 1) as f64,
                        z_min + ((i + 1) as f64) * (z_max - z_min) / (row - 1) as f64,
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
                        x_min + (j as f64) * (x_max - x_min) / (col - 1) as f64,
                        x_min + ((j + 1) as f64) * (x_max - x_min) / (col - 1) as f64,
                    );
                    vf_x[[i, j]] =
                        (velocity_potential[[i, j + 1]] - velocity_potential[[i, j]]) / (x2 - x1)
                } else if j == col - 1 {
                    // backward
                    let (x1, x2) = (
                        x_min + ((j - 1) as f64) * (x_max - x_min) / (col - 1) as f64,
                        x_min + (j as f64) * (x_max - x_min) / (col - 1) as f64,
                    );
                    vf_x[[i, j]] =
                        (velocity_potential[[i, j]] - velocity_potential[[i, j - 1]]) / (x2 - x1);
                } else {
                    // Centered
                    let (x1, x2) = (
                        x_min + ((j - 1) as f64) * (x_max - x_min) / (col - 1) as f64,
                        x_min + ((j + 1) as f64) * (x_max - x_min) / (col - 1) as f64,
                    );
                    vf_x[[i, j]] = (velocity_potential[[i, j + 1]] - velocity_potential[[i, j - 1]])
                        / (2.0 * (x2 - x1))
                }
            }
        }
    }

    (vf_x, vf_z)
}

/// compute the relative acoustic radiation potential Ũ in Gor’kov theory
///
/// The output type is a complex field BUT the field itself is 100% real
pub fn compute_relative_potential_field(
    pressure: &Field,
    sp: impl Into<SimulationParameters>,
) -> Field {
    let vf = compute_velocity_field(pressure, sp);
    assert_eq!(pressure.shape(), vf.0.shape());

    // compute the mean square amplitude of the pressure & velocity field
    let (row, col) = match *pressure.shape() {
        [r, c] => (r, c),
        _ => panic!(),
    };

    let mut relative_acoustic_potential: Array2<f64> = Array2::zeros(pressure.raw_dim());
    for i in 0..row {
        for j in 0..col {
            let mut nb_of_neighbor: usize = 0;
            let (mut p_sum, mut v_sum) = (0.0, 0.0);

            const SQUARE_SIZE: isize = 1; // only odd positive number
            for k in 0..SQUARE_SIZE {
                for l in 0..SQUARE_SIZE {
                    let (k, l) = (k - (SQUARE_SIZE - 1) / 2, l - (SQUARE_SIZE - 1) / 2);
                    let (n_row, n_col) = (i as isize + k, j as isize + l);
                    if n_row < 0
                        || n_row > row as isize - 1
                        || n_col < 0
                        || n_col > col as isize - 1
                    {
                        continue;
                    }

                    nb_of_neighbor += 1;
                    let (n_row, n_col) = (n_row as usize, n_col as usize);
                    // p_sum += pressure[[n_row, n_col]].re.powi(2);
                    // v_sum += vf.0[[n_row, n_col]].re.powi(2) + vf.1[[n_row, n_col]].re.powi(2);
                    p_sum += pressure[[n_row, n_col]].norm_sqr();
                    v_sum += vf.0[[n_row, n_col]].norm_sqr() + vf.1[[n_row, n_col]].norm_sqr();
                }
            }
            let (pmsa, vmsa) = (p_sum / nb_of_neighbor as f64, v_sum / nb_of_neighbor as f64);
            relative_acoustic_potential[[i, j]] =
                (pmsa / (3.0 * RHO.re * C.re.powi(2))) - ((RHO.re * vmsa) / 2.0);
        }
    }

    relative_acoustic_potential.map(|ap| Complex::new(*ap, 0.0))
}
