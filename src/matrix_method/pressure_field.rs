use ndarray::{Array, Array1, Array2};

use super::*;

// Array indices
const N: usize = ((2.0 * R) / DISC) as usize;
const L: usize = ((X_MAX - X_MIN) / DISC) as usize;
const Z: usize = ((Z_MAX - Z_MIN) / DISC) as usize;
const M: usize = L * Z;

pub fn compute_pressure_field() -> Array2<Complex<f64>> {
    // Creation of grid, indices, and distance arrays
    let x: Array1<f64> = Array::linspace(X_MIN, X_MAX, L);
    let z: Array1<f64> = Array::linspace(Z_MIN, Z_MAX, Z);
    let transducer: Array1<f64> = Array::linspace(-R, R, N);

    // Create zeroed distance arrays in memory
    let mut r_nm: Array2<f64> = Array::zeros((N, M));
    let mut r_im: Array2<f64> = Array::zeros((L, M));
    let mut r_in: Array2<f64> = Array::zeros((N, L));

    // Calculate distance array r_nm
    for i in 0..N {
        let mut q = 0;
        for j in 0..L {
            for k in 0..Z {
                r_nm[[i, q]] = ((transducer[i] - x[j]).powi(2) + (Z_MAX - z[k]).powi(2)).sqrt();
                q += 1;
            }
        }
    }

    // Calculate distance array r_im
    for i in 0..L {
        let mut q = 0;
        for j in 0..L {
            for k in 0..Z {
                r_im[[i, q]] = ((x[i] - x[j]).powi(2) + (Z_MIN - z[k]).powi(2)).sqrt();
                q += 1;
            }
        }
    }

    // Calculate distance array r_in
    for i in 0..N {
        for j in 0..L {
            r_in[[i, j]] = ((transducer[i] - x[j]).powi(2) + (Z_MAX - Z_MIN).powi(2)).sqrt();
        }
    }

    let r_ni = r_in.t();

    // Creation of cell discretization and transfer matrices calculations
    let cells: Complex<f64> = Complex::new(100.0, 0.0); // Number of discrete cells
    let sn = TRANSDUCER_AREA / cells; // Unit cell area for transducer
    let si = REFLECTOR_AREA / (cells * 4.0); // Unit cell area for reflector
    let sh = HOLE_AREA / cells; // Unit cell area for transducer hole

    // Create zeroed transfer matrices in memory
    let mut t_tm: Array2<Complex<f64>> = Array2::zeros((N, M));
    let mut t_tr: Array2<Complex<f64>> = Array2::zeros((N, L));
    let mut t_rt: Array2<Complex<f64>> = Array2::zeros((L, N));
    let mut t_rm: Array2<Complex<f64>> = Array2::zeros((L, M));

    // Calculate transfer matrices
    const I: Complex<f64> = Complex::I;

    for i in 0..N {
        for j in 0..M {
            t_tm[[i, j]] = ((sn * (-I * WAVENUMBER * r_nm[[i, j]]).exp()) / r_nm[[i, j]])
                - ((sh * (-I * WAVENUMBER * r_nm[[i, j]]).exp()) / r_nm[[i, j]]);
        }
        for k in 0..L {
            t_tr[[i, k]] = ((sn * (-I * WAVENUMBER * r_in[[i, k]]).exp()) / r_in[[i, k]])
                - ((sh * (-I * WAVENUMBER * r_in[[i, k]]).exp()) / r_in[[i, k]])
        }
    }

    for i in 0..L {
        for j in 0..N {
            t_rt[[i, j]] = (si * (-I * WAVENUMBER * r_ni[[i, j]]).exp()) / r_ni[[i, j]]
        }
    }

    for i in 0..L {
        for j in 0..M {
            if r_im[[i, j]] == 0.0 {
                continue;
            }
            t_rm[[i, j]] = (si * (-I * WAVENUMBER * r_im[[i, j]]).exp()) / r_im[[i, j]]
        }
    }

    // Rotation for column notation
    let t_tm = t_tm.t();
    let t_tr = t_tr.t();
    let t_rt = t_rt.t();
    let t_rm = t_rm.t();

    // Boundary conditions of transducer
    #[allow(non_snake_case)]
    let mut U: Array2<Complex<f64>> = Array2::zeros((N, 1));
    for i in 0..N {
        U[[i, 0]] = U_0 * (I * OMEGA).exp()
    }

    // Calculation of pressure, where each line is an order of approximation
    let pressure_base = D * t_tm.dot(&U)
        + (D * E) * (t_rm.dot(&t_tr)).dot(&U)
        + (D * E.powi(2)) * (t_tm.dot(&t_rt)).dot(&(t_tr.dot(&U)))
        + (D * E.powi(3)) * t_rm.dot(&((t_tr.dot(&t_rt)).dot(&(t_tr.dot(&U)))))
        + (D * E.powi(4)) * (t_tm.dot(&t_rt)).dot(&((t_tr.dot(&t_rt)).dot(&(t_tr.dot(&U)))));

    pressure_base
        .into_shape_with_order((L, Z))
        .expect("Failed to reshape pressure field")
        .t()
        .to_owned()
}
