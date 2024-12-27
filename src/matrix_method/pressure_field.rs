use ndarray::{Array1, Array2};

use super::*;

// Array indices
const N: usize = ((2.0 * R) / DISC) as usize;
const L: usize = ((X_MAX - X_MIN) / DISC) as usize;
const Z: usize = ((Z_MAX - Z_MIN) / DISC) as usize;
const M: usize = L * Z;

/// compute the pressure field for each points.
/// The greater the number of reflections, the more precise the result will be. However it comes at a computational cost
pub fn compute_pressure_field(nb_of_reflection: u8) -> Array2<Complex<f64>> {
    // Creation of grid, indices, and distance arrays
    let x: Array1<f64> = Array1::linspace(X_MIN, X_MAX, L);
    let z: Array1<f64> = Array1::linspace(Z_MIN, Z_MAX, Z);
    let transducer: Array1<f64> = Array1::linspace(-R, R, N);

    // Create zeroed distance arrays in memory
    let mut r_nm: Array2<f64> = Array2::zeros((N, M));
    let mut r_im: Array2<f64> = Array2::zeros((L, M));
    let mut r_in: Array2<f64> = Array2::zeros((N, L));

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
    let t_tm = t_tm.t().to_owned();
    let t_tr = t_tr.t().to_owned();
    let t_rt = t_rt.t().to_owned();
    let t_rm = t_rm.t().to_owned();

    // Boundary conditions of transducer
    #[allow(non_snake_case)]
    let mut U: Array2<Complex<f64>> = Array2::zeros((N, 1));
    for i in 0..N {
        U[[i, 0]] = U_0 * (I * OMEGA).exp()
    }

    // Calculation of pressure, where each line is an order of approximation
    let mut pressure_base = D * t_tm.dot(&U);
    for n in 1..=nb_of_reflection {
        pressure_base = pressure_base
            + (D * E.powu(n as u32)) * nth_reflection(n, [&t_tm, &t_tr, &t_rt, &t_rm, &U]);
    }

    pressure_base
        .into_shape_with_order((L, Z))
        .expect("Failed to reshape pressure field")
        .t()
        .to_owned()
}

fn nth_reflection(n: u8, [tm, tr, rt, rm, u]: [&Array2<Complex<f64>>; 5]) -> Array2<Complex<f64>> {
    let mut list_of_reflection = Vec::with_capacity(n as usize + 2);
    list_of_reflection.push(match n % 2 == 0 {
        true => tm,
        false => rm,
    });
    for i in 1..=n {
        let r = match (n % 2 == 0, i % 2 == 0) {
            (true, true) => tr,
            (true, false) => rt,
            (false, true) => rt,
            (false, false) => tr,
        };
        list_of_reflection.push(r);
    }
    list_of_reflection.push(u);
    multiple_matrix_product(list_of_reflection)
}
