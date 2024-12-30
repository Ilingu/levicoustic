use ndarray::{array, Array1, Array2};

use super::*;

#[allow(non_snake_case)]
/// global_simulation_params are non transducer params : grid dimension, nb reflection, disc, freq, u_0
///
/// simulations_params are the transducers params: offset, radius, hole, phase, inclination
pub fn compute_pressure_fields(
    global_simulation_params: impl Into<SimulationParameters>,
    simulations_params: Vec<impl Into<SimulationParameters>>,
) -> Field {
    let SimulationParameters {
        x_min,
        x_max,
        z_min,
        z_max,

        nb_reflection,
        disc,
        u_0,

        reflector_area,

        omega,
        wavenumber,
        e,
        d,
        ..
    } = global_simulation_params.into();
    let cells = 100.0; // Number of discrete cells
    let si = reflector_area / (cells * 4.0); // Unit cell area for reflector

    let L: usize = ((x_max - x_min) / disc) as usize;
    let Z: usize = ((z_max - z_min) / disc) as usize;
    let M: usize = L * Z;

    let mut t_rm: Option<Array2<Complex<f64>>> = None;
    let mut transfer_matrices: Vec<[Array2<Complex<f64>>; 4]> =
        Vec::with_capacity(simulations_params.len());
    for sp in simulations_params {
        let SimulationParameters {
            offset,
            radius: R,
            transducer_area,
            phase,
            inclination,
            hole_area,
            ..
        }: SimulationParameters = sp.into();
        // convert deg to rad
        let (inclination, phase) = (inclination.to_radians(), phase.to_radians());

        // Array indices
        let N: usize = ((2.0 * R) / disc) as usize;

        // Creation of grid, indices, and distance arrays
        let x: Array1<f64> = Array1::linspace(x_min, x_max, L);
        let z: Array1<f64> = Array1::linspace(z_min, z_max, Z);
        let transducer: Array1<(f64, f64)> = Array1::linspace(-R, R, N)
            .map(|x| rotation_matrix(inclination).dot(&array![[*x], [0.0]]))
            .map(|x| {
                (
                    x[[0, 0]] + offset,
                    x[[1, 0]] + R * inclination.abs().sin() + z_max,
                )
            });

        // Create zeroed distance arrays in memory
        let mut r_nm: Array2<f64> = Array2::zeros((N, M));
        let mut r_in: Array2<f64> = Array2::zeros((N, L));

        // Calculate distance array r_nm
        for i in 0..N {
            let mut q = 0;
            for j in 0..L {
                for k in 0..Z {
                    let (tx, tz) = transducer[i];
                    r_nm[[i, q]] = ((tx - x[j]).powi(2) + (tz - z[k]).powi(2)).sqrt();
                    q += 1;
                }
            }
        }

        // Calculate distance array r_in
        for i in 0..N {
            for j in 0..L {
                let (tx, tz) = transducer[i];
                r_in[[i, j]] = ((tx - x[j]).powi(2) + (tz - z_min).powi(2)).sqrt();
            }
        }

        let r_ni = r_in.t();

        // Creation of cell discretization and transfer matrices calculations
        let sn = transducer_area / cells; // Unit cell area for transducer
        let sh = hole_area / cells; // Unit cell area for transducer hole

        // Create zeroed transfer matrices in memory
        let mut t_tm: Array2<Complex<f64>> = Array2::zeros((N, M));
        let mut t_tr: Array2<Complex<f64>> = Array2::zeros((N, L));
        let mut t_rt: Array2<Complex<f64>> = Array2::zeros((L, N));

        // Calculate transfer matrices
        for i in 0..N {
            for j in 0..M {
                t_tm[[i, j]] = ((sn * (-I * wavenumber * r_nm[[i, j]]).exp()) / r_nm[[i, j]])
                    - ((sh * (-I * wavenumber * r_nm[[i, j]]).exp()) / r_nm[[i, j]]);
            }
            for k in 0..L {
                t_tr[[i, k]] = ((sn * (-I * wavenumber * r_in[[i, k]]).exp()) / r_in[[i, k]])
                    - ((sh * (-I * wavenumber * r_in[[i, k]]).exp()) / r_in[[i, k]])
            }
        }

        for i in 0..L {
            for j in 0..N {
                t_rt[[i, j]] = (si * (-I * wavenumber * r_ni[[i, j]]).exp()) / r_ni[[i, j]]
            }
        }

        // Rotation for column notation
        let t_tm = t_tm.t().to_owned();
        let t_tr = t_tr.t().to_owned();
        let t_rt = t_rt.t().to_owned();

        // t_rm special case - it's a global var that is computed once
        if t_rm.is_none() {
            let mut r_im: Array2<f64> = Array2::zeros((L, M));
            // Calculate distance array r_im
            for i in 0..L {
                let mut q = 0;
                for j in 0..L {
                    for k in 0..Z {
                        r_im[[i, q]] = ((x[i] - x[j]).powi(2) + (z_min - z[k]).powi(2)).sqrt();
                        q += 1;
                    }
                }
            }
            let mut t_rm_l: Array2<Complex<f64>> = Array2::zeros((L, M));
            for i in 0..L {
                for j in 0..M {
                    if r_im[[i, j]] == 0.0 {
                        continue;
                    }
                    t_rm_l[[i, j]] = (si * (-I * wavenumber * r_im[[i, j]]).exp()) / r_im[[i, j]]
                }
            }
            let t_rm_l = t_rm_l.t().to_owned();
            t_rm = Some(t_rm_l)
        }

        // Boundary conditions of transducer
        #[allow(non_snake_case)]
        let mut U_n: Array2<Complex<f64>> = Array2::zeros((N, 1));
        for i in 0..N {
            U_n[[i, 0]] = u_0 * (I * (omega + phase)).exp()
        }

        transfer_matrices.push([t_tm, t_tr, t_rt, U_n]);
    }

    let t_rm = t_rm.expect("No possible in theory");
    let [t_ntm, t_ntr, t_nrt] = transfer_matrices.iter().skip(1).fold(
        {
            let [u, d, t, _] = transfer_matrices[0].clone();
            [u, d, t]
        },
        |[t_ntm, t_ntr, t_nrt], [t_tm_i, t_tr_i, t_rt_i, _]| {
            [
                add_matrices(vec![&t_ntm, t_tm_i]),
                add_matrices(vec![&t_ntr, t_tr_i]),
                add_matrices(vec![&t_nrt, t_rt_i]),
            ]
        },
    );

    let mut pressure_fields: Vec<Array2<Complex<f64>>> =
        Vec::with_capacity(transfer_matrices.len());
    for [t_tm_i, t_tr_i, _, U_i] in transfer_matrices {
        // Calculation of pressure, where each reflection is an order of approximation
        let mut pressure_base = Complex::new(d, 0.0) * t_tm_i.dot(&U_i);
        for n in 1..=nb_reflection {
            pressure_base = pressure_base
                + (d * e.powu(n as u32))
                    * nth_reflection(n, [&t_ntm, &t_ntr, &t_nrt, &t_rm, &t_tr_i.dot(&U_i)], true);
        }

        let pf = pressure_base
            .into_shape_with_order((L, Z))
            .expect("Failed to reshape pressure field")
            .t()
            .to_owned();
        pressure_fields.push(pf);
    }

    add_owned_matrices(pressure_fields) // global pressure field
}

/// compute the pressure field for each points thanks to the matrix method - **only for 1 transducer**
///
/// Output a **complex** field, only the real part of the field correspond to the actual pressure field
#[allow(non_snake_case)]
pub fn compute_pressure_field(simulation_params: impl Into<SimulationParameters>) -> Field {
    let SimulationParameters {
        x_min,
        x_max,
        z_min,
        z_max,

        nb_reflection,
        disc,

        offset,
        radius: R,
        transducer_area,
        hole_area,
        phase,
        u_0,
        inclination,

        reflector_area,

        omega,
        wavenumber,
        e,
        d,
        ..
    }: SimulationParameters = simulation_params.into();
    // convert deg to rad
    let (inclination, phase) = (inclination.to_radians(), phase.to_radians());

    // Array indices
    let N: usize = ((2.0 * R) / disc) as usize;
    let L: usize = ((x_max - x_min) / disc) as usize;
    let Z: usize = ((z_max - z_min) / disc) as usize;
    let M: usize = L * Z;

    // Creation of grid, indices, and distance arrays
    let x: Array1<f64> = Array1::linspace(x_min, x_max, L);
    let z: Array1<f64> = Array1::linspace(z_min, z_max, Z);
    let transducer: Array1<(f64, f64)> = Array1::linspace(-R, R, N)
        .map(|x| rotation_matrix(inclination).dot(&array![[*x], [0.0]]))
        .map(|x| {
            (
                x[[0, 0]] + offset,
                x[[1, 0]] + R * inclination.abs().sin() + z_max,
            )
        });

    // Create zeroed distance arrays in memory
    let mut r_nm: Array2<f64> = Array2::zeros((N, M));
    let mut r_im: Array2<f64> = Array2::zeros((L, M));
    let mut r_in: Array2<f64> = Array2::zeros((N, L));

    // Calculate distance array r_nm
    for i in 0..N {
        let mut q = 0;
        for j in 0..L {
            for k in 0..Z {
                let (tx, tz) = transducer[i];
                r_nm[[i, q]] = ((tx - x[j]).powi(2) + (tz - z[k]).powi(2)).sqrt();
                q += 1;
            }
        }
    }

    // Calculate distance array r_im
    for i in 0..L {
        let mut q = 0;
        for j in 0..L {
            for k in 0..Z {
                r_im[[i, q]] = ((x[i] - x[j]).powi(2) + (z_min - z[k]).powi(2)).sqrt();
                q += 1;
            }
        }
    }

    // Calculate distance array r_in
    for i in 0..N {
        for j in 0..L {
            let (tx, tz) = transducer[i];
            r_in[[i, j]] = ((tx - x[j]).powi(2) + (tz - z_min).powi(2)).sqrt();
        }
    }

    let r_ni = r_in.t();

    // Creation of cell discretization and transfer matrices calculations
    let cells = 100.0; // Number of discrete cells
    let sn = transducer_area / cells; // Unit cell area for transducer
    let si = reflector_area / (cells * 4.0); // Unit cell area for reflector
    let sh = hole_area / cells; // Unit cell area for transducer hole

    // Create zeroed transfer matrices in memory
    let mut t_tm: Array2<Complex<f64>> = Array2::zeros((N, M));
    let mut t_tr: Array2<Complex<f64>> = Array2::zeros((N, L));
    let mut t_rt: Array2<Complex<f64>> = Array2::zeros((L, N));
    let mut t_rm: Array2<Complex<f64>> = Array2::zeros((L, M));

    // Calculate transfer matrices
    for i in 0..N {
        for j in 0..M {
            t_tm[[i, j]] = ((sn * (-I * wavenumber * r_nm[[i, j]]).exp()) / r_nm[[i, j]])
                - ((sh * (-I * wavenumber * r_nm[[i, j]]).exp()) / r_nm[[i, j]]);
        }
        for k in 0..L {
            t_tr[[i, k]] = ((sn * (-I * wavenumber * r_in[[i, k]]).exp()) / r_in[[i, k]])
                - ((sh * (-I * wavenumber * r_in[[i, k]]).exp()) / r_in[[i, k]])
        }
    }

    for i in 0..L {
        for j in 0..N {
            t_rt[[i, j]] = (si * (-I * wavenumber * r_ni[[i, j]]).exp()) / r_ni[[i, j]]
        }
    }

    for i in 0..L {
        for j in 0..M {
            if r_im[[i, j]] == 0.0 {
                continue;
            }
            t_rm[[i, j]] = (si * (-I * wavenumber * r_im[[i, j]]).exp()) / r_im[[i, j]]
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
        U[[i, 0]] = u_0 * (I * (omega + phase)).exp()
    }

    // Calculation of pressure, where each reflection is an order of approximation
    let mut pressure_base = Complex::new(d, 0.0) * t_tm.dot(&U);
    for n in 1..=nb_reflection {
        pressure_base = pressure_base
            + (d * e.powu(n as u32)) * nth_reflection(n, [&t_tm, &t_tr, &t_rt, &t_rm, &U], false);
    }

    pressure_base
        .into_shape_with_order((L, Z))
        .expect("Failed to reshape pressure field")
        .t()
        .to_owned()
}

/* HELPERS */

fn rotation_matrix(theta: f64) -> Array2<f64> {
    array![[theta.cos(), -theta.sin()], [theta.sin(), theta.cos()]]
}

/// list all the transfer matrices necessary for the nth reflection of the wave and return their product
fn nth_reflection(
    n: u8,
    [tm, tr, rt, rm, u]: [&Array2<Complex<f64>>; 5],
    multiple_simu: bool,
) -> Array2<Complex<f64>> {
    let mut list_of_reflection = Vec::with_capacity(n as usize + 2);
    list_of_reflection.push(match n % 2 == 0 {
        true => tm,
        false => rm,
    });
    for i in 1..=if multiple_simu { n - 1 } else { n } {
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
