use std::fs;

use matrix_method::{
    potential_field::acoustic_radiation_potential, pressure_field, SimulationParametersArgs, MM,
};
use num_complex::Complex;
use plot::plot_pressure_field;

mod matrix_method;
mod plot;

fn main() {
    let _ = fs::create_dir_all("./out");
    let _ = fs::create_dir_all("./tmp");

    let simulation_parameters = SimulationParametersArgs {
        x_min: -20.0 * MM,
        x_max: 20.0 * MM,
        z_min: -18.0 * MM,
        z_max: 0.0 * MM,
        nb_of_reflection: 2,
        disc: 0.1 * MM,
        radius: 5.0 * MM,
        hole_radius: 0.0,
        freq: Complex::new(20_000.0, 0.0),
        u_0: 0.000001,
        saturation: 3.0,
    };
    // let simulation_parameters = SimulationParametersArgs::default();

    let pressure = pressure_field::compute_pressure_field(simulation_parameters);
    plot_pressure_field(&pressure, simulation_parameters, "./out/pressure_field.png")
        .expect("Failed to draw pf");

    let acoustic_radiation_field =
        acoustic_radiation_potential(&pressure, simulation_parameters, None);

    plot_pressure_field(
        &acoustic_radiation_field.map(|&a| Complex::new(a, 0.0)),
        simulation_parameters,
        "./out/acoustic_field.png",
    )
    .expect("Failed to draw pf");
}
