use std::fs;

use matrix_method::{pressure_field, SimulationParametersArgs};
use plot::plot_pressure_field;

mod matrix_method;
mod plot;

fn main() {
    let _ = fs::create_dir_all("./out");
    let _ = fs::create_dir_all("./tmp");

    let simulation_parameters = SimulationParametersArgs::default();

    let pressure = pressure_field::compute_pressure_field(simulation_parameters);
    plot_pressure_field(&pressure, simulation_parameters, "./out/pressure_field.png")
        .expect("Failed to draw pf");
}
