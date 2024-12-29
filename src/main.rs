use std::fs;

use matrix_method::{
    potential_field::{compute_acoustic_radiation_field, compute_velocity_field},
    pressure_field, Field, FieldType, SimulationParametersArgs, MM,
};
use ndarray::Array2;
use num_complex::Complex;
use plot::plot_field;

mod matrix_method;
mod plot;

fn main() {
    let _ = fs::create_dir_all("./out");
    let _ = fs::create_dir_all("./tmp");

    let simulation_parameters = SimulationParametersArgs {
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
    };
    // let simulation_parameters = SimulationParametersArgs::default();

    let pressure = pressure_field::compute_pressure_field(simulation_parameters);
    plot_field(
        (&pressure, FieldType::Pressure),
        (None, 4.0),
        simulation_parameters,
        "./out/pressure_field.png",
    )
    .expect("Failed to draw pressure field");

    {
        let velocity_field = compute_velocity_field(&pressure, simulation_parameters);
        plot_field(
            (&velocity_field.0, FieldType::Velocity),
            (None, 5000.0),
            simulation_parameters,
            "./out/velocity_field_x.png",
        )
        .expect("Failed to draw velocity field");
        plot_field(
            (&velocity_field.1, FieldType::Velocity),
            (None, 5.0),
            simulation_parameters,
            "./out/velocity_field_z.png",
        )
        .expect("Failed to draw velocity field");
    }

    let acoustic_radiation_field =
        compute_acoustic_radiation_field(&pressure, simulation_parameters, None);
    plot_field(
        (&acoustic_radiation_field, FieldType::RadiationPotential),
        (Some(-75.0..75.0), 1.0),
        simulation_parameters,
        "./out/acoustic_radiation_potential.png",
    )
    .expect("Failed to draw acoustic field");
}
