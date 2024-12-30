use std::fs;

use matrix_method::{
    potential_field::compute_relative_potential_field, pressure_field::compute_pressure_field,
    radiation_force::compute_radiation_force_field, FieldType, SimulationParametersArgs, MM,
};
use num_complex::Complex;
use plot::{field::plot_field, graph::graph_field};

mod matrix_method;
mod plot;

fn main() {
    let _ = fs::create_dir_all("./out/field");
    let _ = fs::create_dir_all("./out/graph");
    let _ = fs::create_dir_all("./tmp");

    let simulation_parameters = SimulationParametersArgs {
        x_min: -50.0 * MM,
        x_max: 50.0 * MM,
        z_min: 0.0 * MM,
        z_max: 50.0 * MM,
        nb_of_reflection: 4,
        disc: 0.4 * MM,
        radius: 15.0 * MM,
        hole_radius: 2.0 * MM,
        freq: Complex::new(56000.0, 0.0),
        u_0: 0.0000060,
        inclination: 17.0,
        curvature: 0.0,
        sphere_radius: 0.1 * MM, // wavelength=6.1mm
    };
    // let simulation_parameters = SimulationParametersArgs::default();

    // Pressure plotting
    let pressure = compute_pressure_field(simulation_parameters);
    plot_field(
        (&pressure, FieldType::Pressure),
        (None, 4.0),
        simulation_parameters,
        "./out/field/pressure_field.png",
    )
    .expect("Failed to draw pressure field");
    graph_field(
        (&pressure, FieldType::Pressure),
        (0.0, Some(0.005)),
        simulation_parameters,
        "./out/graph/pressure_graph.png",
    )
    .expect("Failed to graph pressure field");

    // Potential plotting
    let relative_potential_field =
        compute_relative_potential_field(&pressure, simulation_parameters);
    plot_field(
        (&relative_potential_field, FieldType::RadiationPotential),
        (None, 50.0),
        simulation_parameters,
        "./out/field/acoustic_radiation_potential.png",
    )
    .expect("Failed to draw acoustic field");

    // Force plotting
    let radiation_force_field =
        compute_radiation_force_field(&relative_potential_field, simulation_parameters);
    plot_field(
        (&radiation_force_field, FieldType::RadiationForce),
        (None, 6.0),
        simulation_parameters,
        "./out/field/radiation_force.png",
    )
    .expect("Failed to draw force field");
    graph_field(
        (&radiation_force_field, FieldType::RadiationForce),
        (0.0, Some(0.005)),
        simulation_parameters,
        "./out/graph/force_graph.png",
    )
    .expect("Failed to graph force field");
}
