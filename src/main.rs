use std::fs;

use matrix_method::{
    potential_field::compute_relative_potential_field,
    pressure_field::{compute_pressure_field, compute_pressure_fields},
    radiation_force::compute_radiation_force_field,
    Field, FieldType, SimulationParametersArgs, MM,
};
use plot::{field::plot_field, graph::graph_field};

mod matrix_method;
mod plot;

fn main() {
    let _ = fs::create_dir_all("./out/field");
    let _ = fs::create_dir_all("./out/graph");
    let _ = fs::create_dir_all("./tmp");

    let spg = SimulationParametersArgs {
        x_min: -30.0 * MM,
        x_max: 30.0 * MM,
        z_min: 0.0 * MM,
        z_max: 25.0 * MM,

        nb_reflection: 4,
        disc: 0.25 * MM,

        freq: 37900.0,
        u_0: 0.0000060,

        sphere_radius: 0.1 * MM,
        ..Default::default()
    };

    let sp1 = SimulationParametersArgs {
        offset: -15.0 * MM,
        radius: 10.0 * MM,
        hole_radius: 0.0 * MM,
        phase: 0.0,
        inclination: 17.0,
        ..Default::default()
    };
    let sp2 = SimulationParametersArgs {
        offset: 15.0 * MM,
        radius: 10.0 * MM,
        hole_radius: 0.0 * MM,
        phase: 0.0,
        inclination: -17.0,
        ..Default::default()
    };

    let pressure = compute_pressure_fields(spg, vec![sp1, sp2]);
    plot_and_graph(pressure, (spg, true));
}

fn plot_and_graph(pressure: Field, (spg, multiple_simulation): (SimulationParametersArgs, bool)) {
    // Pressure plotting
    plot_field(
        (&pressure, FieldType::Pressure),
        (None, 4.0),
        (spg, multiple_simulation),
        "./out/field/pressure_field.png",
    )
    .expect("Failed to draw pressure field");
    graph_field(
        (&pressure, FieldType::Pressure),
        (0.0, Some(0.005)),
        (spg, multiple_simulation),
        "./out/graph/pressure_graph.png",
    )
    .expect("Failed to graph pressure field");

    // Potential plotting
    let relative_potential_field = compute_relative_potential_field(&pressure, spg);
    plot_field(
        (&relative_potential_field, FieldType::RadiationPotential),
        (None, 30.0),
        (spg, multiple_simulation),
        "./out/field/acoustic_radiation_potential.png",
    )
    .expect("Failed to draw acoustic field");

    // Force plotting
    let radiation_force_field = compute_radiation_force_field(&relative_potential_field, spg);
    plot_field(
        (&radiation_force_field, FieldType::RadiationForce),
        (None, 6.0),
        (spg, multiple_simulation),
        "./out/field/radiation_force.png",
    )
    .expect("Failed to draw force field");
    graph_field(
        (&radiation_force_field, FieldType::RadiationForce),
        (0.0, Some(0.005)),
        (spg, multiple_simulation),
        "./out/graph/force_graph.png",
    )
    .expect("Failed to graph force field");
}
