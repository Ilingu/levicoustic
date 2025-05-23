use std::fs;

use matrix_method::{
    potential_field::{
        compute_limit_density, compute_relative_potential_field, compute_sound_intensity,
    },
    pressure_field::{compute_pressure_field, compute_pressure_fields},
    radiation_force::compute_radiation_force_field,
    Field, FieldType, SimulationParametersArgs, MM,
};
use plot::{
    field::plot_field,
    graph::{graph_field, CutOption},
};

mod matrix_method;
mod plot;

fn main() {
    let _ = fs::create_dir_all("./out/field");
    let _ = fs::create_dir_all("./out/graph");
    let _ = fs::create_dir_all("./tmp");

    let sp = SimulationParametersArgs {
        x_min: -30.0 * MM,
        x_max: 30.0 * MM,
        z_min: 0.0 * MM,
        z_max: 49.0 * MM,

        nb_reflection: 0,
        disc: 0.3 * MM,

        freq: 28_000.0,
        u_0: 0.0000060,

        sphere_radius: 1.0 * MM,
        hole_radius: 0.0 * MM,
        inclination: 0.0,
        offset: 0.0,
        phase: 0.0,
        radius: 28.5 * MM,
    };

    let pressure = compute_pressure_field(sp);
    plot_and_graph(pressure, (sp, false));
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
        (0.0, Some(vec![(0.005, CutOption::Both)])),
        (spg, multiple_simulation),
        "./out/graph/pressure_graph.png",
    )
    .expect("Failed to graph pressure field");

    // Potential plotting
    let relative_potential_field = compute_relative_potential_field(&pressure, spg);
    plot_field(
        (&relative_potential_field, FieldType::RadiationPotential),
        (None, 25.0),
        (spg, multiple_simulation),
        "./out/field/acoustic_radiation_potential.png",
    )
    .expect("Failed to draw acoustic field");

    // Force plotting
    let radiation_force_field = compute_radiation_force_field(&relative_potential_field, spg);
    plot_field(
        (&radiation_force_field, FieldType::RadiationForce),
        (None, 8.0),
        (spg, multiple_simulation),
        "./out/field/radiation_force.png",
    )
    .expect("Failed to draw force field");
    graph_field(
        (&radiation_force_field, FieldType::RadiationForce),
        (0.0, Some(vec![(0.005, CutOption::Both)])),
        (spg, multiple_simulation),
        "./out/graph/force_graph.png",
    )
    .expect("Failed to graph force field");

    return;
    // Density limit plotting
    let density_limit = compute_limit_density(&pressure, spg);
    plot_field(
        (&density_limit, FieldType::LimitDensity),
        (Some(0.0..5000.0), 1.0),
        (spg, multiple_simulation),
        "./out/field/density_limit.png",
    )
    .expect("Failed to draw density field");
    graph_field(
        (&density_limit, FieldType::LimitDensity),
        (
            0.01,
            Some(vec![(0.003, CutOption::Left), (0.0325, CutOption::Right)]),
        ),
        (spg, multiple_simulation),
        "./out/graph/density_graph.png",
    )
    .expect("Failed to graph density");
}
