use std::fs;

use matrix_method::pressure_field;
use plot::plot_pressure_field;

mod matrix_method;
mod plot;

fn main() {
    let _ = fs::create_dir_all("./out");
    let _ = fs::create_dir_all("./tmp");

    let pressure = pressure_field::compute_pressure_field();
    plot_pressure_field(&pressure, "./out/pressure_field.png").expect("Failed to draw pf");
}
