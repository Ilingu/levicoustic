use std::f64::consts::PI;

use plotters::prelude::*;

use crate::matrix_method::{Field, FieldType, SimulationParametersArgs, C, RHO};

pub fn graph_field(
    (field, field_type): (&Field, FieldType),
    (at_x, cut): (f64, Option<f64>),
    sp: SimulationParametersArgs,
    path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let SimulationParametersArgs {
        x_min,
        x_max,
        z_min,
        z_max,
        freq,
        sphere_radius,
        inclination,
        curvature,
        ..
    } = sp;

    let cut = cut.unwrap_or_default();
    let (cz_min, cz_max) = (z_min + cut, z_max - cut);

    assert!(at_x >= x_min && at_x <= x_max);

    let (rm, cm) = match *field.shape() {
        [rm, cm] => (rm, cm),
        _ => panic!(),
    };
    let x_colid = (((at_x - x_min) * ((cm - 1) as f64)) / (x_max - x_min)).round() as usize;
    let data = field
        .column(x_colid)
        .indexed_iter()
        .filter_map(|(i, fv)| {
            let z = z_min + (i as f64) * (z_max - z_min) / ((rm - 1) as f64);
            if z < cz_min || z > cz_max {
                None
            } else {
                Some((z, fv.re))
            }
        })
        .collect::<Vec<_>>();
    let (dmin, dmax) = data.iter().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(min, max), &(_, fv)| (min.min(fv), max.max(fv)),
    );

    let root = BitMapBackend::new(path, (600, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("{field_type} at x={at_x}"),
            ("sans-serif", 14).into_font(),
        )
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(50)
        .build_cartesian_2d(cz_min..cz_max, dmin..dmax)?;

    chart
        .configure_mesh()
        .x_desc("z (m)")
        .y_desc(field_type.to_unit())
        .y_label_formatter(&|x| {
            if x.abs() <= 1e-5 && x != &0.0 {
                format!("{:.2e}", x)
            } else {
                format!("{}", x)
            }
        })
        .draw()?;

    // verification
    if field_type == FieldType::RadiationForce && inclination == 0.0 && curvature == 0.0 {
        let k = 2.0 * PI * freq.re / C.re;
        chart.draw_series(LineSeries::new(
            data.iter().map(|&(z, _)| {
                (
                    z,
                    (5.0 * PI * sphere_radius.powi(3) * k * 2200_f64.powi(2)
                        / (6.0 * RHO.re * C.re.powi(2)))
                        * (2.0 * k * z).sin(),
                )
            }),
            &GREEN,
        ))?;
    }

    chart.draw_series(LineSeries::new(data, &RED))?;

    root.present()?;
    Ok(())
}
