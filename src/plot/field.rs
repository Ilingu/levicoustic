use std::{fs, ops::Range};

use image::{GenericImageView, ImageReader, RgbaImage};
use ndarray::Array1;
use plotters::prelude::*;

use crate::matrix_method::{Field, FieldType, SimulationParametersArgs};

pub fn plot_field(
    field_info: (&Field, FieldType),
    (zoom, saturation): (Option<Range<f64>>, f64),
    simulation_parameters: SimulationParametersArgs,
    save_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    draw_field(
        field_info,
        simulation_parameters,
        (zoom.clone(), saturation),
    )?;
    draw_field_color_map(field_info, (zoom, saturation))?;
    post_processing(save_path)?;
    Ok(())
}

/* HELPERS */

const TMP_FIELD_PATH: &str = "./tmp/pressure_field.png";
const TMP_COLORMAP_PATH: &str = "./tmp/colormap.png";

fn draw_field(
    (field, field_type): (&Field, FieldType),
    SimulationParametersArgs {
        x_min,
        x_max,
        z_min,
        z_max,
        nb_of_reflection,
        disc,
        ..
    }: SimulationParametersArgs,
    (zoom, saturation): (Option<Range<f64>>, f64),
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(TMP_FIELD_PATH, (600, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(50)
        .caption(
            format!("{field_type} - {nb_of_reflection} reflections - discretization = {disc}m"),
            ("sans-serif", 14),
        )
        .build_cartesian_2d(x_min..x_max, z_min..z_max)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_desc("x (m)")
        .y_desc("z (m)")
        .draw()?;

    let plotting_area = chart.plotting_area();
    let plta = plotting_area.get_pixel_range();
    let (pw, ph) = (plta.0.end - plta.0.start, plta.1.end - plta.1.start);
    let (rm, cm) = match *field.shape() {
        [rm, cm] => (rm, cm),
        _ => panic!(),
    };

    // min max of the field zoom, aka where I want to watch the field
    let (minlvl, maxlvl) = match zoom {
        Some(range) => (range.start, range.end),
        None => best_field_zoom(field, saturation),
    };

    // plot
    for row_id in 0..ph {
        for col_id in 0..pw {
            let (x, z) = (
                x_min + (col_id as f64) * (x_max - x_min) / (pw - 1) as f64,
                z_min + (row_id as f64) * (z_max - z_min) / (ph - 1) as f64,
            );

            // if data point doesn't exist, take the nearest one
            let (mut i, mut j) = (
                (row_id as f64 * rm as f64 / ph as f64).round() as usize,
                (col_id as f64 * cm as f64 / pw as f64).round() as usize,
            );
            if i == rm {
                i -= 1;
            }
            if j == cm {
                j -= 1;
            }
            let p = field[[i, j]];

            plotting_area.draw_pixel(
                (x, z),
                &VulcanoHSL::get_color_normalized(p.re, minlvl, maxlvl),
            )?;
        }
    }

    root.present()?;
    Ok(())
}

fn draw_field_color_map(
    (field, field_type): (&Field, FieldType),
    (zoom, saturation): (Option<Range<f64>>, f64),
) -> Result<(), Box<dyn std::error::Error>> {
    // min max of the field zoom, aka where I want to watch the field
    let (minlvl, maxlvl) = match zoom {
        Some(range) => (range.start, range.end),
        None => best_field_zoom(field, saturation),
    };

    let root = BitMapBackend::new(TMP_COLORMAP_PATH, (150, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .right_y_label_area_size(50)
        .build_cartesian_2d(0_f32..1.0, minlvl as f32..maxlvl as f32)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .disable_x_axis()
        .y_desc(field_type.to_unit())
        .y_label_formatter(&|x| {
            if x.abs() <= 1e-5 && x != &0.0 {
                format!("{:.2e}", x)
            } else {
                format!("{}", x)
            }
        })
        .draw()?;

    let rect_points = Array1::linspace(minlvl, maxlvl, 1000);
    chart.draw_series(rect_points.iter().enumerate().map(|(i, &y)| {
        Rectangle::new(
            [
                (0.5, y as f32),
                (
                    1.0,
                    rect_points[if i == rect_points.len() - 1 { i } else { i + 1 }] as f32,
                ),
            ],
            VulcanoHSL
                .get_color_normalized(y as f32, minlvl as f32, maxlvl as f32)
                .filled(),
        )
    }))?;

    root.present()?;
    Ok(())
}

fn post_processing(save_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let pf_img = ImageReader::open(TMP_FIELD_PATH)?.decode()?;
    let cm_img = ImageReader::open(TMP_COLORMAP_PATH)?.decode()?;

    // image processing
    const CROP_BY: u32 = 25;
    assert!(cm_img.width() > CROP_BY);

    // let pf_img = pf_img.fast_blur(0.3);
    let cm_img = cm_img.crop_imm(CROP_BY, 0, cm_img.width() - CROP_BY, cm_img.height());

    let (w1, w2, h) = (pf_img.width(), cm_img.width(), pf_img.height());
    assert_eq!(h, cm_img.height());

    let mut merged_img = RgbaImage::new(w1 + w2, h);

    for x in 0..w1 {
        for y in 0..h {
            let pixel = pf_img.get_pixel(x, y);
            merged_img.put_pixel(x, y, pixel);
        }
    }
    for x in 0..w2 {
        for y in 0..h {
            let pixel = cm_img.get_pixel(x, y);
            merged_img.put_pixel(x + w1, y, pixel);
        }
    }

    merged_img.save(save_path).unwrap();
    fs::remove_file(TMP_FIELD_PATH)?;
    fs::remove_file(TMP_COLORMAP_PATH)?;
    Ok(())
}

/* UTILS */
fn best_field_zoom(field: &Field, saturation: f64) -> (f64, f64) {
    let (min_pressure, max_pressure) = field
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), &p| {
            (min.min(p.re), max.max(p.re))
        });

    if min_pressure.abs() >= max_pressure.abs() {
        (
            -max_pressure.abs() / saturation,
            max_pressure.abs() / saturation,
        )
    } else {
        (
            -min_pressure.abs() / saturation,
            min_pressure.abs() / saturation,
        )
    }
}
