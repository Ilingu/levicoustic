use std::fs;

use image::{GenericImageView, ImageReader, RgbaImage};
use ndarray::Array2;
use num_complex::Complex;
use plotters::prelude::*;

use crate::matrix_method::{X_MAX, X_MIN, Z_MAX, Z_MIN};

pub fn plot_pressure_field(
    pressure: &Array2<Complex<f64>>,
    path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    draw_pressure_field(pressure)?;
    draw_color_map(pressure)?;
    post_processing(path)?;
    Ok(())
}

/* UTILS */
const SATURATION: f64 = 3.0;
fn pressure_field_levels(pressure: &Array2<Complex<f64>>, saturation: f64) -> (f64, f64) {
    let min_pressure = pressure.iter().fold(f64::INFINITY, |min, &p| min.min(p.re));
    (min_pressure / saturation, -min_pressure / saturation)
}

/* HELPERS */

fn draw_pressure_field(pressure: &Array2<Complex<f64>>) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("./tmp/pressure_field.png", (600, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(50)
        .caption("Pressure field", ("sans-serif", 20))
        .build_cartesian_2d(X_MIN..X_MAX, Z_MIN..Z_MAX)?;

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
    let (rm, cm) = match *pressure.shape() {
        [rm, cm] => (rm, cm),
        _ => panic!(),
    };

    let (minlvl, maxlvl) = pressure_field_levels(pressure, SATURATION);
    for row_id in 0..ph {
        for col_id in 0..pw {
            let (x, z) = (
                X_MIN + (col_id as f64) * (X_MAX - X_MIN) / pw as f64,
                Z_MIN + (row_id as f64) * (Z_MAX - Z_MIN) / ph as f64,
            );
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
            let p = pressure[[i, j]];

            plotting_area.draw_pixel(
                (x, z),
                &VulcanoHSL::get_color_normalized(p.re, minlvl, maxlvl),
            )?;
        }
    }

    root.present()?;
    Ok(())
}

fn draw_color_map(pressure: &Array2<Complex<f64>>) -> Result<(), Box<dyn std::error::Error>> {
    let (minlvl, maxlvl) = pressure_field_levels(pressure, SATURATION);

    let root = BitMapBackend::new("./tmp/colormap.png", (150, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .y_label_area_size(50)
        .build_cartesian_2d(0_f32..1.0, minlvl as f32..maxlvl as f32)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .disable_x_axis()
        .y_desc("Acoustic pressure (Pa)")
        .draw()?;

    chart.draw_series((minlvl as isize..maxlvl as isize).map(|y| {
        Rectangle::new(
            [(0.0, y as f32), (0.5, y as f32 + 1.0)],
            VulcanoHSL
                .get_color_normalized(y as f32, minlvl as f32, maxlvl as f32)
                .filled(),
        )
    }))?;

    root.present()?;
    Ok(())
}

fn post_processing(path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let pf_img = ImageReader::open("./tmp/pressure_field.png")?.decode()?;
    let cm_img = ImageReader::open("./tmp/colormap.png")?.decode()?;

    // image processing
    const CROP_BY: u32 = 25;
    assert!(cm_img.width() > CROP_BY);

    // let pf_img = pf_img.fast_blur(0.3);
    let cm_img = cm_img.crop_imm(0, 0, cm_img.width() - CROP_BY, cm_img.height());

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

    merged_img.save(path).unwrap();
    fs::remove_file("./tmp/pressure_field.png")?;
    fs::remove_file("./tmp/colormap.png")?;
    Ok(())
}
