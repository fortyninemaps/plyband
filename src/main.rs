mod error;
mod swath;
mod transform;
mod types;
mod validation;

use std::path::Path;

use gdal::raster::dataset::Buffer;
use gdal::raster::types::GdalType;
use gdal::raster::{Dataset, Driver, RasterBand};

use gdal_sys::GDALDataType;

use gdal_typed_rasterband::typed_rasterband::{GdalFrom, TypedRasterBand};

use clap::{App, Arg};

use error::Error;
use swath::Swath;

struct OutputOptions {
    filename: String,
    format: String,
}

// Convert a string like `"parent/file.tif:2"` into `("parent/file.tif", 2)`
fn split_path_band<'a>(input: &'a str) -> (&'a str, isize) {
    let mut split = input.split(':');
    let path = split.nth(0).unwrap();
    let band = split
        .nth(0)
        .map(|s| str::parse::<isize>(s).unwrap())
        .unwrap_or(1);

    (path, band)
}

fn ply_bands<T: Copy + GdalType + GdalFrom<f64>>(
    output: &OutputOptions,
    sw: Swath,
    projection: String,
    band1: &RasterBand,
    band2: &RasterBand,
    band3: &RasterBand,
) -> Result<Dataset, Error> {
    let driver = Driver::get(&output.format).unwrap();

    let ds = driver
        .create_with_band_type::<T>(&output.filename, sw.nx.abs(), sw.ny.abs(), 3)
        .expect("failed to create output dataset");

    ds.set_projection(&projection)?;
    ds.set_geo_transform(&sw.gt)?;

    let buf: Buffer<T> = TypedRasterBand::from_rasterband(band1)
        .map_err(|e| Error::from(e))
        .and_then(|b| swath::extract(&sw, &b).map_err(|e| Error::from(e)))?;
    ds.write_raster(
        1,
        (0, 0),
        (sw.nx.abs() as usize, sw.ny.abs() as usize),
        &buf,
    )?;

    let buf: Buffer<T> = TypedRasterBand::from_rasterband(band2)
        .map_err(|e| Error::from(e))
        .and_then(|b| swath::extract(&sw, &b).map_err(|e| Error::from(e)))?;
    ds.write_raster(
        2,
        (0, 0),
        (sw.nx.abs() as usize, sw.ny.abs() as usize),
        &buf,
    )?;

    let buf: Buffer<T> = TypedRasterBand::from_rasterband(band3)
        .map_err(|e| Error::from(e))
        .and_then(|b| swath::extract(&sw, &b).map_err(|e| Error::from(e)))?;
    ds.write_raster(
        3,
        (0, 0),
        (sw.nx.abs() as usize, sw.ny.abs() as usize),
        &buf,
    )?;

    Ok(ds)
}

fn main() {
    let cli = App::new("plyband")
        .version("0.1.0")
        .about("Combines satellite bands into false-colour images")
        .author("Nat Wilson")
        .arg(
            Arg::with_name("red")
                .short("r")
                .long("red")
                .value_name("INPUT")
                .help("Source of red channel, with optional band split by a colon (:)")
                .required(true)
                .display_order(1),
        )
        .arg(
            Arg::with_name("green")
                .short("g")
                .long("green")
                .value_name("INPUT")
                .help("Source of green channel, with optional band split by a colon (:)")
                .required(true)
                .display_order(2),
        )
        .arg(
            Arg::with_name("blue")
                .short("b")
                .long("blue")
                .value_name("INPUT")
                .help("Source of blue channel, with optional band split by a colon (:)")
                .required(true)
                .display_order(3),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("INPUT")
                .help("Output file")
                .required(false),
        )
        .arg(
            Arg::with_name("output_format")
                .long("output-format")
                .value_name("FORMAT")
                .help("Output format driver"),
        )
        .get_matches();

    let (red_input, rb) = split_path_band(cli.value_of("red").unwrap());
    let (green_input, gb) = split_path_band(cli.value_of("green").unwrap());
    let (blue_input, bb) = split_path_band(cli.value_of("blue").unwrap());

    let output = cli.value_of("output").unwrap_or("out");
    let output_format = cli.value_of("output_format").unwrap_or("GTIFF");

    let red_ds = Dataset::open(Path::new(red_input)).expect("failed to open red dataset!");
    let green_ds = Dataset::open(Path::new(green_input)).expect("failed to open green dataset!");
    let blue_ds = Dataset::open(Path::new(blue_input)).expect("failed to open blue dataset!");

    let datasets = &[&red_ds, &green_ds, &blue_ds];

    // Input validation
    let checks = validation::have_same_projection(datasets)
        .and_then(|_| validation::have_compatible_geotransforms(datasets));

    match checks {
        Ok(()) => (),
        Err(msg) => {
            eprintln!("Validation failed:\n{}", msg);
            std::process::exit(1);
        }
    };

    let red_band = red_ds.rasterband(rb).unwrap();
    let green_band = green_ds.rasterband(gb).unwrap();
    let blue_band = blue_ds.rasterband(bb).unwrap();

    let proj = red_ds.projection();

    let sw = swath::intersection(&[&red_band, &green_band, &blue_band])
        .expect("Failed to compute intersection between bands");

    let output_options = OutputOptions {
        filename: output.to_string(),
        format: output_format.to_string(),
    };

    let pixel_type = red_band.band_type();

    let result = match pixel_type {
        GDALDataType::GDT_Byte => ply_bands::<u8>(
            &output_options,
            sw,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_UInt16 => ply_bands::<u16>(
            &output_options,
            sw,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_UInt32 => ply_bands::<u32>(
            &output_options,
            sw,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_Int16 => ply_bands::<i16>(
            &output_options,
            sw,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_Int32 => ply_bands::<i32>(
            &output_options,
            sw,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_Float32 => ply_bands::<f32>(
            &output_options,
            sw,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_Float64 => ply_bands::<f64>(
            &output_options,
            sw,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        _ => Err(Error::from_string(
            format!("Unhandled band type {}", pixel_type).to_string(),
        )),
    };

    std::process::exit(match result {
        Ok(_) => 0,
        Err(error) => {
            eprintln!("{}", error);
            1
        }
    });
}
