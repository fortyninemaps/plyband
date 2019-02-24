use std::cmp::Ordering;
use std::path::Path;

use gdal::raster::dataset::{Buffer, GeoTransform};
use gdal::raster::types::GdalType;
use gdal::raster::{Dataset, Driver, RasterBand};

use gdal_sys::GDALDataType;

use gdal_typed_rasterband::typed_rasterband::{GdalFrom, TypeError, TypedRasterBand};

use clap::{App, Arg};

#[derive(PartialEq, PartialOrd, Clone, Copy)]
struct RealF64 {
    v: f64,
}

impl Ord for RealF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.v.is_nan() {
            Ordering::Equal
        } else {
            self.partial_cmp(other).unwrap_or(Ordering::Equal)
        }
    }
}

impl Eq for RealF64 {}

type Point = (RealF64, RealF64);

struct OutputOptions {
    filename: String,
    format: String,
}

#[derive(Debug)]
struct Error {
    msg: String,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.msg)
    }
}

impl std::error::Error for Error {
    fn description(&self) -> &str {
        &self.msg
    }

    fn cause(&self) -> Option<&std::error::Error> {
        None
    }
}

impl std::convert::From<TypeError> for Error {
    fn from(_error: TypeError) -> Error {
        Error {
            msg: "band TypeError".to_string(),
        }
    }
}

impl std::convert::From<gdal::errors::Error> for Error {
    fn from(error: gdal::errors::Error) -> Error {
        Error {
            msg: format!("GdalError {}", error),
        }
    }
}

trait Transform2 {
    fn apply(&self, pos: (isize, isize)) -> (f64, f64);
    fn invert(&self, pt: &Point) -> (f64, f64);

    // Pixel size needed to include a point given a geo_transform
    fn imsize(&self, pt: &Point) -> (isize, isize) {
        let (nx, ny) = self.invert(pt);
        // TODO: ceil incorrect when negative
        (nx.ceil() as isize, ny.ceil() as isize)
    }
}

impl Transform2 for GeoTransform {
    fn apply(&self, pos: (isize, isize)) -> (f64, f64) {
        let nx = pos.0 as f64;
        let ny = pos.1 as f64;
        let x = self[0] + nx * self[1] + ny * self[2];
        let y = self[3] + nx * self[4] + ny * self[5];
        (x, y)
    }

    fn invert(&self, pt: &Point) -> (f64, f64) {
        let a = self[0];
        let b = self[1];
        let c = self[2];
        let d = self[3];
        let e = self[4];
        let f = self[5];

        let x = pt.0.v;
        let y = pt.1.v;

        if b != 0.0 {
            let ny = (b * y - b * d - e * x + e * a) / (b * f - e * c);
            let nx = (x - a - ny * c) / b;
            (nx, ny)
        } else {
            let ny = (e * x - e * a - b * y + b * d) / (e * c - b * f);
            let nx = (y - d - ny * f) / e;
            (nx, ny)
        }
    }
}

#[derive(Debug)]
struct Swath {
    nx: isize,
    ny: isize,
    gt: GeoTransform,
    proj: String,
}

impl Swath {
    fn from_band(band: &RasterBand) -> Swath {
        let size = band.owning_dataset().size();
        let gt = band.owning_dataset().geo_transform().unwrap();
        let proj = band.owning_dataset().projection();

        Swath {
            nx: size.0 as isize,
            ny: size.1 as isize,
            gt,
            proj,
        }
    }

    fn corners(&self) -> Vec<(RealF64, RealF64)> {
        let pt0 = self.gt.apply((0, 0));
        let pt1 = self.gt.apply((self.nx, 0));
        let pt2 = self.gt.apply((self.nx, self.ny));
        let pt3 = self.gt.apply((0, self.ny));
        vec![pt0, pt1, pt2, pt3]
            .iter()
            .map(|pt| (RealF64 { v: pt.0 }, RealF64 { v: pt.1 }))
            .collect()
    }

    fn left_extreme(&self) -> RealF64 {
        let pts = self.corners();
        let extreme = pts.iter().min_by_key(|a| a.0).unwrap();
        extreme.0
    }

    fn right_extreme(&self) -> RealF64 {
        let pts = self.corners();
        let extreme = pts.iter().max_by_key(|a| a.0).unwrap();
        extreme.0
    }

    fn bottom_extreme(&self) -> RealF64 {
        let pts = self.corners();
        let extreme = pts.iter().min_by_key(|a| a.1).unwrap();
        extreme.1
    }

    fn top_extreme(&self) -> RealF64 {
        let pts = self.corners();
        let extreme = pts.iter().max_by_key(|a| a.1).unwrap();
        extreme.1
    }
}

// Return the rectangular swatch representing the intersection of a sequence of
// RasterBands. The orientation will be according to the first RasterBand.
fn intersection(bands: &[&RasterBand]) -> Result<Swath, Error> {
    if bands.len() == 0 {
        return Err(Error {
            msg: "No bands provided".to_string(),
        });
    }

    let swaths: Vec<Swath> = bands.iter().map(|b| Swath::from_band(b)).collect();

    let left: Vec<RealF64> = swaths.iter().map(|b| b.left_extreme()).collect();
    let right: Vec<RealF64> = swaths.iter().map(|b| b.right_extreme()).collect();
    let bottom: Vec<RealF64> = swaths.iter().map(|b| b.bottom_extreme()).collect();
    let top: Vec<RealF64> = swaths.iter().map(|b| b.top_extreme()).collect();

    let rightmost_left = left.iter().max().unwrap();
    let leftmost_right = right.iter().max().unwrap();
    let upper_bottom = bottom.iter().max().unwrap();
    let lower_top = top.iter().max().unwrap();

    if (rightmost_left > leftmost_right) || (upper_bottom > lower_top) {
        Err(Error {
            msg: "No valid intersection between bands".to_string(),
        })
    } else {
        let gt_fst = bands[0].owning_dataset().geo_transform().unwrap();
        let proj_fst = bands[0].owning_dataset().projection();

        let gt: [f64; 6] = [
            rightmost_left.v,
            gt_fst[1],
            gt_fst[2],
            lower_top.v,
            gt_fst[4],
            gt_fst[5],
        ];

        let (nx, ny) = gt.imsize(&(*leftmost_right, *upper_bottom));

        Ok(Swath {
            nx,
            ny,
            gt,
            proj: proj_fst,
        })
    }
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
    swath: Swath,
    projection: String,
    band1: &RasterBand,
    band2: &RasterBand,
    band3: &RasterBand,
) -> Result<Dataset, Error> {
    let driver = Driver::get(&output.format).unwrap();

    let ds = driver
        .create_with_band_type::<T>(&output.filename, swath.nx.abs(), swath.ny.abs(), 3)
        .expect("failed to create output dataset");

    ds.set_projection(&projection)?;
    ds.set_geo_transform(&swath.gt)?;

    let buf: Buffer<T> = TypedRasterBand::from_rasterband(band1)
        .map_err(|e| Error::from(e))
        .and_then(|b| b.read_band().map_err(|e| Error::from(e)))?;
    ds.write_raster(
        1,
        (0, 0),
        (swath.nx.abs() as usize, swath.ny.abs() as usize),
        &buf,
    )?;

    let buf: Buffer<T> = TypedRasterBand::from_rasterband(band2)
        .map_err(|e| Error::from(e))
        .and_then(|b| b.read_band().map_err(|e| Error::from(e)))?;
    ds.write_raster(
        2,
        (0, 0),
        (swath.nx.abs() as usize, swath.ny.abs() as usize),
        &buf,
    )?;

    let buf: Buffer<T> = TypedRasterBand::from_rasterband(band3)
        .map_err(|e| Error::from(e))
        .and_then(|b| b.read_band().map_err(|e| Error::from(e)))?;
    ds.write_raster(
        3,
        (0, 0),
        (swath.nx.abs() as usize, swath.ny.abs() as usize),
        &buf,
    )?;

    Ok(ds)
}

fn all_same<I, T>(things: I) -> bool
where
    I: Iterator<Item = T>,
    T: Eq,
{
    let (result, _) = things.fold((true, None), |(b, left_opt), right| {
        if b {
            (
                left_opt.is_none() || left_opt.unwrap() == right,
                Some(right),
            )
        } else {
            (b, Some(right))
        }
    });

    result
}

fn have_same_projection(datasets: &[&Dataset]) -> Result<(), String> {
    let projs = datasets.iter().map(|ds| ds.projection());
    if all_same(projs) {
        Ok(())
    } else {
        Err("Different projections".to_string())
    }
}

fn gt_compatible(gt1: GeoTransform, gt2: GeoTransform) -> Result<(), String> {
    let spacing_check =
        if gt1[1] == gt2[1] && gt1[2] == gt2[2] && gt1[4] == gt2[4] && gt1[5] == gt2[5] {
            Ok(())
        } else {
            return Err("non-equal spacing".to_string());
        };

    // check for integer solutions to initial offsets
    let offset_check = {
        let origin2 = (RealF64 { v: gt2[0] }, RealF64 { v: gt2[3] });
        let (nx, ny) = gt1.invert(&origin2);

        if (nx % 1.0 == 0.0) && (ny % 1.0 == 0.0) {
            Ok(())
        } else {
            return Err("grids offset".to_string());
        }
    };

    spacing_check.and_then(|_| offset_check)
}

fn have_compatible_geotransforms(datasets: &[&Dataset]) -> Result<(), String> {
    if datasets.len() == 0 {
        return Err("empty dataset list".to_string());
    }

    let (result, _) = datasets.iter().fold(
        (Ok(()), datasets[0].geo_transform().unwrap()),
        |(acc, left), right_ds| {
            if acc.is_ok() {
                let right = right_ds.geo_transform().unwrap();
                (gt_compatible(left, right), right)
            } else {
                (acc, left)
            }
        },
    );

    result
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
    have_same_projection(datasets)
        .and_then(|_| have_compatible_geotransforms(datasets))
        .unwrap();

    let red_band = red_ds.rasterband(rb).unwrap();
    let green_band = green_ds.rasterband(gb).unwrap();
    let blue_band = blue_ds.rasterband(bb).unwrap();

    let proj = red_ds.projection();

    let swath = intersection(&[&red_band, &green_band, &blue_band]).unwrap();

    let output_options = OutputOptions {
        filename: output.to_string(),
        format: output_format.to_string(),
    };

    let pixel_type = red_band.band_type();

    let result = match pixel_type {
        GDALDataType::GDT_Byte => ply_bands::<u8>(
            &output_options,
            swath,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_UInt16 => ply_bands::<u16>(
            &output_options,
            swath,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_UInt32 => ply_bands::<u32>(
            &output_options,
            swath,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_Int16 => ply_bands::<i16>(
            &output_options,
            swath,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_Int32 => ply_bands::<i32>(
            &output_options,
            swath,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_Float32 => ply_bands::<f32>(
            &output_options,
            swath,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        GDALDataType::GDT_Float64 => ply_bands::<f64>(
            &output_options,
            swath,
            proj,
            &red_band,
            &green_band,
            &blue_band,
        ),
        _ => Err(Error {
            msg: format!("Unhandled band type {}", pixel_type).to_string(),
        }),
    };

    result.unwrap();

    // TODO:
    // - windowed read/write
    // - inputs in different projections
}
