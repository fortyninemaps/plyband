use std::cmp::Ordering;
use std::path::Path;

use gdal::raster::{Dataset, Driver, RasterBand};
use gdal::raster::dataset::GeoTransform;
use gdal::raster::types::GdalType;

use gdal_sys::GDALDataType;

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

struct OutputOptions {
    filename: String,
    format: String,
}

#[derive(Debug)]
struct GeoError {
    msg: String,
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

        Swath{
            nx: size.0 as isize,
            ny: size.1 as isize,
            gt,
            proj,
        }
    }

    fn corners(&self) -> Vec<(RealF64, RealF64)> {
        let pt0 = apply_geotransform(self.gt, 0, 0);
        let pt1 = apply_geotransform(self.gt, self.nx, 0);
        let pt2 = apply_geotransform(self.gt, self.nx, self.ny);
        let pt3 = apply_geotransform(self.gt, 0, self.ny);
        vec![pt0, pt1, pt2, pt3].iter().map(|pt| (RealF64{v: pt.0}, RealF64{v: pt.1})).collect()
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

type Point = (RealF64, RealF64);

// Pixel size needed to include a point given a geo_transform
fn invert_geotransform(gt: GeoTransform, pt: &Point) -> (f64, f64) {
    let a = gt[0];
    let b = gt[1];
    let c = gt[2];
    let d = gt[3];
    let e = gt[4];
    let f = gt[5];

    let x = pt.0.v;
    let y = pt.1.v;

    if b != 0.0 {
        let ny = (b*y - b*d - e*x + e*a) / (b*f - e*c);
        let nx = (x - a - ny*c) / b;
        (nx, ny)
    } else {
        let ny = (e*x - e*a - b*y + b*d) / (e*c - b*f);
        let nx = (y - d - ny*f) / e;
        (nx, ny)
    }
}


// Pixel size needed to include a point given a geo_transform
fn imsize(gt: GeoTransform, pt: &Point) -> (isize, isize) {
    let (nx, ny) = invert_geotransform(gt, pt);
    // TODO: ceil incorrect when negative
    (nx.ceil() as isize, ny.ceil() as isize)
}

// Return the rectangular swatch representing the intersection of a sequence of
// RasterBands. The orientation will be according to the first RasterBand.
fn intersection(bands: &[&RasterBand]) -> Result<Swath, GeoError> {

    if bands.len() == 0 {
        return Err(GeoError{msg: "No bands provided".to_string()})
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
        Err(
            GeoError{msg: "No valid intersection between bands".to_string()}
        )
    } else {
        let gt_fst = bands[0].owning_dataset().geo_transform().unwrap();
        let proj_fst = bands[0].owning_dataset().projection();

        let gt: [f64; 6] = [rightmost_left.v, gt_fst[1], gt_fst[2],
                            lower_top.v, gt_fst[4], gt_fst[5]];

        let (nx, ny) = imsize(gt, &(*leftmost_right, *upper_bottom));

        Ok(
            Swath{nx, ny, gt, proj: proj_fst}
        )
    }
}

fn apply_geotransform(gt: GeoTransform, xpx: isize, ypx: isize) -> (f64, f64) {
    let x = gt[0] + xpx as f64 * gt[1] + ypx as f64 * gt[2];
    let y = gt[3] + xpx as f64 * gt[4] + ypx as f64 * gt[5];
    (x, y)
}

fn ply_bands<T: Copy + GdalType>(
            output: &OutputOptions,
            swath: Swath,
            projection: String,
            band1: &RasterBand,
            band2: &RasterBand,
            band3: &RasterBand) -> Result<Dataset, GeoError> {
    let driver = Driver::get(&output.format).unwrap();

    let ds = driver.create_with_band_type::<T>(&output.filename,
                                               swath.nx.abs(),
                                               swath.ny.abs(),
                                               3).expect("failed to create output dataset");

    ds.set_projection(&projection).unwrap();
    ds.set_geo_transform(&swath.gt).unwrap();

    let buf = band1.read_band_as::<T>().unwrap();
    ds.write_raster(1, (0, 0), (swath.nx.abs() as usize, swath.ny.abs() as usize), &buf).unwrap();

    let buf = band2.read_band_as::<T>().unwrap();
    ds.write_raster(2, (0, 0), (swath.nx.abs() as usize, swath.ny.abs() as usize), &buf).unwrap();

    let buf = band3.read_band_as::<T>().unwrap();
    ds.write_raster(3, (0, 0), (swath.nx.abs() as usize, swath.ny.abs() as usize), &buf).unwrap();

    Ok(ds)
}

fn all_same<I, T>(things: I) -> bool
where
    I: Iterator<Item = T>,
    T: Eq,
{
    let mut last: Option<T> = None;
    for thing in things {
        let same = match last {
            Some(t) => t == thing,
            None => true
        };
        if !same {
            return false;
        }
        last = Some(thing);
    }
    true
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
    let spacing_check = if gt1[1] == gt2[1] && gt1[2] == gt2[2] && gt1[4] == gt2[4] && gt1[5] == gt2[5] {
        Ok(())
    } else {
        return Err("non-equal spacing".to_string());
    };

    // check for integer solutions to initial offsets
    let offset_check = {
        let origin2 = (RealF64{v: gt2[0]}, RealF64{v: gt2[3]});
        let (nx, ny) = invert_geotransform(gt1, &origin2);

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
    for (i, dataset) in datasets.iter().enumerate() {
        if i != 0 {
            if gt_compatible(datasets[i-1].geo_transform().unwrap(), dataset.geo_transform().unwrap()).is_err() {
                return Err("mismatching geotransforms".to_string());
            }
        }
    }

    Ok(())
}

fn main() {

    let cli = App::new("plyband")
                .version("0.1.0")
                .about("Combines satellite bands into false-colour images")
                .author("Nat Wilson")
                .arg(Arg::with_name("red")
                     .short("r")
                     .long("red")
                     .value_name("INPUT")
                     .help("Source of red channel")
                     .required(true)
                     .display_order(1))
                .arg(Arg::with_name("green")
                     .short("g")
                     .long("green")
                     .value_name("INPUT")
                     .help("Source of green channel")
                     .required(true)
                     .display_order(2))
                .arg(Arg::with_name("blue")
                     .short("b")
                     .long("blue")
                     .value_name("INPUT")
                     .help("Source of blue channel")
                     .required(true)
                     .display_order(3))
                .arg(Arg::with_name("output")
                     .short("o")
                     .long("output")
                     .value_name("INPUT")
                     .help("Output file")
                     .required(false))
                .arg(Arg::with_name("output_format")
                     .long("output-format")
                     .value_name("FORMAT")
                     .help("Output format driver"))
                .get_matches();

    let red_input = cli.value_of("red").unwrap();
    let green_input = cli.value_of("green").unwrap();
    let blue_input = cli.value_of("blue").unwrap();
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

    let red_band = red_ds.rasterband(1).unwrap();
    let green_band = green_ds.rasterband(1).unwrap();
    let blue_band = blue_ds.rasterband(1).unwrap();

    let proj = red_ds.projection();

    let swath = intersection(&[&red_band, &green_band, &blue_band]).unwrap();

    let output_options = OutputOptions{
        filename: output.to_string(),
        format: output_format.to_string(),
    };

    let pixel_type = red_band.band_type();

    let result = match pixel_type {
        GDALDataType::GDT_Byte =>
            ply_bands::<u8>(&output_options, swath, proj, &red_band, &green_band, &blue_band),
        GDALDataType::GDT_UInt16 =>
            ply_bands::<u16>(&output_options, swath, proj, &red_band, &green_band, &blue_band),
        GDALDataType::GDT_UInt32 =>
            ply_bands::<u32>(&output_options, swath, proj, &red_band, &green_band, &blue_band),
        GDALDataType::GDT_Int16 =>
            ply_bands::<i16>(&output_options, swath, proj, &red_band, &green_band, &blue_band),
        GDALDataType::GDT_Int32 =>
            ply_bands::<i32>(&output_options, swath, proj, &red_band, &green_band, &blue_band),
        GDALDataType::GDT_Float32 =>
            ply_bands::<f32>(&output_options, swath, proj, &red_band, &green_band, &blue_band),
        GDALDataType::GDT_Float64 =>
            ply_bands::<f64>(&output_options, swath, proj, &red_band, &green_band, &blue_band),
        _ => Err(GeoError{msg: format!("Unhandled band type {}", pixel_type).to_string()}),
    };

    result.unwrap();


    // TODO:
    // - windowed read/write
    // - inputs in different projections
}
