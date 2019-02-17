use std::io;
use std::cmp::Ordering;
use std::path::Path;

use gdal::raster::{Dataset, Driver, RasterBand};
use gdal::raster::dataset::GeoTransform;
use gdal::raster::types::GdalType;

use clap::{App, Arg};

#[derive(PartialEq, PartialOrd, Clone, Copy)]
struct ValidF64 {
    v: f64,
}

impl Ord for ValidF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.v.is_nan() {
            Ordering::Equal
        } else {
            self.partial_cmp(other).unwrap_or(Ordering::Equal)
        }
    }
}

impl Eq for ValidF64 {}

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

    fn corners(&self) -> Vec<(ValidF64, ValidF64)> {
        let pt0 = apply_geotransform(self.gt, 0, 0);
        let pt1 = apply_geotransform(self.gt, self.nx, 0);
        let pt2 = apply_geotransform(self.gt, self.nx, self.ny);
        let pt3 = apply_geotransform(self.gt, 0, self.ny);
        vec![pt0, pt1, pt2, pt3].iter().map(|pt| (ValidF64{v: pt.0}, ValidF64{v: pt.1})).collect()
    }

    fn left_extreme(&self) -> ValidF64 {
        let pts = self.corners();
        let extreme = pts.iter().min_by_key(|a| a.0).unwrap();
        extreme.0
    }

    fn right_extreme(&self) -> ValidF64 {
        let pts = self.corners();
        let extreme = pts.iter().max_by_key(|a| a.0).unwrap();
        extreme.0
    }

    fn bottom_extreme(&self) -> ValidF64 {
        let pts = self.corners();
        let extreme = pts.iter().min_by_key(|a| a.1).unwrap();
        extreme.1
    }

    fn top_extreme(&self) -> ValidF64 {
        let pts = self.corners();
        let extreme = pts.iter().max_by_key(|a| a.1).unwrap();
        extreme.1
    }
}

type Point = (ValidF64, ValidF64);

// Pixel size needed to include a point given a geo_transform
fn imsize(gt: GeoTransform, pt: &Point) -> (isize, isize) {
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
        // TODO: ceil incorrect when negative
        (nx.ceil() as isize, ny.ceil() as isize)
    } else {
        let ny = (e*x - e*a - b*y + b*d) / (e*c - b*f);
        let nx = (y - d - ny*f) / e;
        (nx.ceil() as isize, ny.ceil() as isize)
    }
}

// Return the rectangular swatch representing the intersection of a sequence of
// RasterBands. The orientation will be according to the first RasterBand.
fn intersection(bands: &[&RasterBand]) -> Result<Swath, GeoError> {

    if bands.len() == 0 {
        return Err(GeoError{msg: "No bands provided".to_string()})
    }

    let swaths: Vec<Swath> = bands.iter().map(|b| Swath::from_band(b)).collect();

    let left: Vec<ValidF64> = swaths.iter().map(|b| b.left_extreme()).collect();
    let right: Vec<ValidF64> = swaths.iter().map(|b| b.right_extreme()).collect();
    let bottom: Vec<ValidF64> = swaths.iter().map(|b| b.bottom_extreme()).collect();
    let top: Vec<ValidF64> = swaths.iter().map(|b| b.top_extreme()).collect();

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

fn ply_bands<T: Copy + GdalType>(filename: &str,
             swath: Swath,
             projection: String,
             band1: &RasterBand,
             band2: &RasterBand,
             band3: &RasterBand) -> Result<Dataset, io::Error> {
    let driver = Driver::get("JPEG").unwrap();

    let ds = driver.create_with_band_type::<T>(filename,
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


fn main() {

    let cli = App::new("plyband")
                .version("0.1.0")
                .about("Combines satellite bands into false-colour images")
                .author("Nat Wilson")
                .arg(Arg::with_name("red")
                     .short("r")
                     .long("red")
                     .value_name("GEOTIFF")
                     .help("Source of red channel")
                     .required(true))
                .arg(Arg::with_name("green")
                     .short("g")
                     .long("green")
                     .value_name("GEOTIFF")
                     .help("Source of green channel")
                     .required(true))
                .arg(Arg::with_name("blue")
                     .short("b")
                     .long("blue")
                     .value_name("GEOTIFF")
                     .help("Source of blue channel")
                     .required(true))
                .arg(Arg::with_name("output")
                     .short("o")
                     .long("output")
                     .value_name("GEOTIFF")
                     .help("Output file")
                     .required(false))
                .get_matches();

    let red_input = cli.value_of("red").unwrap();
    let green_input = cli.value_of("green").unwrap();
    let blue_input = cli.value_of("blue").unwrap();
    let output = cli.value_of("output").unwrap_or("out.tif");

    let red_ds = Dataset::open(Path::new(red_input)).expect("failed to open red dataset!");
    let green_ds = Dataset::open(Path::new(green_input)).expect("failed to open green dataset!");
    let blue_ds = Dataset::open(Path::new(blue_input)).expect("failed to open blue dataset!");

    let red_band = red_ds.rasterband(1).unwrap();
    let green_band = green_ds.rasterband(1).unwrap();
    let blue_band = blue_ds.rasterband(1).unwrap();

    let proj = red_ds.projection();

    let swath = intersection(&[&red_band, &green_band, &blue_band]).unwrap();

    ply_bands::<u16>(output, swath, proj, &red_band, &green_band, &blue_band).unwrap();

    // TODO:
    // - get right output type
    // - windowed read/write
    // - inputs in different projections
}
