use std::cmp::Ordering;

use gdal::raster::dataset::GeoTransform;

#[derive(PartialEq, PartialOrd, Clone, Copy)]
pub struct RealF64 {
    pub v: f64,
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

pub type Point = (RealF64, RealF64);

pub trait Transform2 {
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
