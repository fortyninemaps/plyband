use gdal::raster::dataset::{Buffer, GeoTransform};
use gdal::raster::types::GdalType;
use gdal::raster::RasterBand;

use gdal_typed_rasterband::typed_rasterband::{GdalFrom, TypedRasterBand};

use crate::error::Error;
use crate::transform::Transform2;
use crate::types::RealF64;

#[derive(Debug)]
pub struct Swath {
    pub nx: isize,
    pub ny: isize,
    pub gt: GeoTransform,
    pub proj: String,
}

impl Swath {
    fn from_band(band: &RasterBand) -> Swath {
        let size = band.owning_dataset().size();
        let gt = band
            .owning_dataset()
            .geo_transform()
            .expect("band has no associated geotransform");
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

    pub fn left_extreme(&self) -> RealF64 {
        let pts = self.corners();
        let extreme = pts.iter().min_by_key(|a| a.0).unwrap();
        extreme.0
    }

    pub fn right_extreme(&self) -> RealF64 {
        let pts = self.corners();
        let extreme = pts.iter().max_by_key(|a| a.0).unwrap();
        extreme.0
    }

    pub fn bottom_extreme(&self) -> RealF64 {
        let pts = self.corners();
        let extreme = pts.iter().min_by_key(|a| a.1).unwrap();
        extreme.1
    }

    pub fn top_extreme(&self) -> RealF64 {
        let pts = self.corners();
        let extreme = pts.iter().max_by_key(|a| a.1).unwrap();
        extreme.1
    }
}

pub fn extract<T: Copy + GdalType + GdalFrom<f64>>(
    swath: &Swath,
    band: &TypedRasterBand<T>,
) -> Result<Buffer<T>, Error> {
    let top_left_coord = (RealF64 { v: swath.gt[0] }, RealF64 { v: swath.gt[3] });

    let (ix, iy) = band
        .owning_dataset()
        .geo_transform()?
        .invert(&top_left_coord);
    let top_left_idx = (ix as isize, iy as isize);
    let size = (swath.nx as usize, swath.ny as usize);

    // we require inputs to have the same resolution, so the buffer size will be the same as the
    // window read
    let buf = band.read(top_left_idx, size, size)?;
    Ok(buf)
}

// Return the rectangular swath representing the intersection of a sequence of
// RasterBands. The orientation will be according to the first RasterBand.
pub fn intersection(bands: &[&RasterBand]) -> Result<Swath, Error> {
    if bands.len() == 0 {
        return Err(Error::from_string("No bands provided".to_string()));
    }

    let swaths: Vec<Swath> = bands.iter().map(|b| Swath::from_band(b)).collect();

    let left: Vec<RealF64> = swaths.iter().map(|b| b.left_extreme()).collect();
    let right: Vec<RealF64> = swaths.iter().map(|b| b.right_extreme()).collect();
    let bottom: Vec<RealF64> = swaths.iter().map(|b| b.bottom_extreme()).collect();
    let top: Vec<RealF64> = swaths.iter().map(|b| b.top_extreme()).collect();

    let rightmost_left = left.iter().max().unwrap();
    let leftmost_right = right.iter().min().unwrap();
    let upper_bottom = bottom.iter().max().unwrap();
    let lower_top = top.iter().min().unwrap();

    if (rightmost_left > leftmost_right) || (upper_bottom > lower_top) {
        Err(Error::from_string(
            "No valid intersection between bands".to_string(),
        ))
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
