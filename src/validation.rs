use gdal::raster::{Dataset, dataset::GeoTransform};

use crate::transform::Transform2;
use crate::types::RealF64;

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

pub fn have_same_projection(datasets: &[&Dataset]) -> Result<(), String> {
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

pub fn have_compatible_geotransforms(datasets: &[&Dataset]) -> Result<(), String> {
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
