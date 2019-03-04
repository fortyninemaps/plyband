use std::cmp::Ordering;

#[derive(PartialEq, PartialOrd, Clone, Copy, Debug)]
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
