use crate::math::*;
use crate::query::Ray;
use simba::simd::SimdValue;

/// A structure representing 4 rays in an SIMD SoA fashion.
#[derive(Debug, Copy, Clone)]
pub struct SimdRay {
    /// The origin of the rays represented as a single SIMD point.
    pub origin: SimdPoint,
    /// The direction of the rays represented as a single SIMD vector.
    pub dir: SimdVector,
}

impl SimdRay {
    /// Creates a new SIMD ray with all its lanes filled with the same ray.
    pub fn splat(ray: Ray) -> Self {
        Self {
            origin: SimdPoint::splat(ray.origin.into()),
            dir: SimdVector::splat(ray.dir.into()),
        }
    }
}
