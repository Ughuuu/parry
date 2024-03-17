//! Bounding sphere.

use crate::bounding_volume::BoundingVolume;
use crate::math::*;
use num::Zero;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// A Bounding Sphere.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Pod, bytemuck::Zeroable))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[derive(Debug, PartialEq, Copy, Clone)]
#[repr(C)]
pub struct BoundingSphere {
    pub center: Point,
    pub radius: Real,
}

impl BoundingSphere {
    /// Creates a new bounding sphere.
    pub fn new(center: Point, radius: Real) -> BoundingSphere {
        BoundingSphere { center, radius }
    }

    /// The bounding sphere center.
    #[inline]
    pub fn center(&self) -> &Point {
        &self.center
    }

    /// The bounding sphere radius.
    #[inline]
    pub fn radius(&self) -> Real {
        self.radius
    }

    /// Transforms this bounding sphere by `m`.
    #[inline]
    pub fn transform_by(&self, m: &Isometry) -> BoundingSphere {
        BoundingSphere::new(m.transform_point(&self.center), self.radius)
    }
}

impl BoundingVolume for BoundingSphere {
    #[inline]
    fn center(&self) -> Point {
        *self.center()
    }

    #[inline]
    fn intersects(&self, other: &BoundingSphere) -> bool {
        // FIXME: refactor that with the code from narrow_phase::ball_ball::collide(...) ?
        let delta_pos = other.center - self.center;
        let distance_squared = delta_pos.norm_squared();
        let sum_radius = self.radius + other.radius;

        distance_squared <= sum_radius * sum_radius
    }

    #[inline]
    fn contains(&self, other: &BoundingSphere) -> bool {
        let delta_pos = other.center - self.center;
        let distance = delta_pos.norm();

        distance + other.radius <= self.radius
    }

    #[inline]
    fn merge(&mut self, other: &BoundingSphere) {
        let mut dir = *other.center() - *self.center();
        let norm = dir.normalize_mut();

        if norm.is_zero() {
            if other.radius > self.radius {
                self.radius = other.radius
            }
        } else {
            let s_center_dir = self.center.as_vector().dot(dir);
            let o_center_dir = other.center.as_vector().dot(dir);

            let right;
            let left;

            if s_center_dir + self.radius > o_center_dir + other.radius {
                right = self.center + dir * self.radius;
            } else {
                right = other.center + dir * other.radius;
            }

            if -s_center_dir + self.radius > -o_center_dir + other.radius {
                left = self.center - dir * self.radius;
            } else {
                left = other.center - dir * other.radius;
            }

            self.center = center(left, right);
            self.radius = distance(right, self.center);
        }
    }

    #[inline]
    fn merged(&self, other: &BoundingSphere) -> BoundingSphere {
        let mut res = self.clone();

        res.merge(other);

        res
    }

    #[inline]
    fn loosen(&mut self, amount: Real) {
        assert!(amount >= 0.0, "The loosening margin must be positive.");
        self.radius = self.radius + amount
    }

    #[inline]
    fn loosened(&self, amount: Real) -> BoundingSphere {
        assert!(amount >= 0.0, "The loosening margin must be positive.");
        BoundingSphere::new(self.center, self.radius + amount)
    }

    #[inline]
    fn tighten(&mut self, amount: Real) {
        assert!(amount >= 0.0, "The tightening margin must be positive.");
        assert!(amount <= self.radius, "The tightening margin is to large.");
        self.radius = self.radius - amount
    }

    #[inline]
    fn tightened(&self, amount: Real) -> BoundingSphere {
        assert!(amount >= 0.0, "The tightening margin must be positive.");
        assert!(amount <= self.radius, "The tightening margin is to large.");
        BoundingSphere::new(self.center, self.radius - amount)
    }
}
