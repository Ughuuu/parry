use crate::math::*;

/// A nonlinear motion from a starting isometry traveling at constant translational and rotational velocity.
#[derive(Debug, Copy, Clone)]
pub struct NonlinearRigidMotion {
    /// The starting isometry at `t = 0`.
    pub start: Isometry,
    /// The local-space point at which the rotational part of this motion is applied.
    pub local_center: Point,
    /// The translational velocity of this motion.
    pub linvel: Vector,
    /// The angular velocity of this motion.
    #[cfg(feature = "dim2")]
    pub angvel: Real,
    /// The angular velocity of this motion.
    #[cfg(feature = "dim3")]
    pub angvel: Vector,
}

impl NonlinearRigidMotion {
    /// Initialize a motion from a starting isometry and linear and angular velocities.
    #[cfg(feature = "dim2")]
    pub fn new(start: Isometry, local_center: Point, linvel: Vector, angvel: Real) -> Self {
        NonlinearRigidMotion {
            start,
            local_center,
            linvel,
            angvel,
        }
    }

    /// Initialize a motion from a starting isometry and linear and angular velocities.
    #[cfg(feature = "dim3")]
    pub fn new(start: Isometry, local_center: Point, linvel: Vector, angvel: Vector) -> Self {
        NonlinearRigidMotion {
            start,
            local_center,
            linvel,
            angvel,
        }
    }

    /// Create a `NonlinearRigidMotion` that always returns the identity matrix.
    pub fn identity() -> Self {
        Self::constant_position(Isometry::identity())
    }

    /// Create a `NonlinearRigidMotion` that always return `pos`.
    pub fn constant_position(pos: Isometry) -> Self {
        Self {
            start: pos,
            linvel: Vector::zeros(),
            angvel: AngVector::default(),
            local_center: Point::origin(),
        }
    }

    fn set_start(&mut self, new_start: Isometry) {
        // NOTE: we need to adjust the local_center so that the angular
        // velocity is still expressed wrt. the original center.
        self.local_center =
            new_start.inverse_transform_point(&(self.start.transform_point(&self.local_center)));
        self.start = new_start;
    }

    /// Freezes this motion at the time `t`.
    ///
    /// After calling this, any further calls to `self.position_at_time`
    /// will always return `self.position_at_time(t)` (where `t` is the value given
    /// to this method). This sets the linear velocity and angular velocity
    /// of `self` to zero.
    pub fn freeze(&mut self, t: Real) {
        self.start = self.position_at_time(t);
        self.linvel = Vector::zeros();
        self.angvel = AngVector::default();
    }

    /// Appends a constant translation to this rigid-motion.
    #[must_use]
    pub fn append_translation(&self, tra: Vector) -> Self {
        let mut result = self.clone();
        result.set_start(Translation::from(tra) * result.start);
        result
    }

    /// Prepends a constant translation to this rigid-motion.
    #[must_use]
    pub fn prepend_translation(&self, tra: Vector) -> Self {
        let mut result = self.clone();
        result.set_start(result.start * Translation::from(tra));
        result
    }

    /// Appends a constant isometry to this rigid-motion.
    #[must_use]
    pub fn append(&self, iso: Isometry) -> Self {
        let mut result = self.clone();
        result.set_start(iso * result.start);
        result
    }

    /// Prepends a constant translation to this rigid-motion.
    #[must_use]
    pub fn prepend(&self, iso: Isometry) -> Self {
        let mut result = self.clone();
        result.set_start(result.start * iso);
        result
    }

    /// Computes the position at time `t` of a rigid-body following the motion described by `self`.
    pub fn position_at_time(&self, t: Real) -> Isometry {
        let center = self.start.transform_point(&self.local_center);
        let shift = Translation::from(center.into_vector());
        (shift * Isometry::new(self.linvel * t, self.angvel * t)) * (shift.inverse() * self.start)
    }
}
