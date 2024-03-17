use crate::math::*;
use crate::shape::RoundCone;
use crate::transformation::utils;

impl RoundCone {
    /// Outlines this round cone’s shape using polylines.
    pub fn to_outline(&self, nsubdiv: u32, border_nsubdiv: u32) -> (Vec<Point>, Vec<[u32; 2]>) {
        let r = self.inner_shape.radius;
        let br = self.border_radius;
        let he = self.inner_shape.half_height;

        let mut out_vtx = vec![];
        let mut out_idx = vec![];

        // Compute the profile.
        let center_ab = Point::new(-r, -he, 0.0);
        let center_cd = Point::new(0.0, he, 0.0);
        let side_dir = Vector::new(-2.0 * he, r, 0.0).normalize();

        let a = Point::new(-r, -he - br, 0.0);
        let b = Point::new(-r, -he, 0.0) + side_dir * br;
        let c = Point::new(0.0, he, 0.0) + side_dir * br;
        let d = Point::new(0.0, he + br, 0.0);

        out_vtx.push(a);
        utils::push_arc(center_ab, a, b, border_nsubdiv, &mut out_vtx);
        out_vtx.push(b);
        out_vtx.push(c);
        utils::push_arc(center_cd, c, d, border_nsubdiv, &mut out_vtx);
        out_vtx.push(d);

        let circles = [
            0..1,
            border_nsubdiv..border_nsubdiv + 2,
            border_nsubdiv * 2 + 1..border_nsubdiv * 2 + 2,
        ];
        utils::apply_revolution(false, true, &circles, nsubdiv, &mut out_vtx, &mut out_idx);
        (out_vtx, out_idx)
    }
}
