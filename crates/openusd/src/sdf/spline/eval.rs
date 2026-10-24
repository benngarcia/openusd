//! Spline evaluation, as C++ `Ts_Eval` does it (`pxr/base/ts/eval.cpp`).
//!
//! An evaluation first maps its time out of any loop echo into the authored
//! knots, then interpolates or extrapolates there in `f64` with C++'s
//! formulas, operation for operation, and finally rounds the result to the
//! spline's value type, as C++ `TsSpline::Eval` returns it in a `VtValue`.

use std::f64::consts::PI;

use crate::gf;

use super::{CurveType, Extrapolation, Interpolation, Knot, LoopParams, Spline};

impl Spline {
    /// The value at `time` (C++ `TsSpline::Eval`), or `None` where the spline
    /// has none: it has no knots, or `time` falls in a value-block segment or
    /// extrapolation. At a knot this is the value leaving it.
    pub fn eval(&self, time: f64) -> Option<f64> {
        self.evaluate(time, Aspect::Value, Side::Post)
    }

    /// The value approaching `time` from before (C++
    /// `TsSpline::EvalPreValue`). It differs from [`eval`](Self::eval) only at
    /// a discontinuity: a dual-valued knot, the end of a held segment, a reset
    /// loop's seam.
    pub fn eval_pre_value(&self, time: f64) -> Option<f64> {
        self.evaluate(time, Aspect::Value, Side::Pre)
    }

    /// The slope at `time` (C++ `TsSpline::EvalDerivative`), leaving it at a
    /// knot.
    pub fn eval_derivative(&self, time: f64) -> Option<f64> {
        self.evaluate(time, Aspect::Derivative, Side::Post)
    }

    fn evaluate(&self, time: f64, aspect: Aspect, side: Side) -> Option<f64> {
        if self.knots.is_empty() {
            return None;
        }
        let loops = LoopResolver::new(self, time, aspect, side);
        let result = loops.eval()?;
        let sign = if loops.negate { -1.0 } else { 1.0 };
        Some(self.value_type.round((result + loops.value_offset) * sign))
    }

    /// The inner loop and the index of the knot at its prototype start, if
    /// the loop is valid (C++ `Ts_SplineData::HasInnerLoops`).
    fn inner_loop(&self) -> Option<(LoopParams, usize)> {
        let lp = self.inner_loops?;
        if lp.proto_end <= lp.proto_start || (lp.num_pre_loops == 0 && lp.num_post_loops == 0) {
            return None;
        }
        let at = self.knots.partition_point(|k| k.time < lp.proto_start);
        (at != self.knots.len() && self.knots[at].time == lp.proto_start).then_some((lp, at))
    }
}

/// What an evaluation computes (C++ `Ts_EvalAspect`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Aspect {
    Value,
    Derivative,
}

/// Which side of the evaluation time a discontinuity is taken from (C++
/// `Ts_EvalLocation`; `Post` stands for both `Ts_EvalAtTime` and
/// `Ts_EvalPost`, which C++ evaluates alike).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Pre,
    Post,
}

impl Side {
    fn flipped(self) -> Self {
        match self {
            Self::Pre => Self::Post,
            Self::Post => Self::Pre,
        }
    }
}

/// One knot's parameters in `f64` (C++ `Ts_TypedKnotData<double>`). The
/// default is C++'s default-constructed knot, which stands in for a missing
/// neighbour at either end of the spline.
#[derive(Debug, Clone, Copy, Default)]
struct KnotData {
    time: f64,
    pre_width: f64,
    post_width: f64,
    next_interpolation: Interpolation,
    value: f64,
    /// The pre-value of a dual-valued knot.
    dual_value: Option<f64>,
    pre_slope: f64,
    post_slope: f64,
}

/// Where an evaluation in a loop echo evaluates instead, and what it adds to
/// the result there (C++ `_LoopResolver`).
///
/// Extrapolating loops are resolved first, then inner loops: the reverse of
/// the order in which loops copy knots. The flags and edge knots record the
/// interpolation cases inner loops make, where a segment runs to or from a
/// copy of the first prototype knot.
struct LoopResolver<'a> {
    spline: &'a Spline,
    aspect: Aspect,
    /// The inner loop and the index of its first prototype knot.
    inner: Option<(LoopParams, usize)>,
    /// The time to evaluate at, outside every echo.
    time: f64,
    side: Side,
    /// Added to the value found at [`time`](Self::time).
    value_offset: f64,
    /// Whether to negate the result, for a derivative in a reflected
    /// iteration of an oscillating loop.
    negate: bool,
    /// Whether the first or last knot of the spline is an inner-loop echo.
    first_time_looped: bool,
    last_time_looped: bool,
    between_last_proto_and_end: bool,
    between_pre_unlooped_and_looped: bool,
    between_looped_and_post_unlooped: bool,
    /// The two edge knots an extrapolation beyond an echoed edge runs from.
    extrap_knots: (KnotData, KnotData),
}

impl<'a> LoopResolver<'a> {
    fn new(spline: &'a Spline, time: f64, aspect: Aspect, side: Side) -> Self {
        let inner = spline.inner_loop();
        let mut resolver = Self {
            spline,
            aspect,
            inner,
            time,
            side,
            value_offset: 0.0,
            negate: false,
            first_time_looped: false,
            last_time_looped: false,
            between_last_proto_and_end: false,
            between_pre_unlooped_and_looped: false,
            between_looped_and_post_unlooped: false,
            extrap_knots: Default::default(),
        };

        // Valid inner looping always makes at least two knots.
        let multiple_knots = inner.is_some() || spline.knots.len() > 1;
        let pre_loops = multiple_knots && spline.pre_extrapolation.is_looping();
        let post_loops = multiple_knots && spline.post_extrapolation.is_looping();
        if inner.is_none() && !pre_loops && !post_loops {
            return resolver;
        }

        // The first and last knot times, authored or echoed.
        let mut first = spline.knots[0].time;
        let mut last = spline.knots[spline.knots.len() - 1].time;
        if let Some((lp, _)) = inner {
            let looped = lp.looped_interval();
            if looped.min() < first {
                first = looped.min();
                resolver.first_time_looped = true;
            }
            if looped.max() > last {
                last = looped.max();
                resolver.last_time_looped = true;
            }
        }

        if pre_loops || post_loops {
            resolver.resolve_extrapolation(first, last, pre_loops, post_loops);
        }
        if let Some((lp, first_proto)) = inner {
            resolver.resolve_inner(lp, first_proto, first, last);
        }
        resolver
    }

    /// Hops an evaluation in a looping extrapolation into the knots (C++
    /// `_ResolveExtrap` and `_DoExtrap`).
    fn resolve_extrapolation(&mut self, first: f64, last: f64, pre_loops: bool, post_loops: bool) {
        // One end is closed, the other open; which depends on the side.
        let knots = gf::Interval::new(first, last, self.side != Side::Pre, self.side == Side::Pre);
        if knots.contains(self.time) {
            return;
        }
        let (extrapolation, offset, is_pre) = if pre_loops && self.time < last {
            (self.spline.pre_extrapolation, first - self.time, true)
        } else if post_loops && self.time > first {
            (self.spline.post_extrapolation, self.time - last, false)
        } else {
            return;
        };

        let span = last - first;
        // A loop of no length extrapolates as held.
        if span == 0.0 {
            (self.time, self.side) = if is_pre { (first, Side::Pre) } else { (last, Side::Post) };
            return;
        }

        // Hop one more than the whole iterations the offset covers, or exactly
        // that many when the time is on an iteration boundary and evaluated
        // from the side nearer the knots. The integer arithmetic wraps where
        // C++'s overflows, so a time beyond `i32::MAX` iterations is
        // meaningless in both.
        let iterations_frac = offset / span;
        let iterations_trunc = iterations_frac as i32;
        let boundary = f64::from(iterations_trunc) == iterations_frac;
        let short = boundary && (is_pre == (self.side != Side::Pre));
        let iterations = if short {
            iterations_trunc
        } else {
            iterations_trunc.wrapping_add(1)
        };
        let hop = if is_pre { iterations } else { iterations.wrapping_neg() };
        self.time += f64::from(hop) * span;

        match extrapolation {
            // Each iteration starts where the last ended; the value offset
            // runs opposite to the hop.
            Extrapolation::LoopRepeat if self.aspect != Aspect::Derivative => {
                self.value_offset -= f64::from(hop) * self.extrapolation_value_offset(is_pre);
            }
            // Every other iteration is reflected in time.
            Extrapolation::LoopOscillate if hop % 2 != 0 => {
                self.time = first + (span - (self.time - first));
                self.side = self.side.flipped();
                if self.aspect == Aspect::Derivative {
                    self.negate = true;
                }
            }
            // A reset loop's seam falls out of the boundary hop alone.
            _ => {}
        }
    }

    /// The value change over one repeat iteration (C++
    /// `_ComputeExtrapValueOffset`): from the earliest knot's value to the
    /// latest's, either of which may be an inner-loop echo.
    fn extrapolation_value_offset(&self, is_pre: bool) -> f64 {
        let knots = &self.spline.knots;
        let side_value = |k: &Knot| {
            if is_pre {
                k.pre_value.unwrap_or(k.value)
            } else {
                k.value
            }
        };
        let first_value = match self.inner {
            Some((lp, first_proto)) if self.first_time_looped => {
                side_value(&knots[first_proto]) - f64::from(lp.num_pre_loops) * lp.value_offset
            }
            _ => side_value(&knots[0]),
        };
        let last_value = match self.inner {
            // The final echo of the prototype start knot, offset in both its
            // values.
            Some((lp, first_proto)) if self.last_time_looped => {
                side_value(&knots[first_proto]) + f64::from(lp.num_post_loops.wrapping_add(1)) * lp.value_offset
            }
            _ => side_value(&knots[knots.len() - 1]),
        };
        last_value - first_value
    }

    /// Hops an evaluation in an inner-loop echo into the prototype, and notes
    /// the special interpolation cases around the looped region (C++
    /// `_ResolveInner`).
    fn resolve_inner(&mut self, lp: LoopParams, first_proto: usize, first: f64, last: f64) {
        let looped = lp.looped_interval();
        let proto = lp.prototype_interval();

        if looped.contains(self.time) && !proto.contains(self.time) {
            let span = proto.size();
            if self.time < lp.proto_start {
                let iteration = ((lp.proto_start - self.time) / span).ceil() as i32;
                self.time += f64::from(iteration) * span;
                if self.aspect == Aspect::Value {
                    self.value_offset -= f64::from(iteration) * lp.value_offset;
                }
            } else {
                let iteration = (((self.time - lp.proto_end) / span) as i32).wrapping_add(1);
                self.time -= f64::from(iteration) * span;
                if self.aspect == Aspect::Value {
                    self.value_offset += f64::from(iteration) * lp.value_offset;
                }
            }
        }

        let knots = &self.spline.knots;
        // The index of the first knot at or after the prototype end; the one
        // before it is the last prototype knot, which may be the first.
        let proto_end_index = first_proto + knots[first_proto..].partition_point(|k| k.time < lp.proto_end);

        if proto.contains(self.time) {
            // Between the last prototype knot and the prototype end.
            if self.time > knots[proto_end_index - 1].time {
                self.between_last_proto_and_end = true;
            }
        } else if self.time < first {
            if self.first_time_looped {
                // The first knot is a copy of the first prototype knot, and
                // the second a copy of the second, or of the first again when
                // the prototype has only one.
                let pre = lp.num_pre_loops.wrapping_neg();
                let second = match knots.get(first_proto + 1) {
                    Some(k) if proto.contains(k.time) => self.copy_proto_knot(first_proto + 1, pre),
                    _ => self.copy_proto_knot(first_proto, pre.wrapping_add(1)),
                };
                self.extrap_knots = (self.copy_proto_knot(first_proto, pre), second);
            }
        } else if self.time > last {
            if self.last_time_looped {
                // The last knot is a copy of the first prototype knot, and
                // the one before it a copy of the last.
                self.extrap_knots = (
                    self.copy_proto_knot(first_proto, lp.num_post_loops.wrapping_add(1)),
                    self.copy_proto_knot(proto_end_index - 1, lp.num_post_loops),
                );
            }
        } else if self.time < looped.min() {
            // After the last authored knot before the looped region, if any.
            let at = knots[..first_proto].partition_point(|k| k.time < looped.min());
            if at != 0 && self.time > knots[at - 1].time {
                self.between_pre_unlooped_and_looped = true;
            }
        } else if self.time > looped.max() {
            // Before the first authored knot after the looped region, if any.
            let at = first_proto + 1 + knots[first_proto + 1..].partition_point(|k| k.time <= looped.max());
            if at != knots.len() && self.time < knots[at].time {
                self.between_looped_and_post_unlooped = true;
            }
        }
    }

    /// A prototype knot shifted by whole iterations (C++
    /// `_CopyProtoKnotData`).
    fn copy_proto_knot(&self, index: usize, iterations: i32) -> KnotData {
        let (lp, _) = self.inner.expect("only inner loops copy prototype knots");
        let mut knot = KnotData::from(&self.spline.knots[index]);
        knot.time += f64::from(iterations) * lp.prototype_interval().size();
        if self.aspect == Aspect::Value {
            let offset = f64::from(iterations) * lp.value_offset;
            knot.value += offset;
            if let Some(dual) = &mut knot.dual_value {
                *dual += offset;
            }
        }
        knot
    }

    /// Evaluates at the resolved time, before the value offset (C++
    /// `_EvalMain`).
    fn eval(&self) -> Option<f64> {
        let spline = self.spline;
        let knots = &spline.knots;
        let count = knots.len();
        let time = self.time;
        let data = |i: usize| knots.get(i).map(KnotData::from).unwrap_or_default();

        // The bracketing knots, as indices that are `count` where C++'s
        // iterators are `end()`.
        let lower = knots.partition_point(|k| k.time < time);
        let at_knot = lower != count && knots[lower].time == time;
        let prev = if lower != 0 { lower - 1 } else { count };
        let next = if at_knot { lower + 1 } else { lower };
        let multiple_knots = count > 1;
        let (knot_data, mut prev_data, mut next_data) = (data(lower), data(prev), data(next));

        if at_knot {
            return self.eval_at_knot(lower, knot_data, prev_data, next_data);
        }

        if next == 0 {
            if spline.pre_extrapolation == Extrapolation::ValueBlock {
                return None;
            }
            let mut next_data2 = data(1);
            if self.first_time_looped {
                (next_data, next_data2) = self.extrap_knots;
            }
            let slope = extrapolation_slope(
                spline.pre_extrapolation,
                multiple_knots,
                next_data,
                next_data2,
                Side::Pre,
            )?;
            return Some(match self.aspect {
                Aspect::Derivative => slope,
                Aspect::Value => extrapolate_linear(next_data, slope, time, Side::Pre),
            });
        }

        if !self.between_last_proto_and_end && prev == count - 1 {
            if spline.post_extrapolation == Extrapolation::ValueBlock {
                return None;
            }
            let mut prev_data2 = if prev != 0 { data(prev - 1) } else { KnotData::default() };
            if self.last_time_looped {
                (prev_data, prev_data2) = self.extrap_knots;
            }
            let slope = extrapolation_slope(
                spline.post_extrapolation,
                multiple_knots,
                prev_data,
                prev_data2,
                Side::Post,
            )?;
            return Some(match self.aspect {
                Aspect::Derivative => slope,
                Aspect::Value => extrapolate_linear(prev_data, slope, time, Side::Post),
            });
        }

        // Between knots, one of which may be a copy of the first prototype
        // knot (C++ `ReplaceBoundaryKnots`).
        if let Some((lp, first_proto)) = self.inner {
            if self.between_last_proto_and_end {
                next_data = self.copy_proto_knot(first_proto, 1);
            } else if self.between_pre_unlooped_and_looped {
                next_data = self.copy_proto_knot(first_proto, lp.num_pre_loops.wrapping_neg());
            } else if self.between_looped_and_post_unlooped {
                prev_data = self.copy_proto_knot(first_proto, lp.num_post_loops.wrapping_add(1));
            }
        }
        self.interpolate(prev_data, next_data)
    }

    /// Evaluates exactly at the knot at `index`.
    fn eval_at_knot(&self, index: usize, knot: KnotData, prev: KnotData, next: KnotData) -> Option<f64> {
        let spline = self.spline;
        let at_first = index == 0;
        let at_last = index == spline.knots.len() - 1;
        let multiple_knots = spline.knots.len() > 1;
        match (self.aspect, self.side) {
            (Aspect::Value, Side::Pre) => {
                if at_first {
                    if spline.pre_extrapolation == Extrapolation::ValueBlock {
                        return None;
                    }
                } else {
                    match prev.next_interpolation {
                        Interpolation::ValueBlock => return None,
                        // The end of a held segment holds the previous value.
                        Interpolation::Held => return Some(prev.value),
                        Interpolation::Linear | Interpolation::Curve => {}
                    }
                }
                Some(knot.pre_value())
            }
            (Aspect::Value, Side::Post) => {
                let blocked = if at_last {
                    spline.post_extrapolation == Extrapolation::ValueBlock
                } else {
                    knot.next_interpolation == Interpolation::ValueBlock
                };
                (!blocked).then_some(knot.value)
            }
            (Aspect::Derivative, Side::Pre) => {
                if at_first {
                    return extrapolation_slope(spline.pre_extrapolation, multiple_knots, knot, next, Side::Pre);
                }
                match prev.next_interpolation {
                    Interpolation::ValueBlock => None,
                    Interpolation::Held => Some(0.0),
                    Interpolation::Linear => Some(segment_slope(prev, knot)),
                    Interpolation::Curve => Some(knot.pre_slope),
                }
            }
            (Aspect::Derivative, Side::Post) => {
                if at_last {
                    return extrapolation_slope(spline.post_extrapolation, multiple_knots, knot, prev, Side::Post);
                }
                match knot.next_interpolation {
                    Interpolation::ValueBlock => None,
                    Interpolation::Held => Some(0.0),
                    Interpolation::Linear => Some(segment_slope(knot, next)),
                    Interpolation::Curve => Some(knot.post_slope),
                }
            }
        }
    }

    /// Interpolates the segment from `begin` to `end` (C++ `_Interpolate`).
    fn interpolate(&self, begin: KnotData, end: KnotData) -> Option<f64> {
        let time = self.time;
        match begin.next_interpolation {
            Interpolation::ValueBlock => None,
            Interpolation::Curve => Some(match self.spline.curve_type {
                CurveType::Bezier => eval_bezier(begin, end, time, self.aspect),
                CurveType::Hermite => eval_hermite(begin, end, time, self.aspect),
            }),
            Interpolation::Held => Some(match self.aspect {
                Aspect::Value => begin.value,
                Aspect::Derivative => 0.0,
            }),
            Interpolation::Linear => {
                let slope = segment_slope(begin, end);
                Some(match self.aspect {
                    Aspect::Derivative => slope,
                    Aspect::Value => extrapolate_linear(begin, slope, time, Side::Post),
                })
            }
        }
    }
}

impl KnotData {
    /// The value approaching the knot (C++ `GetPreValue`).
    fn pre_value(&self) -> f64 {
        self.dual_value.unwrap_or(self.value)
    }

    /// The pre-tangent's rise from its end to the knot (C++
    /// `GetPreTanHeight`).
    fn pre_height(&self) -> f64 {
        -self.pre_width * self.pre_slope
    }

    /// The post-tangent's rise from the knot to its end (C++
    /// `GetPostTanHeight`).
    fn post_height(&self) -> f64 {
        self.post_width * self.post_slope
    }
}

impl From<&Knot> for KnotData {
    fn from(knot: &Knot) -> Self {
        Self {
            time: knot.time,
            pre_width: knot.pre_tangent.width,
            post_width: knot.post_tangent.width,
            next_interpolation: knot.next_interpolation,
            value: knot.value,
            dual_value: knot.pre_value,
            pre_slope: knot.pre_tangent.slope,
            post_slope: knot.post_tangent.slope,
        }
    }
}

impl LoopParams {
    /// The prototype's span, closed at its start and open at its end, or an
    /// empty interval when the loop makes no copies (C++
    /// `Ts_LoopParams::GetPrototypeInterval`).
    fn prototype_interval(&self) -> gf::Interval {
        if self.proto_end > self.proto_start && (self.num_pre_loops > 0 || self.num_post_loops > 0) {
            gf::Interval::new(self.proto_start, self.proto_end, true, false)
        } else {
            gf::Interval::open(0.0, 0.0)
        }
    }

    /// The prototype's span with its echoes, closed at both ends, or an
    /// empty interval when the loop makes no copies (C++
    /// `Ts_LoopParams::GetLoopedInterval`).
    fn looped_interval(&self) -> gf::Interval {
        if self.proto_end > self.proto_start && (self.num_pre_loops > 0 || self.num_post_loops > 0) {
            let span = self.proto_end - self.proto_start;
            gf::Interval::closed(
                self.proto_start - f64::from(self.num_pre_loops) * span,
                self.proto_end + f64::from(self.num_post_loops) * span,
            )
        } else {
            gf::Interval::open(0.0, 0.0)
        }
    }
}

/// Evaluates a Bezier segment (C++ `_EvalBezier`), first de-regressing a
/// copy of it as evaluation always does.
fn eval_bezier(mut begin: KnotData, mut end: KnotData, time: f64, aspect: Aspect) -> f64 {
    keep_ratio(&mut begin, &mut end);
    let t = bezier_parameter(begin, end, time);
    let value = Cubic::from_points(
        begin.value,
        begin.value + begin.post_height(),
        end.pre_value() + end.pre_height(),
        end.pre_value(),
    );
    match aspect {
        Aspect::Value => value.eval(t),
        // dy/dx as dy/dt over dx/dt.
        Aspect::Derivative => {
            let time = Cubic::from_points(
                begin.time,
                begin.time + begin.post_width,
                end.time - end.pre_width,
                end.time,
            );
            value.derivative().eval(t) / time.derivative().eval(t)
        }
    }
}

/// The parameter in `[0, 1]` at which a Bezier segment's time curve reaches
/// `time` (C++ `_FindBezierParameter`).
fn bezier_parameter(begin: KnotData, end: KnotData, time: f64) -> f64 {
    // Offset by the evaluation time so that the answer is a zero.
    let cubic = Cubic::from_points(
        begin.time - time,
        begin.time + begin.post_width - time,
        end.time - end.pre_width - time,
        end.time - time,
    );
    cubic.monotonic_zero().clamp(0.0, 1.0)
}

/// Evaluates a Hermite segment (C++ `_EvalHermite`) over the slopes alone.
fn eval_hermite(begin: KnotData, end: KnotData, time: f64, aspect: Aspect) -> f64 {
    let (t0, v0, m0) = (begin.time, begin.value, begin.post_slope);
    let (t1, v1, m1) = (end.time, end.pre_value(), end.pre_slope);
    // A segment of no length, or with a NaN end, has its start value.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    if !(t0 < t1) {
        return v0;
    }

    // The cubic over u in [0, 1], with the slopes scaled to match.
    let dt = t1 - t0;
    let dv = v1 - v0;
    let u = (time - t0) / dt;
    let um0 = m0 * dt;
    let um1 = m1 * dt;
    let a = -2.0 * dv + um0 + um1;
    let b = 3.0 * dv - 2.0 * um0 - um1;
    let c = um0;
    let d = v0;
    match aspect {
        Aspect::Derivative => (u * (u * 3.0 * a + 2.0 * b) + c) / dt,
        Aspect::Value => u * (u * (u * a + b) + c) + d,
    }
}

/// The slope of the straight line from one knot to the next (C++
/// `_GetSegmentSlope`).
fn segment_slope(begin: KnotData, end: KnotData) -> f64 {
    (end.pre_value() - begin.value) / (end.time - begin.time)
}

/// The slope beyond an edge knot, from `edge` and the knot `adjacent` to it
/// (C++ `_GetExtrapolationSlope`), or `None` for a value-block extrapolation.
fn extrapolation_slope(
    extrapolation: Extrapolation,
    multiple_knots: bool,
    edge: KnotData,
    adjacent: KnotData,
    side: Side,
) -> Option<f64> {
    match extrapolation {
        Extrapolation::ValueBlock => return None,
        Extrapolation::Held => return Some(0.0),
        Extrapolation::Sloped(slope) => return Some(slope),
        Extrapolation::Linear | Extrapolation::LoopRepeat | Extrapolation::LoopReset | Extrapolation::LoopOscillate => {
        }
    }
    // One knot, a loop already resolved away, or a dual-valued edge knot is
    // flat.
    if !multiple_knots || extrapolation != Extrapolation::Linear || edge.dual_value.is_some() {
        return Some(0.0);
    }
    // Continue the edge segment.
    let segment = match side {
        Side::Pre => edge.next_interpolation,
        Side::Post => adjacent.next_interpolation,
    };
    Some(match (segment, side) {
        (Interpolation::Held | Interpolation::ValueBlock, _) => 0.0,
        (Interpolation::Linear, Side::Pre) => segment_slope(edge, adjacent),
        (Interpolation::Linear, Side::Post) => segment_slope(adjacent, edge),
        (Interpolation::Curve, Side::Pre) => edge.post_slope,
        (Interpolation::Curve, Side::Post) => edge.pre_slope,
    })
}

/// The line through `knot` with `slope`, at `time` on `side` of it (C++
/// `_ExtrapolateLinear`).
fn extrapolate_linear(knot: KnotData, slope: f64, time: f64, side: Side) -> f64 {
    match side {
        Side::Pre => knot.pre_value() - slope * (knot.time - time),
        Side::Post => knot.value + slope * (time - knot.time),
    }
}

/// The amount by which de-regression shortens tangents past the exact
/// limit, so that the result is non-regressive under imprecise processing.
const WRITE_PADDING: f64 = 1e-5;

/// The amount by which a segment must be over-fixed to count as
/// non-regressive; smaller than [`WRITE_PADDING`], so that de-regressed output
/// passes.
const READ_PADDING: f64 = 1e-6;

/// The longest normalized tangent width that is regressive under no
/// circumstances.
const CONTAINED_MAX: f64 = 1.0;

/// Shortens a regressive Bezier segment's tangents so that its time curve
/// never runs backward, keeping their ratio (C++
/// `Ts_RegressionPreventerBatchAccess::ProcessSegment` under
/// `TsAntiRegressionKeepRatio`, the start knot active).
///
/// In widths normalized to the segment, the non-regressive pairs are those
/// inside an ellipse; the line of constant ratio through the original pair
/// meets its edge at the solution with the longer tangents. A width is written
/// back only where the solution changes it.
fn keep_ratio(start: &mut KnotData, end: &mut KnotData) {
    if start.next_interpolation != Interpolation::Curve {
        return;
    }
    let interval = end.time - start.time;
    let interval = if interval > 0.0 { interval } else { 1.0 };
    let active = start.post_width / interval;
    let opposite = end.pre_width / interval;
    if !are_widths_regressive(active, opposite) {
        return;
    }

    let (active_width, opposite_width) = if active < READ_PADDING {
        (active, CONTAINED_MAX - WRITE_PADDING)
    } else if opposite < READ_PADDING {
        (CONTAINED_MAX - WRITE_PADDING, opposite)
    } else {
        let ratio = active / opposite;
        let adjusted_opposite = (ratio.sqrt() + ratio + 1.0) / (ratio * ratio + ratio + 1.0);
        (
            ratio * adjusted_opposite - WRITE_PADDING,
            adjusted_opposite - WRITE_PADDING,
        )
    };
    if active_width != active {
        start.post_width = active_width * interval;
    }
    if opposite_width != opposite {
        end.pre_width = opposite_width * interval;
    }
}

/// Whether normalized tangent widths make a Bezier segment regressive (C++
/// `_AreTanWidthsRegressive`): whether, padded, they lie outside the ellipse.
fn are_widths_regressive(width1: f64, width2: f64) -> bool {
    // Widths inside the unit square never are, even outside the ellipse.
    if width1 <= CONTAINED_MAX && width2 <= CONTAINED_MAX {
        return false;
    }
    let w1 = width1 + READ_PADDING;
    let w2 = width2 + READ_PADDING;
    (w1 * w1) + (w2 * w2) - 2.0 * (w1 + w2) + (w1 * w2) + 1.0 > 0.0
}

/// A quadratic in power form, `at² + bt + c` (C++ `_Quadratic`).
#[derive(Debug, Clone, Copy)]
struct Quadratic {
    a: f64,
    b: f64,
    c: f64,
}

/// One of a Bezier's two cubics, time or value, in power form, `at³ + bt² +
/// ct + d` (C++ `_Cubic`).
#[derive(Debug, Clone, Copy)]
struct Cubic {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
}

impl Quadratic {
    fn eval(&self, t: f64) -> f64 {
        t * (t * self.a + self.b) + self.c
    }

    /// The zero in `[0, 1]` of a quadratic monotonically increasing there,
    /// by the quadratic formula.
    fn monotonic_zero(&self) -> f64 {
        let discriminant = (self.b.powi(2) - 4.0 * self.a * self.c).sqrt();
        let root0 = (-self.b - discriminant) / (2.0 * self.a);
        let root1 = (-self.b + discriminant) / (2.0 * self.a);
        nearest_half(&[root0, root1])
    }
}

impl Cubic {
    /// The cubic through Bezier control points: the start `p0`, the start
    /// tangent's end `p1`, the end tangent's end `p2` and the end `p3`.
    fn from_points(p0: f64, p1: f64, p2: f64, p3: f64) -> Self {
        Self {
            a: -p0 + 3.0 * p1 - 3.0 * p2 + p3,
            b: 3.0 * p0 - 6.0 * p1 + 3.0 * p2,
            c: -3.0 * p0 + 3.0 * p1,
            d: p0,
        }
    }

    fn eval(&self, t: f64) -> f64 {
        t * (t * (t * self.a + self.b) + self.c) + self.d
    }

    fn derivative(&self) -> Quadratic {
        Quadratic {
            a: 3.0 * self.a,
            b: 2.0 * self.b,
            c: self.c,
        }
    }

    /// The zero in `[0, 1]` of a cubic monotonically increasing there,
    /// falling back to the quadratic or linear formula when the leading
    /// coefficients are negligible.
    fn monotonic_zero(&self) -> f64 {
        // C++'s tininess constant, not tuned carefully.
        const EPSILON: f64 = 1e-10;
        let negligible = |x: f64| x.abs() < EPSILON;
        let (a_zero, b_zero, c_zero) = (negligible(self.a), negligible(self.b), negligible(self.c));
        if a_zero && b_zero && c_zero {
            return 0.0;
        }
        if a_zero && b_zero {
            return -self.d / self.c;
        }
        if a_zero {
            return Quadratic {
                a: self.b,
                b: self.c,
                c: self.d,
            }
            .monotonic_zero();
        }
        // Scaled so that the leading coefficient is 1.
        monic_cubic_zero(self.b / self.a, self.c / self.a, self.d / self.a)
    }
}

/// The real zero in `[0, 1]` of the monotonically increasing `t³ + bt² + ct +
/// d`, by Cardano's method (C++ `_FindMonotonicZero(b, c, d)`; see
/// <https://pomax.github.io/bezierinfo/#yforx>).
fn monic_cubic_zero(b: f64, c: f64, d: f64) -> f64 {
    let p = (3.0 * c - b * b) / 3.0;
    let p3 = p / 3.0;
    let p33 = p3 * p3 * p3;
    let q = (2.0 * b * b * b - 9.0 * b * c + 27.0 * d) / 27.0;
    let q2 = q / 2.0;
    let discriminant = q2 * q2 + p33;
    let b3 = b / 3.0;

    if discriminant < 0.0 {
        // Three real roots.
        let r = (-p33).sqrt();
        let t = -q / (2.0 * r);
        let phi = t.clamp(-1.0, 1.0).acos();
        let t1 = 2.0 * r.cbrt();
        let root1 = t1 * (phi / 3.0).cos() - b3;
        let root2 = t1 * ((phi + 2.0 * PI) / 3.0).cos() - b3;
        let root3 = t1 * ((phi + 4.0 * PI) / 3.0).cos() - b3;
        nearest_half(&[root1, root2, root3])
    } else if discriminant == 0.0 {
        // Two real roots.
        let u1 = -q2.cbrt();
        nearest_half(&[2.0 * u1 - b3, -u1 - b3])
    } else {
        // One real root.
        let sd = discriminant.sqrt();
        let u1 = (sd - q2).cbrt();
        let v1 = (sd + q2).cbrt();
        u1 - v1 - b3
    }
}

/// The root nearest 0.5, the first among ties (C++ `_FilterZeros`).
fn nearest_half(roots: &[f64]) -> f64 {
    let mut result = roots[0];
    let mut min_error = (roots[0] - 0.5).abs();
    for &root in &roots[1..] {
        let error = (root - 0.5).abs();
        if error < min_error {
            result = root;
            min_error = error;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use crate::sdf::{SplineValueType, Tangent};

    fn knot(time: f64, value: f64, next_interpolation: Interpolation) -> Knot {
        Knot {
            time,
            value,
            pre_value: None,
            next_interpolation,
            pre_tangent: Tangent::default(),
            post_tangent: Tangent::default(),
            custom_data: HashMap::new(),
        }
    }

    fn spline(knots: impl IntoIterator<Item = Knot>) -> Spline {
        let mut spline = Spline::new(SplineValueType::Double);
        for knot in knots {
            spline.set_knot(knot).unwrap();
        }
        spline
    }

    #[test]
    fn segments_and_sides() {
        let mut dual = knot(2.0, 5.0, Interpolation::ValueBlock);
        dual.pre_value = Some(4.0);
        let spline = spline([
            knot(0.0, 2.0, Interpolation::Linear),
            dual,
            knot(4.0, 1.0, Interpolation::Held),
            knot(6.0, 3.0, Interpolation::Held),
        ]);
        assert_eq!(spline.eval(1.0), Some(3.0));
        assert_eq!(spline.eval_derivative(1.0), Some(1.0));
        assert_eq!((spline.eval_pre_value(2.0), spline.eval(2.0)), (Some(4.0), None));
        assert_eq!(spline.eval(3.0), None);
        assert_eq!((spline.eval_pre_value(4.0), spline.eval(4.0)), (None, Some(1.0)));
        assert_eq!((spline.eval_pre_value(6.0), spline.eval(6.0)), (Some(1.0), Some(3.0)));
        assert_eq!(spline.eval(-10.0), Some(2.0));
        assert_eq!(Spline::new(SplineValueType::Double).eval(0.0), None);
    }

    #[test]
    fn linear_extrapolation() {
        let mut spline = spline([
            knot(0.0, 0.0, Interpolation::Linear),
            knot(2.0, 1.0, Interpolation::Held),
        ]);
        spline.pre_extrapolation = Extrapolation::Linear;
        spline.post_extrapolation = Extrapolation::Sloped(-2.0);
        assert_eq!(spline.eval(-4.0), Some(-2.0));
        assert_eq!(spline.eval(3.0), Some(-1.0));
        spline.post_extrapolation = Extrapolation::ValueBlock;
        assert_eq!(
            (spline.eval_pre_value(2.0), spline.eval(2.0), spline.eval(3.0)),
            (Some(1.0), None, None)
        );
    }

    #[test]
    fn extrapolating_loops() {
        let mut spline = spline([
            knot(0.0, 0.0, Interpolation::Linear),
            knot(2.0, 1.0, Interpolation::Held),
        ]);
        spline.post_extrapolation = Extrapolation::LoopRepeat;
        spline.pre_extrapolation = Extrapolation::LoopRepeat;
        assert_eq!(spline.eval(5.0), Some(2.5));
        assert_eq!(spline.eval(-3.0), Some(-1.5));
        assert_eq!(spline.eval_derivative(5.0), Some(0.5));
        spline.post_extrapolation = Extrapolation::LoopReset;
        assert_eq!((spline.eval_pre_value(4.0), spline.eval(4.0)), (Some(1.0), Some(0.0)));
        spline.post_extrapolation = Extrapolation::LoopOscillate;
        assert_eq!(spline.eval(3.0), Some(0.5));
        assert_eq!(spline.eval_derivative(3.0), Some(-0.5));
        assert_eq!(spline.eval(4.5), Some(0.25));
    }

    #[test]
    fn inner_loop_echoes() {
        let mut spline = spline([
            knot(0.0, 0.0, Interpolation::Linear),
            knot(1.0, 2.0, Interpolation::Linear),
            knot(10.0, 7.0, Interpolation::Held),
        ]);
        spline.inner_loops = Some(LoopParams {
            proto_start: 0.0,
            proto_end: 2.0,
            num_pre_loops: 1,
            num_post_loops: 2,
            value_offset: 10.0,
        });
        // The prototype's last segment runs to an echo of its first knot.
        assert_eq!(spline.eval(1.5), Some(6.0));
        assert_eq!(spline.eval(2.5), Some(11.0));
        assert_eq!(spline.eval(4.5), Some(21.0));
        assert_eq!(spline.eval(-1.5), Some(-9.0));
        // After the looped region, from the echo at its end.
        assert_eq!(spline.eval(8.0), Some(18.5));
        assert_eq!(spline.eval(-3.0), Some(-10.0));
    }

    #[test]
    fn regressive_bezier_kept_in_order() {
        let (mut start, mut end) = (KnotData::default(), KnotData::default());
        start.next_interpolation = Interpolation::Curve;
        (start.post_width, end.time, end.pre_width) = (9.0, 4.0, 6.0);
        keep_ratio(&mut start, &mut end);
        assert!(!are_widths_regressive(start.post_width / 4.0, end.pre_width / 4.0));
        assert!((start.post_width / end.pre_width - 1.5).abs() < 1e-3);

        // A zero width on one side clamps the other to just under the segment.
        (start.post_width, end.pre_width) = (6.0, 0.0);
        keep_ratio(&mut start, &mut end);
        assert_eq!((start.post_width, end.pre_width), ((1.0 - WRITE_PADDING) * 4.0, 0.0));
    }

    #[test]
    fn float_result_rounded() {
        let mut spline = spline([
            knot(0.0, 0.0, Interpolation::Linear),
            knot(3.0, 1.0, Interpolation::Held),
        ]);
        spline.value_type = SplineValueType::Float;
        assert_eq!(spline.eval(1.0), Some(f64::from(1.0_f32 / 3.0)));
    }
}
