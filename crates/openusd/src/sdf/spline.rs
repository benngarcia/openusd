//! A spline-valued attribute's authored curve (C++ `TsSpline`, spec §7.4.2.4).
//!
//! [`Spline`] holds what a layer stores: the knots, the curve type, the
//! extrapolation on either side and any inner loop. It does not evaluate the
//! curve. The usda reader and writer spell it in the native spline grammar,
//! and the usdc reader and writer encode it as C++ `Ts_BinaryDataAccess` does
//! (`pxr/base/ts/binary.cpp`), version 1: the encoding crate version 0.12.0
//! carries, which has no tangent algorithms, loop boundary times or time-code
//! values.

use std::collections::HashMap;

use strum::FromRepr;

use crate::gf::f16;

use super::Value;

/// The authored data of a spline: its knots in ascending time, and how the
/// curve runs between and beyond them.
#[derive(Debug, Clone, PartialEq)]
pub struct Spline {
    pub value_type: SplineValueType,
    pub curve_type: CurveType,
    pub pre_extrapolation: Extrapolation,
    pub post_extrapolation: Extrapolation,
    /// A repeated stretch of knots inside the curve (C++ inner loop
    /// parameters, deprecated in OpenUSD but still read and written).
    pub inner_loops: Option<LoopParams>,
    knots: Vec<Knot>,
}

/// The type of a spline's values, which is its attribute's type. Each value
/// is held in an `f64`, rounded to this type when it is set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplineValueType {
    Double,
    Float,
    Half,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, FromRepr)]
#[repr(u8)]
pub enum CurveType {
    #[default]
    Bezier = 0,
    /// Tangent widths are fixed at a third of the segment, so only slopes are
    /// authored.
    Hermite = 1,
}

/// How a spline continues before its first knot or after its last.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub enum Extrapolation {
    /// No value.
    ValueBlock,
    /// The edge knot's value.
    #[default]
    Held,
    /// The line through the edge knot along its tangent.
    Linear,
    /// The line through the edge knot with this slope.
    Sloped(f64),
    /// The knots repeated, each copy offset so the ends meet.
    LoopRepeat,
    /// The knots repeated exactly.
    LoopReset,
    /// The knots repeated, every other copy reversed.
    LoopOscillate,
}

/// How a segment runs from a knot to the next.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, FromRepr)]
#[repr(u8)]
pub enum Interpolation {
    /// No value.
    ValueBlock = 0,
    /// The knot's value, up to the next knot.
    #[default]
    Held = 1,
    Linear = 2,
    /// A Bezier or Hermite curve, as the spline's curve type says.
    Curve = 3,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LoopParams {
    pub proto_start: f64,
    pub proto_end: f64,
    pub num_pre_loops: i32,
    pub num_post_loops: i32,
    pub value_offset: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Knot {
    pub time: f64,
    pub value: f64,
    /// The value approached from before the knot, when it differs from
    /// [`value`](Self::value) (a dual-valued knot).
    pub pre_value: Option<f64>,
    /// How the segment after this knot runs.
    pub next_interpolation: Interpolation,
    pub pre_tangent: Tangent,
    pub post_tangent: Tangent,
    pub custom_data: HashMap<String, Value>,
}

/// A knot's tangent on one side: its length in time and its slope.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Tangent {
    pub width: f64,
    pub slope: f64,
}

impl Spline {
    pub fn new(value_type: SplineValueType) -> Self {
        Self {
            value_type,
            curve_type: CurveType::default(),
            pre_extrapolation: Extrapolation::default(),
            post_extrapolation: Extrapolation::default(),
            inner_loops: None,
            knots: Vec::new(),
        }
    }

    /// The knots, in ascending time.
    pub fn knots(&self) -> &[Knot] {
        &self.knots
    }

    /// Adds a knot, replacing any at the same time, with its values rounded
    /// to the spline's value type. Like C++ `TsKnot`, it refuses a value or
    /// slope the type cannot hold finitely and a negative tangent width.
    pub fn set_knot(&mut self, mut knot: Knot) -> Result<(), String> {
        let round = |v: f64| {
            let rounded = self.value_type.round(v);
            match rounded.is_finite() {
                true => Ok(rounded),
                false => Err(format!(
                    "knot at time {}: {v} is not a finite {:?}",
                    knot.time, self.value_type
                )),
            }
        };
        let value = round(knot.value)?;
        let pre_value = knot.pre_value.map(round).transpose()?;
        let pre_slope = round(knot.pre_tangent.slope)?;
        let post_slope = round(knot.post_tangent.slope)?;
        if knot.pre_tangent.width < 0.0 || knot.post_tangent.width < 0.0 {
            return Err(format!("knot at time {}: a tangent width is negative", knot.time));
        }
        (knot.value, knot.pre_value) = (value, pre_value);
        (knot.pre_tangent.slope, knot.post_tangent.slope) = (pre_slope, post_slope);
        match self.knots.binary_search_by(|k| k.time.total_cmp(&knot.time)) {
            Ok(at) => self.knots[at] = knot,
            Err(at) => self.knots.insert(at, knot),
        }
        Ok(())
    }

    /// Whether any segment between two knots is a curve (C++
    /// `TsKnotMap::HasCurveSegments`); the last knot's interpolation leads
    /// nowhere.
    pub fn has_curve_segments(&self) -> bool {
        let segments = self.knots.len().saturating_sub(1);
        self.knots[..segments]
            .iter()
            .any(|k| k.next_interpolation == Interpolation::Curve)
    }

    /// Encodes the spline as C++ `Ts_BinaryDataAccess::GetBinaryData` does,
    /// version 1. Knot custom data is not part of it: the crate format stores
    /// that beside the blob, keyed by knot time.
    pub(crate) fn to_binary(&self) -> Vec<u8> {
        let mut out = Vec::new();
        let hermite = self.curve_type == CurveType::Hermite;
        out.push(BINARY_VERSION | (self.value_type.descriptor() << 4) | ((self.curve_type as u8) << 7));
        out.push(
            self.pre_extrapolation.mode()
                | (self.post_extrapolation.mode() << 3)
                | (u8::from(self.inner_loops.is_some()) << 6),
        );
        for extrapolation in [self.pre_extrapolation, self.post_extrapolation] {
            if let Extrapolation::Sloped(slope) = extrapolation {
                out.extend(slope.to_le_bytes());
            }
        }
        if let Some(lp) = self.inner_loops {
            out.extend(lp.proto_start.to_le_bytes());
            out.extend(lp.proto_end.to_le_bytes());
            out.extend(lp.num_pre_loops.to_le_bytes());
            out.extend(lp.num_post_loops.to_le_bytes());
            out.extend(lp.value_offset.to_le_bytes());
        }
        out.extend((self.knots.len() as u32).to_le_bytes());
        for knot in &self.knots {
            out.push(
                u8::from(knot.pre_value.is_some())
                    | ((knot.next_interpolation as u8) << 1)
                    | ((self.curve_type as u8) << 3),
            );
            out.extend(knot.time.to_le_bytes());
            self.value_type.put(&mut out, knot.value);
            if let Some(pre_value) = knot.pre_value {
                self.value_type.put(&mut out, pre_value);
            }
            if !hermite {
                out.extend(knot.pre_tangent.width.to_le_bytes());
                out.extend(knot.post_tangent.width.to_le_bytes());
            }
            self.value_type.put(&mut out, knot.pre_tangent.slope);
            self.value_type.put(&mut out, knot.post_tangent.slope);
        }
        out
    }

    /// Decodes a blob [`to_binary`](Self::to_binary) writes, as C++
    /// `Ts_BinaryDataAccess::CreateSplineFromBinaryData` does. An empty blob is
    /// a spline with no knots, of the attribute's value type.
    pub(crate) fn from_binary(bytes: &[u8], value_type: SplineValueType) -> Result<Self, String> {
        let mut spline = Spline::new(value_type);
        if bytes.is_empty() {
            return Ok(spline);
        }
        let mut blob = Blob(bytes);
        let header = blob.u8()?;
        let version = header & 0x0f;
        if version != BINARY_VERSION {
            return Err(format!(
                "spline data version {version} (only version {BINARY_VERSION} is supported)"
            ));
        }
        if header & 0x40 != 0 {
            return Err("time-valued spline (not supported)".into());
        }
        spline.value_type = SplineValueType::from_descriptor((header & 0x30) >> 4)?;
        spline.curve_type = CurveType::from_repr((header & 0x80) >> 7).expect("one bit");
        let hermite = spline.curve_type == CurveType::Hermite;

        let modes = blob.u8()?;
        spline.pre_extrapolation = Extrapolation::from_mode(modes & 0x07, &mut blob)?;
        spline.post_extrapolation = Extrapolation::from_mode((modes & 0x38) >> 3, &mut blob)?;
        if modes & 0x40 != 0 {
            spline.inner_loops = Some(LoopParams {
                proto_start: blob.f64()?,
                proto_end: blob.f64()?,
                num_pre_loops: blob.i32()?,
                num_post_loops: blob.i32()?,
                value_offset: blob.f64()?,
            });
        }

        let count = blob.u32()?;
        for _ in 0..count {
            let flags = blob.u8()?;
            let next_interpolation = Interpolation::from_repr((flags & 0x06) >> 1).expect("two bits");
            let time = blob.f64()?;
            let value = blob.value(spline.value_type)?;
            let pre_value = if flags & 0x01 != 0 {
                Some(blob.value(spline.value_type)?)
            } else {
                None
            };
            let (pre_width, post_width) = if hermite {
                (0.0, 0.0)
            } else {
                (blob.f64()?, blob.f64()?)
            };
            let pre_slope = blob.value(spline.value_type)?;
            let post_slope = blob.value(spline.value_type)?;
            spline.knots.push(Knot {
                time,
                value,
                pre_value,
                next_interpolation,
                pre_tangent: Tangent {
                    width: pre_width,
                    slope: pre_slope,
                },
                post_tangent: Tangent {
                    width: post_width,
                    slope: post_slope,
                },
                custom_data: HashMap::new(),
            });
        }
        if !blob.0.is_empty() {
            return Err(format!("{} bytes after the last knot", blob.0.len()));
        }
        Ok(spline)
    }

    /// Attaches each knot's custom data, which the crate format stores apart
    /// from the knots, keyed by knot time.
    pub(crate) fn set_custom_data(&mut self, time: f64, custom_data: HashMap<String, Value>) -> Result<(), String> {
        let knot = self
            .knots
            .iter_mut()
            .find(|k| k.time == time)
            .ok_or_else(|| format!("custom data for time {time}, where the spline has no knot"))?;
        knot.custom_data = custom_data;
        Ok(())
    }
}

impl LoopParams {
    /// C++'s default `TsLoopParams`, which means no inner loop.
    pub const NONE: Self = Self {
        proto_start: 0.0,
        proto_end: 0.0,
        num_pre_loops: 0,
        num_post_loops: 0,
        value_offset: 0.0,
    };
}

/// The only spline encoding crate version 0.12.0 carries.
const BINARY_VERSION: u8 = 1;

impl SplineValueType {
    /// The value type of a spline on an attribute of this type, if a spline
    /// can hold it.
    pub fn from_type_name(name: &str) -> Option<Self> {
        match name {
            "double" => Some(Self::Double),
            "float" => Some(Self::Float),
            "half" => Some(Self::Half),
            _ => None,
        }
    }

    /// The value as this type holds it.
    pub fn round(self, value: f64) -> f64 {
        match self {
            Self::Double => value,
            Self::Float => f64::from(value as f32),
            // Through a float, as C++ narrows a double to `GfHalf`.
            Self::Half => f16::from_f32(value as f32).to_f64(),
        }
    }

    fn descriptor(self) -> u8 {
        match self {
            Self::Double => 1,
            Self::Float => 2,
            Self::Half => 3,
        }
    }

    fn from_descriptor(descriptor: u8) -> Result<Self, String> {
        match descriptor {
            0 | 1 => Ok(Self::Double),
            2 => Ok(Self::Float),
            3 => Ok(Self::Half),
            other => Err(format!("spline value type descriptor {other}")),
        }
    }

    fn put(self, out: &mut Vec<u8>, value: f64) {
        match self {
            Self::Double => out.extend(value.to_le_bytes()),
            Self::Float => out.extend((value as f32).to_le_bytes()),
            Self::Half => out.extend(f16::from_f32(value as f32).to_le_bytes()),
        }
    }
}

impl Extrapolation {
    fn mode(self) -> u8 {
        match self {
            Self::ValueBlock => 0,
            Self::Held => 1,
            Self::Linear => 2,
            Self::Sloped(_) => 3,
            Self::LoopRepeat => 4,
            Self::LoopReset => 5,
            Self::LoopOscillate => 6,
        }
    }

    fn from_mode(mode: u8, blob: &mut Blob<'_>) -> Result<Self, String> {
        Ok(match mode {
            0 => Self::ValueBlock,
            1 => Self::Held,
            2 => Self::Linear,
            3 => Self::Sloped(blob.f64()?),
            4 => Self::LoopRepeat,
            5 => Self::LoopReset,
            6 => Self::LoopOscillate,
            other => return Err(format!("spline extrapolation mode {other}")),
        })
    }

    /// Whether the knots repeat beyond the edge.
    pub fn is_looping(self) -> bool {
        matches!(self, Self::LoopRepeat | Self::LoopReset | Self::LoopOscillate)
    }
}

/// The unread rest of a spline blob, little-endian as C++ writes it on every
/// platform OpenUSD supports.
struct Blob<'a>(&'a [u8]);

impl Blob<'_> {
    fn take<const N: usize>(&mut self) -> Result<[u8; N], String> {
        let (head, rest) = self.0.split_first_chunk::<N>().ok_or("spline data ends early")?;
        self.0 = rest;
        Ok(*head)
    }

    fn u8(&mut self) -> Result<u8, String> {
        Ok(self.take::<1>()?[0])
    }

    fn u32(&mut self) -> Result<u32, String> {
        Ok(u32::from_le_bytes(self.take()?))
    }

    fn i32(&mut self) -> Result<i32, String> {
        Ok(i32::from_le_bytes(self.take()?))
    }

    fn f64(&mut self) -> Result<f64, String> {
        Ok(f64::from_le_bytes(self.take()?))
    }

    fn value(&mut self, value_type: SplineValueType) -> Result<f64, String> {
        Ok(match value_type {
            SplineValueType::Double => self.f64()?,
            SplineValueType::Float => f64::from(f32::from_le_bytes(self.take()?)),
            SplineValueType::Half => f16::from_le_bytes(self.take()?).to_f64(),
        })
    }
}

/// The shape of the AOUSD file-format baselines
/// (`file_formats/tests/assets/text/baseline/splines.json`): held
/// extrapolation and a missing inner loop are `null`, and knot custom data is
/// keyed by knot time apart from the knots.
#[cfg(feature = "serde")]
impl serde::Serialize for Spline {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        #[derive(serde::Serialize)]
        #[serde(rename_all = "camelCase")]
        struct Baseline<'a> {
            curve_type: &'static str,
            pre_extrapolation: Option<ExtrapolationEntry>,
            post_extrapolation: Option<ExtrapolationEntry>,
            loop_parameters: Option<LoopEntry>,
            knots: Vec<KnotEntry>,
            knot_custom_data: HashMap<String, &'a HashMap<String, Value>>,
        }
        #[derive(serde::Serialize)]
        struct ExtrapolationEntry {
            mode: &'static str,
            slope: f64,
        }
        #[derive(serde::Serialize)]
        #[serde(rename_all = "camelCase")]
        struct LoopEntry {
            proto_start: f64,
            proto_end: f64,
            num_pre_loops: i32,
            num_post_loops: i32,
            value_offset: f64,
        }
        #[derive(serde::Serialize)]
        #[serde(rename_all = "camelCase")]
        struct KnotEntry {
            time: f64,
            value: f64,
            pre_value: f64,
            pre_tangent_slope: f64,
            pre_tangent_width: f64,
            post_tangent_slope: f64,
            post_tangent_width: f64,
            next_interpolation_mode: &'static str,
        }

        let extrapolation = |e: Extrapolation| match e {
            Extrapolation::Held => None,
            Extrapolation::Sloped(slope) => Some(ExtrapolationEntry { mode: "sloped", slope }),
            other => Some(ExtrapolationEntry {
                mode: other.keyword(),
                slope: 0.0,
            }),
        };
        Baseline {
            curve_type: self.curve_type.keyword(),
            pre_extrapolation: extrapolation(self.pre_extrapolation),
            post_extrapolation: extrapolation(self.post_extrapolation),
            loop_parameters: self.inner_loops.map(|lp| LoopEntry {
                proto_start: lp.proto_start,
                proto_end: lp.proto_end,
                num_pre_loops: lp.num_pre_loops,
                num_post_loops: lp.num_post_loops,
                value_offset: lp.value_offset,
            }),
            knots: self
                .knots
                .iter()
                .map(|k| KnotEntry {
                    time: k.time,
                    value: k.value,
                    pre_value: k.pre_value.unwrap_or(0.0),
                    pre_tangent_slope: k.pre_tangent.slope,
                    pre_tangent_width: k.pre_tangent.width,
                    post_tangent_slope: k.post_tangent.slope,
                    post_tangent_width: k.post_tangent.width,
                    next_interpolation_mode: k.next_interpolation.keyword(),
                })
                .collect(),
            knot_custom_data: self
                .knots
                .iter()
                .filter(|k| !k.custom_data.is_empty())
                .map(|k| (time_key(k.time), &k.custom_data))
                .collect(),
        }
        .serialize(serializer)
    }
}

#[cfg(feature = "serde")]
fn time_key(time: f64) -> String {
    if time.fract() == 0.0 && time.is_finite() {
        format!("{}", time as i64)
    } else {
        format!("{time}")
    }
}

impl CurveType {
    /// The keyword usda spells it with.
    pub fn keyword(self) -> &'static str {
        match self {
            Self::Bezier => "bezier",
            Self::Hermite => "hermite",
        }
    }
}

impl Extrapolation {
    /// The keyword usda spells its mode with (C++
    /// `Sdf_FileIOUtility::Stringify(TsExtrapMode)`).
    pub fn keyword(self) -> &'static str {
        match self {
            Self::ValueBlock => "none",
            Self::Held => "held",
            Self::Linear => "linear",
            Self::Sloped(_) => "sloped",
            Self::LoopRepeat => "loop repeat",
            Self::LoopReset => "loop reset",
            Self::LoopOscillate => "loop oscillate",
        }
    }
}

impl Interpolation {
    /// The keyword usda spells it with.
    pub fn keyword(self) -> &'static str {
        match self {
            Self::ValueBlock => "none",
            Self::Held => "held",
            Self::Linear => "linear",
            Self::Curve => "curve",
        }
    }
}
