//! Splines against OpenUSD's own output. `fixtures/splines.usda` and
//! `fixtures/splines.usdc` are one layer OpenUSD 26.11 wrote in both formats:
//! double, float and half values, Bezier and Hermite curves, every
//! extrapolation mode, an inner loop, a dual-valued knot and knot custom data.

use std::path::Path;

use openusd::sdf::AbstractData;
use openusd::usda::{self, TextWriter};
use openusd::usdc::CrateData;

/// Each `.spline` statement of a layer's text, sorted, since the crate reader
/// does not keep property order.
fn spline_statements(text: &str) -> Vec<String> {
    let mut statements = Vec::new();
    let mut lines = text.lines();
    while let Some(line) = lines.next() {
        if line.contains(".spline = {") {
            let mut statement = vec![line];
            statement.extend(lines.by_ref().take_while(|l| l.trim() != "}"));
            statements.push(statement.join("\n"));
        }
    }
    statements.sort();
    statements
}

fn written(data: &dyn AbstractData) -> Vec<String> {
    spline_statements(&TextWriter::write_to_string(data).expect("write usda"))
}

#[test]
fn text_is_written_as_openusd_writes_it() {
    let expected = std::fs::read_to_string("fixtures/splines.usda").unwrap();
    let data = usda::read_file(Path::new("fixtures/splines.usda")).unwrap();
    assert_eq!(written(&data), spline_statements(&expected));
}

#[test]
fn crate_file_reads_as_openusd_wrote_it() {
    let expected = std::fs::read_to_string("fixtures/splines.usda").unwrap();
    let data = CrateData::open(std::fs::read("fixtures/splines.usdc").unwrap(), true).unwrap();
    assert_eq!(written(&data), spline_statements(&expected));
}

// A uniform attribute declared only by its `.spline`, `.timeSamples` or
// `.connect` statement keeps its variability, which each statement carries
// as C++ `Sdf_WriteAttribute` writes them.
#[test]
fn uniform_is_kept_on_every_statement() {
    let text = "#usda 1.0\n\ndef \"S\"\n{\n    uniform double s.spline = {\n        1: 0; post held,\n    }\n    uniform double t.timeSamples = {\n        1: 2,\n    }\n    uniform token u.connect = </S.c>\n}\n";
    let data = usda::parse(text).unwrap();
    let written = TextWriter::write_to_string(&data).unwrap();
    for statement in [
        "uniform double s.spline",
        "uniform double t.timeSamples",
        "uniform token u.connect",
    ] {
        assert!(written.contains(statement), "{statement} in:\n{written}");
    }
}

fn round_trip(text: &str) -> String {
    TextWriter::write_to_string(&usda::parse(text).unwrap()).unwrap()
}

fn spline_layer(body: &str) -> String {
    format!("#usda 1.0\n\ndef \"S\"\n{{\n{body}}}\n\n")
}

// A knot of b01-027's golden, whose value and tangent sit halfway between two
// shortest spellings; OpenUSD writes the even one.
#[test]
fn tie_digits_are_written_as_openusd_writes_them() {
    let knot = "        96: 3.7029800415039062; pre (0.33333587646484375, -0.022218534709665606); post curve (0.33333587646484375, -0.022218534709665606),\n";
    let text = spline_layer(&format!("    double x.spline = {{\n        bezier,\n{knot}    }}\n"));
    assert!(round_trip(&text).contains(knot));

    let knot = "        1: 0.043945312; pre (0, 0); post held,\n";
    let text = spline_layer(&format!("    float x.spline = {{\n{knot}    }}\n"));
    assert!(round_trip(&text).contains(knot));
}

// OpenUSD keeps no inner loop of all zeros, and ignores a negative loop count.
#[test]
fn loops_are_kept_as_openusd_keeps_them() {
    let text =
        spline_layer("    double x.spline = {\n        loop: (0, 0, -1, 0, 0),\n        1: 0; post held,\n    }\n");
    assert!(!round_trip(&text).contains("loop:"));
    let text =
        spline_layer("    double x.spline = {\n        loop: (1, 2, -1, 3, 0),\n        1: 0; post held,\n    }\n");
    assert!(round_trip(&text).contains("loop: (1, 2, 0, 3, 0),"));
}

// OpenUSD narrows a half knot value through float, so a value just above 1
// becomes 1, where narrowing the double directly gives 1.00098.
#[test]
fn half_values_narrow_through_float() {
    let text = spline_layer("    half x.spline = {\n        1: 1.0004882812509095; pre (0, 0); post held,\n    }\n");
    assert!(round_trip(&text).contains("        1: 1; pre (0, 0); post held,\n"));
}

// A namespaced key in knot custom data is quoted; OpenUSD refuses it bare.
#[test]
fn namespaced_custom_data_keys_are_quoted() {
    let text =
        spline_layer("    double x.spline = {\n        1: 0; pre (0, 0); post held; { int \"a:b\" = 3 },\n    }\n");
    assert!(round_trip(&text).contains(r#"1: 0; pre (0, 0); post held; { int "a:b" = 3 },"#));
}
