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
