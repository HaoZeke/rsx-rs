// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

#[test]
fn generated_c_header_matches_the_r_package_copy() {
    let core_header = include_str!("../include/rsx.h");
    let r_header = include_str!("../../rsx-r/src/rsx.h");

    assert_eq!(core_header, r_header);
    assert!(core_header.contains("posterior_linked_probability"));
    assert!(core_header.contains("posterior_null_probability"));
}
