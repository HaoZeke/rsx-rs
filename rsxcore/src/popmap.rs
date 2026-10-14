// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Population map: maps individual names to groups (e.g. M/F).

use indexmap::IndexMap;
use std::collections::HashMap;
use std::io::{self, BufRead};
use std::path::Path;

/// Population map storing individual-to-group assignments and group counts.
/// Uses IndexMap for group_counts to preserve insertion order (matching C++
/// unordered_map behavior where first-seen group becomes group1).
#[derive(Debug, Clone, Default)]
pub struct Popmap {
    /// individual_name -> group_name
    pub individual_groups: HashMap<String, String>,
    /// group_name -> count of individuals (insertion-ordered)
    pub group_counts: IndexMap<String, u32>,
    /// Total number of individuals
    pub n_individuals: u16,
}

/// Configuration for group comparison.
#[derive(Debug, Clone, Default)]
pub struct GroupConfig {
    pub group1: String,
    pub group2: String,
}

impl Popmap {
    /// Load a popmap from a TSV file (individual\tgroup per line).
    pub fn from_file(path: &Path) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let reader = io::BufReader::new(file);
        let mut popmap = Popmap::default();

        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;
            let line = line.trim_end_matches('\r');
            if line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() == 2 {
                let individual = fields[0].to_string();
                let group = fields[1].to_string();
                if popmap.individual_groups.contains_key(&individual) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "duplicate individual \"{individual}\" in popmap (line {})",
                            line_num + 1
                        ),
                    ));
                }
                *popmap.group_counts.entry(group.clone()).or_insert(0) += 1;
                popmap.individual_groups.insert(individual, group);
                popmap.n_individuals = popmap.n_individuals.checked_add(1).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "popmap exceeds 65535 individuals",
                    )
                })?;
            } else {
                log::warn!(
                    "Could not process line {} in popmap (expected format: \"individual\\tgroup\")",
                    line_num + 1
                );
            }
        }

        Ok(popmap)
    }

    /// Resolve group1/group2 from the popmap when not specified by the user.
    /// Returns an error message if validation fails.
    pub fn resolve_groups(&self, config: &mut GroupConfig) -> Result<(), String> {
        let n_groups = self.group_counts.len();

        if n_groups < 2 {
            return Err(format!(
                "Found {} group(s) in popmap ({}) but at least two are required",
                n_groups,
                self.print_groups(false)
            ));
        }

        if n_groups > 2 && (config.group1.is_empty() || config.group2.is_empty()) {
            return Err(format!(
                "Found {} groups in popmap ({}) but groups to compare were not defined (use --groups group1,group2)",
                n_groups,
                self.print_groups(false)
            ));
        }

        if n_groups == 2 && (config.group1.is_empty() || config.group2.is_empty()) {
            let mut iter = self.group_counts.keys();
            config.group1 = iter
                .next()
                .expect("n_groups==2 guarantees two keys")
                .clone();
            config.group2 = iter
                .next()
                .expect("n_groups==2 guarantees two keys")
                .clone();
        }

        if !config.group1.is_empty()
            && !config.group2.is_empty()
            && (!self.group_counts.contains_key(&config.group1)
                || !self.group_counts.contains_key(&config.group2))
        {
            return Err(format!(
                "Groups specified with --groups (\"{}\", \"{}\") were not found in popmap groups ({})",
                config.group1,
                config.group2,
                self.print_groups(false)
            ));
        }

        if !config.group1.is_empty() && config.group1 == config.group2 {
            return Err(format!(
                "Groups to compare must be different, got \"{}\" twice",
                config.group1
            ));
        }

        Ok(())
    }

    /// Reject a comparison whose marker-table columns outnumber the popmap.
    ///
    /// Duplicate sample names set more than one presence bit for one popmap
    /// row. The resulting presence count cannot be a cell of a table whose
    /// margin is `get_count`, and `distrib` would index off the end of that
    /// margin.
    pub fn check_compared_columns(
        &self,
        group1: &str,
        columns1: u32,
        group2: &str,
        columns2: u32,
    ) -> Result<(), String> {
        for (group, columns) in [(group1, columns1), (group2, columns2)] {
            let total = self.get_count(group);
            if columns > total {
                return Err(format!(
                    "group \"{group}\" is assigned to {columns} marker-table columns but the popmap lists {total}; duplicate sample names make presence counts exceed the group total"
                ));
            }
        }
        Ok(())
    }

    /// Get the group for an individual.
    pub fn get_group(&self, individual: &str) -> Option<&str> {
        self.individual_groups.get(individual).map(|s| s.as_str())
    }

    /// Get the count of individuals in a group.
    pub fn get_count(&self, group: &str) -> u32 {
        self.group_counts.get(group).copied().unwrap_or(0)
    }

    /// Format groups for display.
    pub fn print_groups(&self, with_counts: bool) -> String {
        let mut parts: Vec<String> = self
            .group_counts
            .iter()
            .map(|(g, c)| {
                if with_counts {
                    format!("\"{g}\": {c}")
                } else {
                    format!("\"{g}\"")
                }
            })
            .collect();
        parts.sort(); // deterministic output
        parts.join(", ")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp_popmap(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    #[test]
    fn test_basic_popmap() {
        let f = write_temp_popmap("ind1\tM\nind2\tF\nind3\tM\n");
        let popmap = Popmap::from_file(f.path()).unwrap();
        assert_eq!(popmap.n_individuals, 3);
        assert_eq!(popmap.get_count("M"), 2);
        assert_eq!(popmap.get_count("F"), 1);
        assert_eq!(popmap.get_group("ind1"), Some("M"));
    }

    #[test]
    fn test_windows_line_endings() {
        let f = write_temp_popmap("ind1\tM\r\nind2\tF\r\n");
        let popmap = Popmap::from_file(f.path()).unwrap();
        assert_eq!(popmap.n_individuals, 2);
        assert_eq!(popmap.get_group("ind2"), Some("F"));
    }

    #[test]
    fn test_resolve_groups_two() {
        let f = write_temp_popmap("ind1\tM\nind2\tF\n");
        let popmap = Popmap::from_file(f.path()).unwrap();
        let mut config = GroupConfig::default();
        popmap.resolve_groups(&mut config).unwrap();
        assert!(!config.group1.is_empty());
        assert!(!config.group2.is_empty());
        assert_ne!(config.group1, config.group2);
    }

    #[test]
    fn test_resolve_groups_too_few() {
        let f = write_temp_popmap("ind1\tM\nind2\tM\n");
        let popmap = Popmap::from_file(f.path()).unwrap();
        let mut config = GroupConfig::default();
        assert!(popmap.resolve_groups(&mut config).is_err());
    }
}
