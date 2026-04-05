// GPL-3.0-or-later
// Copyright 2024--present radsex-rs developers

//! Population map: maps individual names to groups (e.g. M/F).

use std::collections::HashMap;
use std::io::{self, BufRead};
use std::path::Path;

/// Population map storing individual-to-group assignments and group counts.
#[derive(Debug, Clone, Default)]
pub struct Popmap {
    /// individual_name -> group_name
    pub individual_groups: HashMap<String, String>,
    /// group_name -> count of individuals
    pub group_counts: HashMap<String, u32>,
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
                *popmap.group_counts.entry(group.clone()).or_insert(0) += 1;
                popmap.individual_groups.insert(individual, group);
                popmap.n_individuals += 1;
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

        if n_groups > 2
            && (!self.group_counts.contains_key(&config.group1)
                || !self.group_counts.contains_key(&config.group2))
        {
            return Err(format!(
                "Groups specified with --groups (\"{}\", \"{}\") were not found in popmap groups ({})",
                config.group1, config.group2,
                self.print_groups(false)
            ));
        }

        if n_groups == 2 && (config.group1.is_empty() || config.group2.is_empty()) {
            let mut iter = self.group_counts.keys();
            config.group1 = iter.next().unwrap().clone();
            config.group2 = iter.next().unwrap().clone();
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
