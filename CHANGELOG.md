# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- `Atom.implicit_hydrogens` and `Bond.order` fields, so molecular graphs mirror the
  SMILES/InChI they are read from (stereochemistry, formal charge, and isotopes aside).
- `graph.smiles` to convert a molecular graph to a SMILES string. The output is not
  canonicalized: it is written in graph-key order so it tracks the graph it came from.
- `graph.total_valences`.

### Changed

- `from_smiles` / `from_inchi` now read hydrogens as implicit by default; only standalone
  `[H]` atoms become nodes. Bonds are Kekulized to integer orders on input.
- Renamed `open_valences` to `unpaired_electrons`, redefined as
  bonding capacity minus total valence.

## [0.0.0] - YYYY-MM-DD

### Fixed

- Fix 1
- Fix 2...

### Changed

- Change 1
- Change 2...
