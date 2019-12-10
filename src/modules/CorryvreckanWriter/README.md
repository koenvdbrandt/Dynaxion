# CorryvreckanWriter
**Maintainer**: Simon Spannagel (<simon.spannagel@cern.ch>), Daniel Hynds (<daniel.hynds@cern.ch>)  
**Status**: Functional  
**Input**: PixelHit

### Description
Takes all digitised pixel hits and converts them into Corryvreckan pixel format. These are then written to an output file in the expected format to be read in by the reconstruction software. Will optionally write out the MC Truth information, storing the MC particle class from Corryvreckan. It is noted that the time resolution is hard-coded as `5ns` for all detectors due to time structure of written out events: events of length `5ns`, with a gap of `10ns` in between events.

This module writes output compatible with Corryvreckan 0.9.12 and later.

### Parameters
* `file_name` : Output filename (file extension `.root` will be appended if not present). Defaults to `corryvreckanOutput.root`
* `geometry_file` : Name of the output geometry file in the Corryvreckan format. Defaults to `corryvreckanGeometry.conf`
* `reference`: Name of the detector used as reference in the reconstruction.
* `dut`: List of detector names to be treated as device under test in the reconstruction. Defaults to an empty list.
* `output_mctruth` : Flag to write out MCParticle information for each hit. Defaults to true.

### Usage
Typical usage is:

```ini
[CorryvreckanWriter]
file_name = corryvreckan
output_mctruth = true
reference = "telescope_plane0"
```
