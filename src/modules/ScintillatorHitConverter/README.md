# ScintillatorHitConverter
**Maintainer**: Koen van den Brandt (<kbrandt@nikhef.nl>)
**Status**: Functional
**Input**: ScintHit
**Output**: DepositedCharge

### Description
Converts the optical photons hitting the sensor of a scintillator to photo-electrons on a photo-cathode.
The modules takes the wavelength of a photon that hit the sensor and will use a provided quantum efficiency file to check if the photon will cause the release of an photo-electron.

### Parameters
* `quantum_efficiency` : Location of file containing the quantum efficiecny of the photocathode.

### Plotting parameters
* `output_plots` : Determines if simple output plots should be generated for a monitoring of the simulation flow. Disabled by default.
* `output_scale_wavelength` : Defines the scale of the plot axis to plot the wavelength of both the incoming optical photons and the optical photons that create a photo-electron. Only used if `output_plots` is enabled. Defaults to 1000nm.
* `output_scale_photo_electrons` : Defines the scale of the plot axis to plot the total number of photo-electrons created per event. Only used if `output_plots` is enabled.  Defaults to 10000 photo-electrons.
### Usage
A example of scintillator hit convter which takes an input file quantum_efficiency.mac.

```toml
[ScintillatorHitConverter]
quantum_efficiency = "quantum_efficiency.mac"
```

