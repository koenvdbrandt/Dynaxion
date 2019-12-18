# ScintillatorPropagation
**Maintainer**: Koen van den Brandt (<kbrandt@nikhef.nl>)
**Status**: Work in Progress
**Input**: DepositedCharge  
**Output**: PropagatedCharge

### Description
Simulates the propagation of photo-electrons through the photo-multiplier stage of a scintillator detector. 
The module mimics the charge multiplication of the gain stages of a given photo-multiplier without having to calculate all trajectories of the photo-electrons. It will also propagate the photo-electrons in time and location relative to the photocathode.

### Parameters
* `gain` : Defines the multiplication of the photo-electrons for each gain stage. Defaults to 10.
* `gain_smearing` : Defines the width of the spread of the gain. Defaults to 0.
* `gain_stages` : Defines the number of gain stages. Defaults to 10.
* `transit_time` : Defines the time it takes for the signal of the photo-multiplier to reach the treshold intensity. Defaults to 50ns.
* `transit_time_spread` : Defines the spread of the transit time. Defaults to 10ns.
* `rise_time` : Defines the time it takes for the signal of the photo-multiplier to go from no signal to its maximum value. Defaults to 10ns
* `rise_time_spread` : Defines the spread in the rise time. Defaults to 2ns.
* `pm_prop` : Defines the propagation in space relative to the center of the photocathode. Defaults to 0mm 0mm 0mm.


### Plotting parameters
* `output_plots` : Determines if simple output plots should be generated for a monitoring of the simulation flow. Disabled by default.
* `output_scale_electrons` : Defines the scale of the plot axis to plot the total number of photo-electrons that reach the end of the photo-multiplier stage. Only used if `output_plots` is enabled.  Defaults to (10000 * **gain** * **gain-stages**) photo-electrons.



### Usage
A example of scintillator propagation for detector _detector1_ for a photo-multiplier with a gain of 2 and 10 gain stages. The photo-multiplier has a transit time of 100ns, a rise-time of 10ns and all photo-electrons will be deposited on a position which is 50mm away from the middle of the photocathode in local coordinates:

```ini
[ScintillatorPropagation]
name = "detector1"
gain = 2
gain_stages = 10
transit_time = 100ns
rise_time = 10ns
pm_prop = 0mm 0mm 50mm


```


