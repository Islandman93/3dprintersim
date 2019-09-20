# 3dprintersim

The purpose of this simulator is to enable testing of new ideas for 3d printer motion, ranging from move planning to acceleration profiles. It uses a simple damped spring model for independent X & Y axis motion to determine the toolhead deviation for a move or set of moves. This simple surrogate model was validated through acceleration data capture discussed in https://github.com/KevinOConnor/klipper/issues/57#issuecomment-479302453 and my own personal test prints. The point is, the simulator has been validated to real world conditions and gives an absolute measure of performance instead of anecdotal judgment of print quality.

## Current support

* Klipper motion planner with square corner velocity algorithm.
* Constant, Smoothstep, & Smootherstep acceleration profiles.
* Dynamic acceleration (like Duet) and jerk limited scurve from [dmbutyugin](https://github.com/dmbutyugin)
* GCode support is limited but should work for testing single layers. 


### Requires
matplotlib, numpy, pymunk

### Usage Examples
python main.py klipper ./examples/test.gcode 50 1000 0.032 0.032 5 -at constant

python main.py klipper ./examples/test.gcode 50 1000 0.032 0.032 5 -at dynamic

python main.py klipper ./examples/test.gcode 50 1000 0.032 0.032 5 -at destructive

