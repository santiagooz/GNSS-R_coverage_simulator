# LEO satellite coverage simulation for GNSS-R

Set of Matlab scripts to simulate the coverage of GNSS-R events with GPS signals from an LEO satellite trajectory.

## Description

This set of scripts calculates the position of the specular reflection points in a bistatic geometry with the GPS constellation based on the input LEO trajectory. Then, modelling the antenna beam as an elliptical cone with the input aperture angles, determines which reflections are being captured, i. e. the specular reflection points that are within the modelled antenna footprint.

The results are the location and number of specular reflection points of the signals captured during the simulation time interval.

Three relevant characteristics of the reflections are given as output as well: incidence angle, path-loss for diffuse reflections and, if an antenna pattern file is available, the receiver antenna gain for each reflection. These results are given for a global coverage and for the Argentinean territory.

It also gives statistics in the presence of signals in limited time registers. This is useful to plan the capture of raw IF intermediate signals where the registers are time limited due to the high number of samples to store and downlink.

A series of plots are generated to illutrate the result and histograms to present statistical characteristics of the captured reflections.

## Getting Started

The input parameters for the simulation are defined in the configuration file, including the simulation time interval and paths to the necessary file. The mandatory input files are two: a rinex file with the GNSS ephemeris data and a file containing the LEO satellite trajectory in ECEF coordinates. Additionaly, a file with the antenna gain pattern can be loaded to extract the antenna aperture angles and gain values for the received signals. When running the script for the first time, it will automatically save a .mat file with the calculated specular reflection points position during the simulated time interval. The sp_file can be loaded when using the same time interval in future occasions.

### Executing program

Edit the input parameters in config.m to fit your simulation.
* Define the desired time interval with the appropiate time resolution. The time resolution has to be an integer value in seconds.
* Set the correct paths to the rinex navigation file and LEO satellite trajectory file.
* Define the length of the time limited register, this is for the time limited satistical analysis.
* If there is a previous sp_locations file for this time interval and resolution, it can be loaded by setting the load_sp flag to 1.
* Also, if there is a receiver antenna gain pattern file for different values of elevation and azimuth, its path has to be set as input as well.

After configuration, run the complete main.m script.


## Acknowledgments

The script SpecularReflectionPoint.m was developed by Lucas Sanz, 2021, UIDET-SENyT, Facultad de Ingenieria, UNLP.