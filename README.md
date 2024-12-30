# Levicoustic 🍃

#### ⟾ Numerical simulation of acoustic levitation thanks to the matrix method

This project is part of my TIPE (travail d'initiative personnelle encadré) on the subject of acoustic levitation.

To test the various aspects of this phenomena and to see what's not visible, a simulation was necessary.

## Features

- Compute the [pressure, velocity, acoustic radiation potential, radiation force on small sphere] fields for any reflection number with good approximation and low compute power
- Handle multiples transducers with their inclination/phase/offset
- Graph field at one point on the x-axis for better visualisation

## Credit

This repository implements the matrix method used in the publication ["Matrix Method for Acoustic Levitation Simulation"](https://www.researchgate.net/publication/224254694_Matrix_Method_for_Acoustic_Levitation_Simulation).

The core of this repo is heavely inspired by [slkiser's](https://github.com/slkiser) implementation of the matrix method done in python, it only covers the section II of the paper (pressure field in the simplest scenario) but is an excellent start to understand how the algorithm works outside of theory.

> [1] S. L. Kiser, Matrix method in Python for acoustic levitation simulations. [Online]. Available: https://github.com/slkiser/acousticLevitation
