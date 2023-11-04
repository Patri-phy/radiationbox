## My Project: Radiation outgoing from a box

### Main topic: time evolution of the intensity of the photons outgoing
### from the square window of a face of a hollow cubic body.
### In the body is present a source of monochromatic electromagnetic
### radiation. The light is impulsive at the time t=0.

# In the problem one supposes that:
1. the source is point-like and isotropic; it is placed in the center of the cube;
2. the inner part of the body is filled with a gas with isotropically and elastically
diffusing atoms;
3. every atom of the inner surface of the body can either scatter (in an isotropic and
elastic way) or absorb photons;
4. each photon behavoir is independent from each other photon behavoir.

The cross section magnitude is fixed such that the mean free path is greater than the length of
the cube edge (lower limit of the mean free path in the case of molecular regime).
Each particle position is updated by a time-dependent Monte Carlo algorithm.

### Secondary topic: calculation of the escape factor via statistical weights MC method.

## Presentation

1. Introduction (theoretical feature: diffusion, scattering and absorption)
2. Definition of the problem (description of the set-up and illustration)
3. Method of implementation (algorithm, illustration of the algorithm)
4. Results and discussion (variation of the parameter (intensity)

## Algorithm structure

1. Parameter settings, variable type declaration, variable initialization
2. Geometry implementation
3. Physical behavoir implementation
4. Data elaboration
