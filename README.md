# Programming Project

This is a programming project for the M2 TCCM unit Computer Programming Project. 
The aim is to implement a simulation of a Lennard-Jones fluid of atoms that runs in parralel.

## Compiling and running the code

    make 
    ./Project -i <input> -o <output [options]

## System of units

| Variable | Length   | Time |  Mass | Energy | Pressure |Temperature | Angle
|--|--|--|--|--|--|--|--|
| **Unit** | bohr  | picosecond | Dalton | Joules | Pascal | Kelvin | °

## Building the input file

  **THE INPUT FILE IS CASE SENSITIVE**
  **EVERYTHING AFTER A "!" IS CONSIDERED AS A COMMENT**

Some of the input options are blocks, each block will end with the **End** keyword.

### Specifying the units

It is not mandatory to provide the informations of the input file in the default units. It is possible to specify the unit by adding one letter after the value.
|Unit| Angstrom | °C | °F |
|--|--|--|--|--|
|**Letter**  | A | C | F |


### Specifying the box size

The keyword for specifying the periodic box dimensions is **Box**. It needs 6 arguments 

 - **a** the lengh of the box
 - **b** the width of the box
 - **c** the height of the box
 - **alpha** the angle between **b** and **c**
 - **beta** the angle between **a** and **c**
 - **gamma** the angle between **a** and **b**

It is not mandatory to provide all the angles, they will be considered to be 90° if they are not provided.

Note that the values of **alpha**, **beta** and **gamma** have to comply to the [allowed  values for the triclinic unit-cell angles ](https://doi.org/10.1107/S0108767310044296)

*Example*

    Box
	    a 3.5512 A
	    c 6.2115
	    b 4.1123
	    alpha 65
	    beta 115
    End
*This input will lead to a box of size 6.7108 x 6.2115 x 4.1123 (bohr)  with angles of 65, 115 and 90°*

### Specifying the atoms

**Random generation**

**Reading from file**

## Theoretical background

### Velocity verlet
$$
\vec{r}(t + \delta t) = \vec{r} (t) + \delta t \vec{v}(t) + \frac{1}{2} \delta t^2 \vec{a}(t) \\
\vec{v}(t + \delta t ) = \vec{v}(t) + \frac{1}{2} \delta t \left[ \vec{a}(t) + \vec{a}(t + \delta t)  \right]
$$

### Radial distribution function

$$
g(r) = \frac{\text d n_r}{4 \pi r^2 \text d r  \rho(r)}
$$
