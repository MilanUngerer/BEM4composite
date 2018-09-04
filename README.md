# BEM_microwire

## Abstract

In the last time composite materials have been subject of multiple studies due to the great utility of the special properties that they present. In this work the study was oriented to a material composed of a dielectric matrix with inclusions of conductive materials called ferromagnetic micro-wires. These wires present the peculiarity of changing their electromagnetic properties in front of an external stimulus. When applying a magnetic field or a mechanical tension the permittivity and permeability change, thus changing the way in which the incident electromagnetic radiation is scattered.

This work is oriented to find the optimal configuration of the micro-wires to maximize the change of signal in order to design a stress sensor for dielectrics. The radiation will be of the microwave type and a mechanic tension stimulus will be used. In order to reach the proposed objectives a program was developed that allows to calculate the electromagnetic field scattered by the composite at any point in which it is desired to position the antenna.

For the calculation of the field the Helmholtz equation was used through its boundary integral formulation which allows a simpler manipulation for the calculate of the field and to realize the interaction between the different obstacles. Because this equation does not have an analytical solution it is necessary to approximate it by means of some numerical method, in this work we used the Boundary Element Method with, the library for Python, BEM ++. 

