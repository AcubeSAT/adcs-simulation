This repository is the Attitude, Orbit Determination and Control Subsystem (ADCS/AOCS) simulation analysis of the AcubeSAT project, part of the SpaceDot team. The general design approach, theoretical implementation and simulation analyses conducted by the subsystem can be found in the Critical Document Review (CDR) of AcubeSAT, in the Design Definition & Justification File of AOCS (DDJF_AOCS). Further documents composed by the AOCS members to describe and analyse certain design and theoretical aspects can be found in our documentation platform.

#### Useful links
+ [SpaceDot Website](https://spacedot.gr/)
+ [AcubeSAT Website](https://acubesat.spacedot.gr/)
+ [ADCS Website](https://acubesat.spacedot.gr/subsystems/attitude-determination-control/)
+ [DDJF AOCS](https://gitlab.com/acubesat/documentation/cdr-public/-/blob/master/DDJF/DDJF_AOCS.pdf)
+ [AOCS Documentation](https://helit.org/mm/docList/public/ADC)

## Design Definition & Justification File of Attitude Orbit Control Subsystem

The [Design Definition Justification File (DDJF)](https://gitlab.com/acubesat/documentation/cdr-public/-/blob/master/DDJF/DDJF_AOCS.pdf) of AOCS is the document aiming at describing in detail all the aspects regarding attitude control of the satellite for the AcubeSAT mission. The mission includes a scientific payload, which is analyzed through photographs, requires a directional patch antenna as explained in [DDJF_PL](https://gitlab.com/acubesat/documentation/cdr-public/-/blob/master/DDJF/DDJF_PL.pdf) and [DDJF_TTC](https://gitlab.com/acubesat/documentation/cdr-public/-/blob/master/DDJF/DDJF_TTC.pdf). Consequently, the AOCS is required to achieve nadir pointing during satellite nominal operations in a safe and orderly fashion. This is achieved by using
active attitude determination and control techniques.

This document provides information on the satellite determination and control approach with regards to both hardware and software according to ESA guidelines. It also contains the analysis approaches, along with their justification and description, as well as the required documentation, all in compliance with the ESA Policy. This document will remain open until the later stages of the project and therefore will be updated on a regular basis during the different FYS phases.

## Physical Architecture
#### Sensors
The AOCS Physical Architecture involves six sensors, two different types of actuators and one MCU board.

##### Magnetometers
Two 3-axis magnetometers are used in our AOCS unit. One of them functions as the main sensor in order to acquire high precision geomagnetic field measurements, while the other is equipped only for redundancy purposes, in case a failure occurs. The RM3100 provided by PNI is selected as the precision magnetometer. It consists of three coil sensors (Sen XY/Z), which are driven by the MagI2C controller in a separate IC. Τhe second module is provided by ISIS and is embedded in the iMTQ magnetorquers
board, with an accuracy of < 3 µT. 

The RM3100 is utilized consistently, while the second magnetometer, provided by ISIS operates periodically, to verify the precision of the RM3100, essentially using them simultaneously using the warm-redundancy concept.

##### Gyroscope
The ADXRS453 three-axis MEMS gyroscope is utilized as the sole gyroscope unit, provided by Analog-Devices.

##### Sun Sensor
Our design uses Coarse sun sensors, which are integrated on the Solar Panels, provided by EnduroSat. However one additional coarse sun sensor, the TSL2540 is utilized, provided by ams AG. Essentially completing the required 6 unit sun sensor configuration, having one on each satellite side.

##### Temperature sensor
Two MCP9808 temperature sensors are used, provided by Microchip Technologies.

#### Actuators
##### Magnetorquers
The ISIS iMTQ magnetorquer board has been selected, containing three magnetorquers, one for each satellite body frame axis. The board is equipped with two magnetorquer rods and one air-core magnetorquer. They are connected to the ADCS/OBC board and are controlled through the AOCS MCU.

##### Reaction Wheel
The RW210 reaction wheel is used, manufactured by Hyperion Technologies.

#### AOCS PCB
##### CAN transceiver
A CAN-FD transceiver is necessary to shift the CAN bus voltages and establish the physical architecture of the bus. The one used for the AOCS and the rest of the subsystems is Texas Instruments TCAN337G. It is compliant with the CAN-FD bus requirement of the SatNOGS COMMS board, and offers protection features including:
+ Overvoltage, undervoltage and overcurrent protection
+ Fault detection output
+ Babbling idiot protection (via bus timeout detection)
Interfacing between AOCS and OBDH is achieved by the on‐board CAN bus, which is used to transfer parameters, events, commands, science and logging data, TT&C packets and synchronisation markers/heartbeats. The CAN bus is described in detail in [DDJF_OBDH](https://gitlab.com/acubesat/documentation/cdr-public/-/blob/master/DDJF/DDJF_OBDH.pdf).

##### MCU
The SAMV71Q21RT radiation-tolerant MCU, manufactured by Microchip Technologies, is placed in the ADCS/OBC board to fulfill the computational needs of the AOCS.

It is a RT microcontroller (which is a drop-in compatible non-RT version) with a high IO count and number of available peripherals.

This model was chosen since it is one of the few commercially available MCU featuring
the modern ARM Cortex-M architecture as of 2020, while supporting the performance and interfacing needs of the subsystem. This solution was chosen against a non-RT MCU due to the increased reliability it offers, while still being compatible with any design changes or constraints requiring a non-RT solution, either by AcubeSAT or any users of its open-source design. 

Features common to the ARM architecture or specific to the SAMV71Q21 include:
+Dedicated low power modes
+ Integrated RTC functionality
+ An internal temperature sensor
+ Interfaces: I2C, SPI, CAN-FD, UART
+ ECC in the embedded flash memory
+ 2 independent internal watchdogs. Due to the Radiation Tolerance of this part, no external watchdogs are used.
+ Full support for debugging features
+ A large degree of hardware & software support, including a dedicated development board (ATSAMV71-XULT), IDE and configuration tools.

The Attitude & Orbit Control System (AOCS) MCU is hosted on the ADCS/OBC board, which is based on the PC/104 standard. Due to a lack of available vertical space, the board hosts both of the aforementioned subsystems, each occupying one side of the PCB.

## Justification of the Design
AcubeSAT's AOCS physical architecture design approach aims for a low-cost, robust and reliable AOCS system, using components with relatively low power consumption, small size and high accuracy.

Regarding the sensor selection, the main focus was being able to dissipate the satellite angular velocity autonomously, without ground segment support, by utilizing the magnetometer. Specifically, the magnetometer is considered reliable, flexible and presents minimal restrictions, due to its capability of measuring the magnetic field with high accuracy, its small size and low power consumption. Owing to the robustness of that sensor, we opted for having a second redundant magnetometer. Additionally, in order to unambiguously define the satellite orientation, measurements from two sensors, namely the magnetometer and the sun sensor, are combined. The sun sensor is deployed only during AOCS nominal mode, complementing the magnetometer measurements. A MEMS gyroscope is deployed, enabling the satellite's accurate angular velocity definition in attitude determination, during nominal mode.

Regarding actuator selection, we aimed to safely and uninterruptedly achieve the desired pointing accuracy using two low-cost and small-size actuators. Magnetorquers are considered a reliable and flexible way of attitude control in LEO, where the magnetic field is considered sufficient for the control of the satellite. During the AOCS detumbling mode, the magnetorquers are solely used. However, the torque of the magnetorquers is perpendicular to the aforementioned magnetic field, depriving the controller of one degree of freedom. In order to tackle this, a reaction wheel was selected to compensate the ancillary torque demanded to be produced by the magnetorquer on the Z axis of the body frame. This allows the satellite to achieve the desired pointing accuracy in AOCS nominal mode.

We opted for using a MCU separate from the OBDH subsystem, enabling us to achieve simpler data distribution, greater memory storage, greater processing power, while also
offering redundancy. The integrated MCU on the ISIS MTQ board is not utilized in our design due to lack of software interface, thus it is used exclusively to calculate any required interface between the commanded torque and the MTQ. The MCU is accompanied with its own internal independent watchdog, which ensures that any latent failure in the software will be reliably detected.
