This folder includes functions that return structs which store all required input parameters for specific functions throughout the repository.

#### msrCookieFinal
Includes the struct variable storing all required inputs for the real satellite measurement model. The satellite model is used to simulate the real attitude at any point during orbit.

#### msrCookieFinalExtended
Includes the struct variable storing all required inputs for the satellite measurement model used for MEKF correct function.

#### stateTransCookieFinalNominal
Includes the struct variable storing all required inputs for the satellite state model. The satellite model is used for the MEKF predict function and to simulate the real attitude at any point during orbit.
