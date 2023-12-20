# Project Title

A self-built c++ fluid simulator with adaption of all kinds of Total Variation Diminishing (TVD) convection schemes.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Demo](#demo)
<!-- [Contributing](#contributing) -->
<!-- [License](#license) -->
<!-- [Acknowledgements](#acknowledgements) -->

## Features

- Any fluxe limiter can be applied within TVD framework.

## Installation

<!-- Provide instructions on how to install and set up your project. Include any dependencies and steps needed for a successful installation. -->

```bash
# You can utilize wmake to compile all c++ files 
cd solver
wmake
```

## Usage

<!-- Explain how to use your project. Include examples or code snippets if applicable. -->
  
```bash
cd run
source source_to_run
cd ../postplot
python newplot.py
```
## Demo
![Vmag_demo](https://github.com/Furkath/TVD_Fluid_Simulator/blob/master/results/320_Re1000_UMIST/pic_mag.png)
