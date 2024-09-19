# CryoSim
## Introduction
CryoSim is a model of a generic cylindrical 4-K cryostat cooled with a commercial pulse tube cryocooler (PTC), extensively parametrised and setup to conduct a joint mechanical and thermal analysis via finite element methods (FEMs).  
Such test cryostat can be used to characterise optical components and full reimaging optical systems dedicated to CMB observations.  
More information can be found in the ralated publication:
- Thomas J. L. J. Gascard, Yi Wang, Jon E. Gudmundsson, Eve M. Vavagiakis, Cody J. Duell, Zachary B. Huber, Lawrence T. Lin, Michael D. Niemack, and Rodrigo G. Freundt "Thermal and mechanical study of a parametrised cryostat model for optical characterisation of upcoming CMB experiments", Proc. SPIE 13102, Millimeter, Submillimeter, and Far-Infrared Detectors and Instrumentation for Astronomy XII, 131022O (16 August 2024); https://doi.org/10.1117/12.3019995

## CAD Model
The fully parametrised CAD model, built in Solidworks 2022 SP5.0, is available [here](https://www.dropbox.com/scl/fo/h18vc017o2419hdnloc0h/AHrIsl6PEw-ja6f75JGTkvk?rlkey=ck93zx9ar5skr54sd0vqw96nw&st=sw0y755q&dl=0):  
Parameters are listed under the "equation" section of the Feature Manager Design tree.  
The overall design geometry is driven by the following parameters:
- the optics tube (OT) length,
- the optics tube diameter,
- the position of the focal plane within,
- the spacing between shells,
- the position of the PTC,
- the position of the DR.  

For instance, changing the OT length in the equation manager will automatically update the 4K, 40K and 300K lengths as well as re-adjust the DR and PTC connections.

## Mechanical analysis
Solidowrks can run FEM mechanical analysis to check the stress and the main accoustic resonnances of the system under study.
Parameters may be used as long as they do not modify the number of fixtures and connections establmished for such process.
Otherwise, the user will need to go through their definition again before trying to run a simulation.  
Regardless, it is highly recommended to work on a separate set of parts and assemblies to conduct such analysis so as to ensure the baseline design is preserved for the thermal analysis.  
An example is provided here: link?

## Thermal analysis
The thermal analysis is conducted in COMSOL Multiphysics 6.1 (Build 357).  
The thermal model for CryoSim and a quick start guide on how to set it up are available in the [comsol](/comsol) folder.  
This analysis is heavily dependant on a database, available under the ["material db"](/material_db/db_files) folder, that provides the various propoerties for the materials involved and the pulse tube colling power curves, in CSV files.  
A [Jupyter notebook](/material_db/thermal_prop_calc.ipynb) provides examples on how the database is built using the [thermal_props](/material_db/thermal_props.py) Python module (Python 3.11.3).  

