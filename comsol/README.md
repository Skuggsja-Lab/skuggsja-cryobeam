### COMSOL Model
An overview of COMSOL capabilities for such thermal analysis and a quick start on the CryoSim model are currently in progress and will be available shortly.  
For user who have experience with COMSOL, the CAD model is available [here](https://www.dropbox.com/scl/fo/h18vc017o2419hdnloc0h/AHrIsl6PEw-ja6f75JGTkvk?rlkey=ck93zx9ar5skr54sd0vqw96nw&st=sw0y755q&dl=0).  
In COMSOL, under the geometry branch, the Solidworks Livelink documlent folder will need to be changed accordingly to where the model was dowloaded. Adjustements in the geometry selections may also be needed.  
Once the geometry is properly defined, another important aspect is to properly setup the study you want to run. Select the first step for the concerned study and run a calculation of the initial values. Once this is done, the full solver setup becomes available.  
Maker sure the solver uses a fully coupled solve, and that the latter uses a highly non-linear Newtonian algorythm with 200 or so iterations for the termination criteria.
