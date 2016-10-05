# iceflow

This module is designed to calculate the ice sheet dynamics (viscosity, velocity, etc.)
given inputs. It is built as a stand alone package that can be used for dynamics
calculations inside of full ice sheet models. 

## Usage

All dynamics information for a given modeling domain is stored inside of 
an iceflow object: 

`type(iceflow_class) :: flow`

The module consists of three public subroutines that control the variables
of the iceflow object, for initialization, updating and termination, respectively: 
```
call iceflow_init(flow,...)
call iceflow_update(flow,dt,...)
call iceflow_end(flow)
```


