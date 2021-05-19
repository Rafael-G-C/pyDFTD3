## Notes about this folder

- For simplicity when refactoring I have removed the Becke-Johnson damping
  scheme but it can be reintroduced.
- I have moved out all parameter and I/O out of the d3 function so that it does
  not have to know anything about functionals and setting the appropriate
  parameters can be done caller-side.
- It turned out to be significantly more efficient to compute all derivatives
  of one order in one shot. [demo.py](demo.py) shows how this can be done.
- Distance-based screening can be done but is not completely trivial and hasn't
  been attempted in this rewrite.
