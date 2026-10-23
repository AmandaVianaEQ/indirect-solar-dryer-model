<h1>Indirect Solar Dryer Model </h1>

Simplified physical model of an indirect natural-convection solar dryer, applied to the drying of banana slices. Implemented in Python as a differential-algebraic (DAE) system.

Developed during a research internship at the LaTEP (Laboratoire de Thermique, Énergétique et Procédés), July–September 2026.

<h2>Context</h2>

The dryer had previously been studied only with CFD simulations, which are accurate but computationally heavy. The goal of this project was to build a lighter physical model that predicts the dryer's behaviour (air temperatures, airflow, product moisture over time) at a fraction of the computational cost, and to validate it against experimental data.

<h2>Model</h2>

The model couples the following blocks:

<table> 
  <tr>
    <th>Block</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>Solar collector</td>
    <td>Energy balance on the absorber plate, including thermal inertia</td>
  </tr>
  
  <tr>
    <td>Chimney effect</td>
    <td>Natural-convection airflow driven by the density difference</td>
  </tr>
  
  <tr>
    <td>Psychrometrics</td>
    <td>Moist-air properties (humidity, enthalpy, saturation)</td>
  </tr>
  
  <tr>
    <td>Drying chamber</td>
    <td>Heat and mass balances on the air and the product</td>
  </tr>
  
  <tr>
    <td>Convective transfer</td>
    <td>Heat transfer coefficients between air, plate and product</td>
  </tr>
  
  <tr>
    <td>Drying kinetics</td>
    <td>Page model for the product moisture ratio</td>
  </tr>

</table>

<h2>Numerical method</h2>
The coupled equations form a DAE system: differential equations for the dynamic states (plate temperature, product moisture) and algebraic equations for the air side.
Time integration with an explicit Euler scheme.
At each time step, the algebraic system is solved with Newton-Raphson.
The Page kinetic parameters were calibrated on experimental drying data using a Gauss-Newton least-squares method.
