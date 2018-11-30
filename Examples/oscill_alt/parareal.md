## FSI Probleme, Fragen
* Teilweise keine Konvergenz wenn θ nicht entsprechend gewählt
* Bereits bei h=0.01 keine Konvergenz für zu kleine Theta
* Probleme werden massiver für größere Teilintervalle
* Generelles Problem von impliziten Verfahren: Verhältnis der Schrittweiten entspricht im Allgemeinen nicht dem Verhältnis vom Rechenaufwand
#### Generelle Ideen
* Weitere Möglichkeiten zur Laufzeitreduktion von gröberem Verfahren
 - gröberes Gitter
 - Toleranz von groben Verfahren
### Parallel time stepping
- encapsulate fine and coarse method over one subinterval
- introduce theta-parareal with theta chosen as deviation from angle

### Numerical examples
- 2d example: works fine, theta stays high

 Finished parareal.
 Elapsed time is	695.128

 Time spent in coarse method: 573.224
 Time spent in fine method: 3522.5
