A sweet code.

1D homogenous Monte Carlo for 22.211 PSet 1

Program Steps

- Particle born
	- Sample Location (0 - 6cm), set region

- Sample Flight Direction
- Sample flight distance

- Move Particle, evaluate
	- if outside -> DONE
	- if in new region, set region

- Sample interaction
	- if absorbtion, tally absorption -> DONE
	- else scatter:
		goto: Sample Flight Direction

