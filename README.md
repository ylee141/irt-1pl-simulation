# irt-1pl-simulation

Simulating and estimating a one-parameter logistic (1PL) item response theory model in R.

The project simulates students' answers to test items from a known 1PL model, then estimates the item difficulties and the spread of ability in two ways:

- **Joint maximum likelihood (JML):** estimates every student's ability together with the item difficulties, with a bias correction.
- **Marginal maximum likelihood (MML):** integrates abilities out over a normal distribution, approximating the integral with simulation draws.

The MML estimates and their standard errors are checked against the `mirt` package.

## Files

| File | What it does |
|---|---|
| `01_simulate_data.R` | Simulates 500 students answering 10 items, with known abilities and difficulties |
| `02_estimate_jml.R` | Estimates the model by JML, with the bias correction and standard errors |
| `03_estimate_mml.R` | Estimates the model by MML, with standard errors, and compares the results with `mirt` |
| `report/1pl-report.qmd` | The report (see below) |

## The report

`report/1pl-report.qmd` is a study note that builds the 1PL model from the problem it solves, explaining why the model takes its form as well as how to estimate it.

- **Introduction:** why raw scores cannot compare students who took different test sets, and how placing ability and difficulty on one scale solves this.
- **Modeling:** the model derived step by step from what a probability of a correct answer should satisfy, the item characteristic curve, and a latent variable formulation that makes the role of luck explicit, with everyday analogies.
- **Assumptions:** each stated as an equation, explained in words, and paired with an example of how it can fail.
- **Identification:** why only the gap between ability and difficulty is identified, and the two conventions for fixing the scale's zero and unit.
- **Estimation:** JML and MML, why JML's difficulty estimates are stretched (the incidental parameter problem), and how MML avoids it.
- **Results and Discussion:** the simulation comparison, and what comes next.
- **Side note:** how the 1PL relates to discrete choice models in economics.

## How to run

`02_estimate_jml.R` and `03_estimate_mml.R` each run `01_simulate_data.R` first to create the simulated data, so either can be run on its own. The check in `03_estimate_mml.R` needs the `mirt` package.

