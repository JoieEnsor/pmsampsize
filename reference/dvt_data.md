# Simulated Deep Vein Thrombosis (DVT) Dataset

A simulated dataset of 9,805 patients from seven studies, used to
illustrate the development of a clinical prediction model for diagnosing
deep vein thrombosis (DVT). All predictors are binary (0/1).

## Usage

``` r
dvt_data
```

## Format

A data frame with 9,805 rows and 12 columns:

- male:

  Male sex (1 = male)

- OCuse:

  Oral contraceptive use (1 = yes)

- malignancy:

  Active malignancy/cancer (1 = yes)

- surgery:

  Recent surgery (1 = yes)

- ab_leg_trauma:

  Recent leg trauma (1 = yes)

- vein_distension:

  Leg vein distension/swelling (1 = yes)

- calf_diff:

  Calf width difference \>3 cm (1 = yes)

- ab_ddimer:

  Abnormal D-dimer blood test result (1 = yes)

- pregnancy:

  Recent pregnancy (1 = yes)

- immobile:

  Recent immobility \>24 hours (1 = yes)

- clotting:

  Recent anticoagulant use (1 = yes)

- dvt:

  DVT diagnosis (1 = DVT present; outcome variable)

## Source

Simulated dataset provided for illustrative purposes.

## Details

Outcome prevalence is approximately 21%. The dataset spans seven studies
of varying size. It is intended as a realistic example of a modest-sized
development dataset for use with the generalised sample size approach
([`pmgen_reference()`](https://joieensor.github.io/pmsampsize/reference/pmgen_reference.md),
[`pmgen_simulate()`](https://joieensor.github.io/pmsampsize/reference/pmgen_simulate.md)).
