#!/bin/bash
# After the expansion solve array finishes (squeue empty), refresh the whole
# downstream chain from the enlarged dataset. Run ON THE CLUSTER head node:
#   ~/cfd/bin/micromamba run -p ~/cfd/env bash refresh_after_expansion.sh
# then locally:
#   scp elazar@100.70.186.114:cfd/disc/cfd/bigsweep/{dataset.csv,coeff_table.csv,surrogate_report.json} cfd/bigsweep/
#   python3 cfd/bigsweep/export_surrogate_js.py      # rewrites cfd_surrogate.js (bias-corrected)
#   python3 cfd/verify_flight.py && node cfd/verify_flight.js   # JS/Python parity
#   python3 scripts/fit_cfd_flight_numbers.py        # refit flight numbers on new surrogate
#   python3 cfd/verify_anchors_flight.py             # anchor-disc trajectory sanity
set -e
BS="$(cd "$(dirname "$0")" && pwd)"
python "$BS/collect.py"
python "$BS/fit_surrogate.py"
echo "done: dataset.csv, coeff_table.csv, surrogate_report.json refreshed"
