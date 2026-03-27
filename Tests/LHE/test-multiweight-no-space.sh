#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
src_root="$(cd "${script_dir}/../.." && pwd)"
install_root="$(cd "${script_dir}/../../../../" && pwd)"
multiweight_dir="${src_root}/Contrib/MultiWeight"

# Prefer the module environment when available; the default repository pulls
# optional plugins from the loaded Herwig stack.
if [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck source=/etc/profile.d/modules.sh
  . /etc/profile.d/modules.sh
  module load herwig/unstable-full >/dev/null 2>&1 || true
fi

if command -v Herwig >/dev/null 2>&1; then
  herwig_bin="$(command -v Herwig)"
else
  herwig_bin="${install_root}/bin/Herwig"
fi

if [ ! -x "${herwig_bin}" ]; then
  echo "Could not find Herwig binary at ${herwig_bin}" >&2
  exit 1
fi

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/herwig-lhe-no-space.XXXXXX")"
cleanup() {
  status=$?
  if [ "${status}" -eq 0 ]; then
    rm -rf "${tmpdir}"
  else
    echo "Regression artifacts kept in ${tmpdir}" >&2
  fi
  exit "${status}"
}
trap cleanup EXIT

make -C "${multiweight_dir}" MultiWeight.so >/dev/null

ln -s "${script_dir}/lhe-weights-no-space.lhe" "${tmpdir}/lhe-weights-no-space.lhe"
ln -s "${script_dir}/LHE-MultiWeight-NoSpace.in" "${tmpdir}/LHE-MultiWeight-NoSpace.in"
ln -s "${multiweight_dir}/MultiWeight.so" "${tmpdir}/MultiWeight.so"

cd "${tmpdir}"
"${herwig_bin}" read LHE-MultiWeight-NoSpace.in > read.log 2>&1
"${herwig_bin}" run LHEWeightsNoSpace.run -N 1 > run.log 2>&1

count="$(grep -c '^[[:space:]]*it = .*->' run.log)"
test "${count}" -eq 3

grep -F 'it = scale_up->' run.log >/dev/null
grep -F 'it = k3_0_k4_0->' run.log >/dev/null
grep -F 'it = k3_1_k4_m1->' run.log >/dev/null

if grep -F 'it = rwgt_' run.log >/dev/null; then
  echo "Found raw rwgt ids in output; expected mapped weight names." >&2
  exit 1
fi

awk -F'->' '/it = scale_up->/ { ok = ($2 + 0 == 1.0) } END { exit ok ? 0 : 1 }' run.log
awk -F'->' '/it = k3_0_k4_0->/ { ok = ($2 + 0 == 2.0) } END { exit ok ? 0 : 1 }' run.log
awk -F'->' '/it = k3_1_k4_m1->/ { ok = ($2 + 0 == 3.0) } END { exit ok ? 0 : 1 }' run.log

echo "PASS: no-space weight tags are parsed and mapped correctly."
