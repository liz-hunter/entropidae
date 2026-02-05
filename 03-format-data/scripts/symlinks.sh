#!/usr/bin/env bash

counts_tsv="${1:-}"
fna_dir="${2:-}"
output_dir="${3:-}"

if [[ -z "${counts_tsv}" || -z "${fna_dir}" || -z "${output_dir}" ]]; then
  echo "ERROR: Missing required arguments." >&2
  echo "Usage: $0 <accession_counts.tsv> <fna_dir> <output_dir>" >&2
  exit 2
fi

if [[ ! -f "${counts_tsv}" ]]; then
  echo "ERROR: counts_tsv not found: ${counts_tsv}" >&2
  exit 2
fi

if [[ ! -d "${fna_dir}" ]]; then
  echo "ERROR: fna_dir not found: ${fna_dir}" >&2
  exit 2
fi

mkdir -p "${output_dir}"

missing_log="${output_dir}/missing_symlinks.tsv"
: > "${missing_log}"
printf "accession\tfilename\texpected_path\n" > "${missing_log}"

# Parse header to find column indices
header="$(head -n 1 "${counts_tsv}")"
IFS=$'\t' read -r -a cols <<< "${header}"

acc_col=0
file_col=0
count_col=0

for i in "${!cols[@]}"; do
  [[ "${cols[$i]}" == "accession" ]] && acc_col=$((i+1))
  [[ "${cols[$i]}" == "filename"  ]] && file_col=$((i+1))
  [[ "${cols[$i]}" == "count"     ]] && count_col=$((i+1))
done

if [[ "${acc_col}" -eq 0 || "${file_col}" -eq 0 || "${count_col}" -eq 0 ]]; then
  echo "ERROR: counts_tsv header must include 'accession', 'filename', and 'count' columns." >&2
  echo "Header was: ${header}" >&2
  exit 2
fi

created=0
missing=0

# Read TSV rows (skip header)
while IFS=$'\t' read -r -a fields; do
  accession="${fields[$((acc_col-1))]}"
  filename="${fields[$((file_col-1))]}"
  count="${fields[$((count_col-1))]}"

  # Check for missing rows
  [[ -z "${accession}" ]] && continue
  [[ -z "${filename}"  ]] && continue
  [[ -z "${count}"     ]] && continue

  # Check that counts are integers
  if ! [[ "${count}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Non-integer count for accession '${accession}': '${count}'" >&2
    exit 2
  fi

  src="${fna_dir}/${filename}"
  if [[ ! -e "${src}" ]]; then
    echo "WARNING: filename '${filename}' for accession '${accession}' was not found in ${fna_dir}" >&2
    printf "%s\t%s\t%s\n" "${accession}" "${filename}" "${src}" >> "${missing_log}"
    missing=$((missing+1))
    continue
  fi

  real_src="$(realpath "${src}")"
  for ((i=1; i<=count; i++)); do
    ln -sf "${real_src}" "${output_dir}/${accession}_dup${i}.fna"
    created=$((created+1))
  done

done < <(tail -n +2 "${counts_tsv}")

echo "Done."
echo "  Output dir: ${output_dir}"
echo "  Symlinks created: ${created}"
echo "  Missing files: ${missing}"
if (( missing > 0 )); then
  echo "  Missing log: ${missing_log}"
  exit 1
fi
