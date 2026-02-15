# utils.sh



# ---------------------------------------
# Get number of CPU cores available
# Priority:
#   1) SLURM_CPUS_PER_TASK
#   2) scontrol show job
#   3) fallback = 1
# ---------------------------------------
get_cpu_cores() {
  local cores=""

  # 1) SLURM: explicit allocation
  if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    if [[ "${SLURM_CPUS_PER_TASK}" =~ ^[0-9]+$ ]]; then
      cores="${SLURM_CPUS_PER_TASK}"
    fi
  fi

  # 2) SLURM: query job info
  if [[ -z "$cores" && -n "${SLURM_JOB_ID:-}" && -x "$(command -v scontrol)" ]]; then
    cores="$(
      scontrol show job -o "$SLURM_JOB_ID" 2>/dev/null \
        | grep -oE 'CPUs/Task=[0-9]+' \
        | cut -d= -f2
    )"
  fi

  # fallback
  if [[ -z "$cores" || "$cores" -lt 1 ]]; then
    cores=1
  fi

  echo "Using ${cores} cores for processing." >&2
  echo "$cores"
}

# ---------------------------------------
# Get usable memory (GB) for processing
# Priority:
#   1) SLURM_MEM_PER_NODE
#   2) scontrol show job (TRES mem)
#   3) /proc/meminfo fallback
# ---------------------------------------
_get_mem_gb() {
  local s="$1"
  [[ -n "$s" ]] || return 1
  s="$(echo "$s" | tr '[:upper:]' '[:lower:]' | xargs)"

  # pure number -> MB (SLURM sometimes uses MB)
  if [[ "$s" =~ ^[0-9]+$ ]]; then
    echo $(( s / 1024 ))
    return 0
  fi

  if [[ "$s" =~ ^([0-9]+(\.[0-9]+)?)([a-z]+)$ ]]; then
    local val="${BASH_REMATCH[1]}"
    local unit="${BASH_REMATCH[3]}"
    case "$unit" in
      g|gb) awk "BEGIN{print int($val)}" ;;
      m|mb) awk "BEGIN{print int($val/1024)}" ;;
      k|kb) awk "BEGIN{print int($val/1024/1024)}" ;;
      t|tb) awk "BEGIN{print int($val*1024)}" ;;
      *) return 1 ;;
    esac
    return 0
  fi
  return 1
}

get_memory() {
  local mem_gb=""

  # 1) SLURM explicit
  if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    mem_gb="$(_get_mem_gb "$SLURM_MEM_PER_NODE" || true)"
  fi

  # 2) SLURM job info
  if [[ -z "$mem_gb" && -n "${SLURM_JOB_ID:-}" && -x "$(command -v scontrol)" ]]; then
    local raw
    raw="$(scontrol show job -o "$SLURM_JOB_ID" 2>/dev/null \
          | grep -oE 'mem=[^, ]+' \
          | cut -d= -f2)"
    [[ -n "$raw" ]] && mem_gb="$(_get_mem_gb "$raw" || true)"
  fi

  # 3) Fallback: /proc/meminfo (available)
  if [[ -z "$mem_gb" && -r /proc/meminfo ]]; then
    mem_gb="$(awk '/MemAvailable:/ {print int($2/1024/1024)}' /proc/meminfo)"
  fi

  if [[ -z "$mem_gb" || "$mem_gb" -le 0 ]]; then
    echo "Warning: Unable to determine system memory, defaulting to 8g" >&2
    echo "8g"
    return 0
  fi

  local usable
  usable="$(awk "BEGIN{print int($mem_gb*0.6)}")"
  [[ "$usable" -ge 1 ]] || usable=1

  echo "Using ${usable}g memory for processing." >&2
  echo "${usable}g"
}

# ---------------------------------------
# Filter out empty files
# ---------------------------------------
drop_empty_files() {
  local -n _paths_ref="$1" 
  local remove_empty="${2:-0}"

  local p
  local -a empty=()

  # helper: check if FASTQ / FASTQ.GZ is empty after decompression
  _is_empty_fastq() {
    local f="$1"
    [[ -e "$f" ]] || return 0
    if [[ "$f" == *.gz ]]; then
      if ! gzip -cd "$f" 2>/dev/null | head -c 1 | grep -q .; then
        return 0
      else
        return 1
      fi
    fi
    [[ ! -s "$f" ]] && return 0 || return 1
  }

  local -a kept=()
  for p in "${_paths_ref[@]}"; do
    if _is_empty_fastq "$p"; then
      empty+=("$p")
    else
      kept+=("$p")
    fi
  done


  _paths_ref=("${kept[@]}")
  if [[ "$remove_empty" -eq 1 && "${#empty[@]}" -gt 0 ]]; then
    echo "Warning: The following files are empty and have been excluded:"
    for p in "${empty[@]}"; do
      echo "  $p"
      rm -f "$p"
    done
  fi
  return 0
}


# ---------------------------------------
# Parse FASTQ filenames and group R1/R2 by prefix
# Output:
#   - FWD_BY_PREFIX[prefix] = R1 fastq
#   - REV_BY_PREFIX[prefix] = R2 fastq
# Return codes:
#   0 success
#   2 duplicated R1/R2
#   3 suffix mismatch
#   4 missing mate
# ---------------------------------------
get_fastq_prefix() {
  local -n _files_ref="$1"     # input: array of fastq paths
  local -n _fwd_ref="$2"       # output: assoc array (prefix -> R1)
  local -n _rev_ref="$3"       # output: assoc array (prefix -> R2)
  local ext="${4:-.fastq.gz}"

  local bad=0
  local b core prefix lane chunk read

  _fwd_ref=()
  _rev_ref=()

  for p in "${_files_ref[@]}"; do
    b="$(basename "$p")"
    [[ "$b" == *"$ext" ]] || { echo "ERROR: $b doesn't match suffix $ext" >&2; bad=1; continue; }
    core="${b%$ext}"

    lane=""
    chunk=""
    read=""

    # chunk: .001 / _001
    if [[ "$core" =~ ^(.+)[._]([0-9]{3})$ ]]; then
      core="${BASH_REMATCH[1]}"
      chunk="${BASH_REMATCH[2]}"
    fi

    # read: R1 / R2 or 1 / 2
    if [[ "$core" =~ ^(.+)[._]R([12])$ ]]; then
      core="${BASH_REMATCH[1]}"
      read="${BASH_REMATCH[2]}"
    elif [[ "$core" =~ ^(.+)[._]([12])$ ]]; then
      core="${BASH_REMATCH[1]}"
      read="${BASH_REMATCH[2]}"
    else
      echo "ERROR: $b missing read token (R1/R2)" >&2
      bad=1
      continue
    fi

    # lane: L001
    if [[ "$core" =~ ^(.+)[._](L[0-9]{3})$ ]]; then
      core="${BASH_REMATCH[1]}"
      lane="${BASH_REMATCH[2]}"
    fi

    prefix="$core"
    [[ -n "$prefix" ]] || { echo "ERROR: empty prefix in $b" >&2; bad=1; continue; }

    if [[ "$read" == "1" ]]; then
      [[ -z "${_fwd_ref[$prefix]:-}" ]] || { echo "ERROR: duplicated R1 for $prefix" >&2; return 2; }
      _fwd_ref["$prefix"]="$p"
    else
      [[ -z "${_rev_ref[$prefix]:-}" ]] || { echo "ERROR: duplicated R2 for $prefix" >&2; return 2; }
      _rev_ref["$prefix"]="$p"
    fi
  done

  [[ "$bad" -eq 0 ]] || return 3

  for prefix in "${!_fwd_ref[@]}"; do
    [[ -n "${_rev_ref[$prefix]:-}" ]] || { echo "ERROR: missing R2 for $prefix" >&2; return 4; }
  done
  for prefix in "${!_rev_ref[@]}"; do
    [[ -n "${_fwd_ref[$prefix]:-}" ]] || { echo "ERROR: missing R1 for $prefix" >&2; return 4; }
  done

  return 0
}