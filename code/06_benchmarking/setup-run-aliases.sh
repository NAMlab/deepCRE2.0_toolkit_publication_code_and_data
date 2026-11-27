#!/bin/bash
# setup-run-aliases.sh - Source this script to enable mpl-run command
#
# Usage: source setup-run-aliases.sh
# Then navigate to a runs/<model> directory and use mpl-run

# Hardcoded cache directory with realpath - works from any directory
CACHE_DIR="$(realpath "$(dirname "${BASH_SOURCE[0]}")/workflows/cwl/containers")"

# Main mpl-run function
mpl-run() {
    local parallel=false
    local singularity=false
    local action="run"
    local extra_args=""

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -p|--parallel)
                parallel=true
                shift
                ;;
            -s|--singularity)
                singularity=true
                shift
                ;;
            -ps|-sp)
                parallel=true
                singularity=true
                shift
                ;;
            reset)
                action="reset"
                shift
                ;;
            -h|--help)
                _mpl_run_help
                return 0
                ;;
            *)
                # Capture any additional arguments to pass through to cwltool
                extra_args="$extra_args $1"
                shift
                ;;
        esac
    done

    # Check if we're in a valid run directory
    if [ ! -f "run.cwl" ]; then
        echo "Error: No run.cwl detected in current folder."
        echo "Are you in a runs/<model> folder?"
        echo ""
        echo "Current directory: $(pwd)"
        return 1
    fi

    # Handle reset action
    if [ "$action" = "reset" ]; then
        echo "This will delete the contents of the folders dataset and logs."
        read -p "Do you want to proceed? (y/N) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Reset cancelled."
            return 0
        fi
        echo "Resetting run directory..."
        # Remove log files
        find logs/* ! -name 'run' -type d -exec rm -rf {} + 2>/dev/null
        find logs/ ! -name '.gitkeep' -type f -exec rm -f {} + 2>/dev/null
        find logs/run/ ! -name '.gitkeep' -type f -exec rm -f {} + 2>/dev/null
        # Empty the dataset directory
        find dataset/* -type d -exec rm -rf {} + 2>/dev/null
        find dataset/ ! -name '.gitkeep' ! -name '.gitignore' -type f -exec rm -f {} + 2>/dev/null
        echo "Reset complete."
        return 0
    fi

    # Build cwltool command
    local cmd="cwltool"

    if [ "$singularity" = true ]; then
        # Set Singularity environment variables
        export SINGULARITY_CACHEDIR="$CACHE_DIR"
        export SINGULARITY_TMPDIR="$CACHE_DIR/tmp"
        export SINGULARITY_PULLDIR="$CACHE_DIR/pull"
        export CWL_SINGULARITY_CACHE="$CACHE_DIR"
        cmd="$cmd --singularity"
    fi

    if [ "$parallel" = true ]; then
        cmd="$cmd --parallel --strict-memory-limit"
    fi

    cmd="$cmd --outdir dataset --log-dir logs --timestamps --disable-color --enable-ext"
    cmd="$cmd$extra_args"
    cmd="$cmd run.cwl run.yml"
    cmd="$cmd > logs/run/run.log 2> logs/run/run.error.log"

    echo "Running: $cmd"
    eval "$cmd"
}

# Help function
_mpl_run_help() {
    cat << EOF
Usage: mpl-run [OPTIONS] [reset]

Run CWL workflows from within a runs/<model> directory.

This command must be run from a directory containing run.cwl and run.yml
(typically runs/deepcre, runs/deepcis, etc.).

Options:
    -p, --parallel          Enable parallel execution
    -s, --singularity       Use Singularity instead of Docker
    -ps, -sp                Enable both parallel and Singularity
    reset                   Clean logs and dataset directories
    -h, --help              Show this help message

Additional arguments:
    Any other arguments are passed directly to cwltool.
    Example: mpl-run -ps --parallel-max 5

Examples:
    cd runs/deepcre
    mpl-run                          # Docker, sequential
    mpl-run -p                       # Docker, parallel
    mpl-run -s                       # Singularity, sequential
    mpl-run -ps                      # Singularity, parallel
    mpl-run reset                    # Clean logs and dataset
    mpl-run -ps --parallel-max 5     # Singularity, maximum 5 parallel jobs (see cwltool documentation)

Convenience aliases (also available):
    mpl-run-p      Same as mpl-run -p
    mpl-run-s      Same as mpl-run -s
    mpl-run-ps     Same as mpl-run -ps
EOF
}

# Convenience aliases
alias mpl-run-p='mpl-run -p'
alias mpl-run-s='mpl-run -s'
alias mpl-run-ps='mpl-run -ps'

# Export function for use in subshells
export -f mpl-run
