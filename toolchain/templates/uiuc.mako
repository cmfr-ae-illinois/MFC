#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>
<%
if not gpu_enabled:
    raise RuntimeError("The UIUC H100 template requires --gpu.")
if not 1 <= tasks_per_node <= 8:
    raise RuntimeError("UIUC H100 nodes support 1 to 8 MPI ranks/GPUs per node.")
%>
% if engine == 'batch':
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${tasks_per_node}
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H100:${tasks_per_node}
#SBATCH --job-name="${name}"
#SBATCH --output="${name}.out"
#SBATCH --error="${name}.err"
#SBATCH --time=${walltime}
#SBATCH --partition=${partition or 'chess-gpu'}
#SBATCH --export=ALL
% if account:
#SBATCH --account=${account}
% endif
% if quality_of_service:
#SBATCH --qos=${quality_of_service}
% endif
% if email:
#SBATCH --mail-user=${email}
#SBATCH --mail-type=BEGIN,END,FAIL
% endif
% endif

${helpers.template_prologue()}

# GPU execution must take place inside an existing Slurm allocation.
if [ -z "$SLURM_JOB_ID" ]; then
    error "UIUC: request a GPU allocation first, or use -e batch."
    exit 1
fi

ok ":) Loading UIUC modules:\n"
cd "${MFC_ROOT_DIR}" || exit 1
. ./mfc.sh load -c uiuc -m g || exit 1
cd - > /dev/null || exit 1

export OMP_NUM_THREADS=1
ulimit -s unlimited

% if engine == 'interactive':
# An srun --pty shell already occupies a job step. MPI may use srun to
# bootstrap its processes, so allow sharing with that interactive shell.
export SLURM_OVERLAP=1
% endif

% for target in targets:
${helpers.run_prologue(target)}

% if mpi:
    (set -x; ${profiler} \
        mpirun -np ${nodes*tasks_per_node} \
               --map-by ppr:${tasks_per_node}:node \
               --bind-to none \
               "${target.get_install_binpath(case)}")
% else:
    (set -x; ${profiler} "${target.get_install_binpath(case)}")
% endif

${helpers.run_epilogue(target)}

echo
% endfor

${helpers.template_epilogue()}