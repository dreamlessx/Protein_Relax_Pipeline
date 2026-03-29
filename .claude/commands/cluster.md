Manage ACCRE and Curium cluster operations.

## Arguments
- $ARGUMENTS: Subcommand. Options:
  - "jobs" - List running/pending SLURM jobs
  - "quota" - Check disk usage and quotas
  - "logs <job_id>" - Tail logs for a specific job
  - "cancel <job_id>" - Cancel a SLURM job (asks confirmation)
  - "ssh" - Show SSH connection commands
  - "sync <direction>" - Sync data between local and ACCRE
  - Empty: show cluster overview

## Cluster Access

| Cluster | Command | VPN Required |
|---------|---------|-------------|
| ACCRE | `ssh accre` | No |
| Curium | `ssh curium` | Yes |

## Steps

### jobs
```bash
ssh accre "squeue -u agarwm5 --format='%.10i %.25j %.8T %.12M %.6D %.20R %.10l' | head -40"
```

### quota
```bash
ssh accre "echo '=== Home ===' && du -sh ~ && echo '=== Data ===' && du -sh /data/p_csb_meiler/agarwm5/ && echo '=== CSB Tmp ===' && du -sh /csbtmp/agarwm5/ 2>/dev/null"
```

### logs
```bash
ssh accre "tail -50 /data/p_csb_meiler/agarwm5/protein_pipeline/logs/<pattern>"
```

### cancel
Confirm with user first, then:
```bash
ssh accre "scancel <job_id>"
```

### sync (local -> ACCRE)
```bash
rsync -avz --progress scripts/ accre:/data/p_csb_meiler/agarwm5/af_work/scripts/
```

### sync (ACCRE -> local)
```bash
rsync -avz --progress accre:/data/p_csb_meiler/agarwm5/af_work/analysis/ analysis/
```

## Key Paths on ACCRE
```
Pipeline working dir:  /data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/
Blue pipeline data:    /data/p_csb_meiler/agarwm5/af_work/
Logs:                  /data/p_csb_meiler/agarwm5/protein_pipeline/logs/
AF2 install:           /sb/apps/alphafold232/
AF2 databases:         /csbtmp/alphafold-data.230/
Rosetta:               /data/p_csb_meiler/apps/rosetta/rosetta-3.15/
```

## SLURM Accounts
- `csb_gpu_acc` - AF GPU jobs (A6000)
- `p_meiler_acc` - Boltz GPU jobs (L40S), standalone AMBER
- `p_csb_meiler` - Rosetta CPU batch jobs
