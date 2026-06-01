# Cloud GPU Compute Plan — SERT Molecular Dynamics
**Fabian A. Alvarez-Primo, Ph.D.**
*Targeting FUTURE_WORK Phases 2–4*

---

## Platform Recommendation: AWS over GCP

Both work. AWS wins on three criteria:
- GROMACS AMIs and Docker images are more mature on AWS
- Spot Instances are cheaper and better documented than GCP Spot VMs for burst HPC
- S3 + EC2 data transfer within the same region is free

GCP is a viable fallback if AWS spot capacity is unavailable in your region.

---

## Instance Selection

| Phase | Workload | Recommended Instance | vCPUs | GPU | GPU RAM | On-Demand $/hr | Spot $/hr (est.) |
|-------|----------|----------------------|-------|-----|---------|----------------|-----------------|
| 2–3 | Ensemble docking, MM-GBSA | `g4dn.xlarge` | 4 | T4 | 16 GB | ~$0.53 | ~$0.16 |
| 4 | 200ns production MD (bilayer) | `p3.2xlarge` | 8 | V100 | 16 GB | ~$3.06 | ~$0.92 |
| 4 (faster) | Same, parallel systems | `p3.8xlarge` | 32 | 4x V100 | 64 GB | ~$12.24 | ~$3.67 |

**Primary recommendation: `p3.2xlarge` Spot for Phase 4.**
V100 is the gold standard for GROMACS — well-benchmarked, CUDA-optimized, 16GB handles a 150k-atom bilayer system comfortably.

**Cost estimate for Phase 4:**
- 200 ns/system × 2 systems = 400 ns total
- V100 throughput for SERT bilayer: ~40–60 ns/day
- Runtime: ~7–10 days of compute per system = 14–20 days total
- Spot cost: ~$0.92/hr × 24hr × 20 days = **~$440 worst case**
- In practice: run both systems sequentially or overlap = **$200–300 realistic**

---

## Environment Setup

### Step 1 — AWS Account Setup
```bash
# Install AWS CLI locally (WSL2)
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip && sudo ./aws/install
aws configure  # enter Access Key, Secret Key, region (us-east-1 recommended)
```

### Step 2 — Launch Instance (Spot)
Use the AWS Console or CLI:
```bash
aws ec2 run-instances \
  --image-id ami-0abcdef1234567890 \  # Deep Learning AMI (Ubuntu 22.04)
  --instance-type p3.2xlarge \
  --key-name your-keypair \
  --instance-market-options MarketType=spot \
  --block-device-mappings DeviceName=/dev/sda1,Ebs={VolumeSize=100} \
  --security-group-ids sg-xxxxxxxx
```

**Recommended AMI:** AWS Deep Learning AMI (Ubuntu 22.04) — comes with CUDA 12, cuDNN, conda pre-installed. GROMACS install becomes a single conda command.

### Step 3 — Install GROMACS on Instance
```bash
ssh -i your-keypair.pem ubuntu@<instance-ip>

# Via conda (easiest, GPU-enabled build)
conda create -n gromacs_env python=3.10
conda activate gromacs_env
conda install -c bioconda -c conda-forge gromacs

# Verify GPU detection
gmx mdrun --version  # should show CUDA device
```

---

## Data Transfer Workflow

### Local → S3 → Instance
```bash
# From WSL2 — upload your project
aws s3 mb s3://sert-md-project
aws s3 sync ./sert_project/ s3://sert-md-project/input/

# On instance — pull from S3 (free, same region)
aws s3 sync s3://sert-md-project/input/ ~/sert_project/
```

### Files to upload:
- `5I6Z_WT_prepared.pdb` and `5I6Z_S438T_prepared.pdb`
- Escitalopram ligand files (`.pdb`, `.mol2`)
- Your existing docking configs and Python scripts
- FUTURE_WORK.md as your run checklist

### Instance → S3 → Local (after run)
```bash
# On instance — push trajectories
aws s3 sync ~/sert_project/output/ s3://sert-md-project/output/

# Local — pull results
aws s3 sync s3://sert-md-project/output/ ./sert_results/
```

---

## Phase 4 Run Protocol (CHARMM-GUI → GROMACS)

### A. System Preparation (local, free web tools)
1. **CGenFF** (https://cgenff.umaryland.edu) — parameterize escitalopram, download `.str` file
2. **CHARMM-GUI Membrane Builder** (https://charmm-gui.org)
   - Input: your prepared SERT PDB
   - Membrane: POPC, 100 lipids/leaflet
   - Solvent: TIP3P, 150 mM NaCl
   - Output: download GROMACS-ready input files (`.gro`, `.top`, `.mdp`)
3. Upload CHARMM-GUI output package to S3

### B. Equilibration on Instance (~1–2 days)
```bash
# Standard GROMACS membrane equilibration sequence
gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 8 -gpu_id 0

gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -ntmpi 1 -ntomp 8 -gpu_id 0

gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt -ntmpi 1 -ntomp 8 -gpu_id 0
```

### C. Production MD (~7–10 days/system)
```bash
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
# Run with nohup so it survives SSH disconnection
nohup gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8 -gpu_id 0 \
  -cpt 30 > md.log 2>&1 &
```

**Critical:** Use `-cpt 30` (checkpoint every 30 min). Spot instances can be interrupted — checkpointing lets you resume from last checkpoint, not restart.

### D. Spot Interruption Protection
```bash
# Monitor for interruption notice (AWS sends 2-min warning)
# Add this to a cron job on the instance:
curl -s http://169.254.169.254/latest/meta-data/spot/termination-time
# Returns empty if safe, timestamp if terminating
# Script: if termination detected, sync output to S3 immediately
```

---

## MM-PBSA Rescoring (Phase 3/4)

After production MD, run MM-GBSA on the trajectory:
```bash
# AmberTools via conda (install alongside GROMACS env)
conda install -c conda-forge ambertools

# Convert GROMACS trajectory to Amber format
# Use gmx_MMPBSA (Python tool, bridges GROMACS + AmberTools)
pip install gmx_MMPBSA

gmx_MMPBSA -O -i mmpbsa.in \
  -cs md.tpr -ci index.ndx \
  -cg 1 13 -ct md_100ns.xtc \
  -o FINAL_RESULTS_MMPBSA.dat
```

---

## Cost Control Strategies

1. **Always use Spot** — 70% cheaper than On-Demand for interruptible workloads
2. **Set a Spot price cap** — bid ~2x current spot price to avoid frequent interruption
3. **Stop vs Terminate** — Stop the instance between work sessions (EBS volume persists); only Terminate when completely done
4. **S3 as primary storage** — EBS costs ~$0.10/GB/month; push trajectories to S3 ($0.023/GB/month) and shrink your EBS volume
5. **Set a billing alarm** — AWS Console → Billing → Alerts → $50 threshold. Non-negotiable.

---

## Realistic Cost Summary

| Phase | Instance | Est. Compute Time | Est. Spot Cost |
|-------|----------|-------------------|----------------|
| 2 — Ensemble docking | g4dn.xlarge | 2–3 days | ~$10–15 |
| 3 — MM-GBSA | g4dn.xlarge | 1–2 days | ~$5–8 |
| 4 — WT MD (200ns) | p3.2xlarge | 7–10 days | ~$155–220 |
| 4 — S438T MD (200ns) | p3.2xlarge | 7–10 days | ~$155–220 |
| **Total** | | | **~$325–465** |

**Minimum viable paper (Phases 1–3 only): ~$15–25 total.**

---

## Local Role (Your Machine Stays Useful)

Your GTX 1070 / Ryzen 7 1700X handles:
- All data preparation (PDB editing, CHARMM-GUI input prep)
- GROMACS trajectory analysis (cpptraj, MDAnalysis) — CPU-bound
- Figure generation (Python, matplotlib, VMD/ChimeraX visualization)
- MM-GBSA result interpretation
- All writing

Cloud handles only the GPU-intensive production MD. Your local machine stays the analysis and writing environment throughout.

---

## First Action Items

1. Create AWS account + enable billing alarm ($50 threshold)
2. Launch `g4dn.xlarge` Spot, install GROMACS via conda, run a 1 ns test simulation to verify the pipeline end-to-end
3. Submit WT system to CHARMM-GUI membrane builder while test runs
4. Upload CHARMM-GUI output to S3, pull to instance, begin equilibration

**Phase 2 and 3 can start within one weekend. Phase 4 follows immediately after.**
