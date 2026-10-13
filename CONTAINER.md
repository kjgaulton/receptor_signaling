# Run the R pipeline in a container

The image contains R, glmmTMB, Matrix, SeuratObject and SingleCellExperiment.
No host R installation or Python environment is required. Input data are mounted
at runtime; the image contains no study data. The earlier Python screen is not
included in this image.

## Build

Install/start Docker Desktop on macOS/Windows, or Docker Engine on Linux:

```bash
docker build -t receptor-signaling:0.1.0 .
docker run --rm receptor-signaling:0.1.0 --self-test
sh tests/test_container.sh receptor-signaling:0.1.0
```

The image targets Linux; Docker Desktop supplies the Linux VM on macOS/Windows.
The base image supports x86-64 (`linux/amd64`) and ARM64 (`linux/arm64`, including
Apple Silicon). Build on the destination architecture or use Buildx/emulation.
The initial source-package compilation can take several minutes and substantial
RAM; subsequent builds reuse the dependency layer. An Internet connection is
required for building, but not for running analyses.

Every build runs the synthetic model tests and RDS reader/CLI checks as a non-root
user. A failing test fails the build. Check installed versions with:

```bash
docker run --rm receptor-signaling:0.1.0 --versions
```

R 4.5.0 and the base image manifest digest are pinned. CRAN uses the fixed
2025-06-01 snapshot; Bioconductor uses release 3.21. Apt packages and maintenance
updates within that Bioconductor release are not fully locked. Each built image
records its actual package versions. For exact redistribution, transfer the built
image or use its registry digest rather than rebuilding at a later date.

## Configure and run (macOS/Linux shell)

Copy `container/config.txt` to `run-config.txt`. Edit receptor names, cell type,
metadata column names and the input RDS filename. Use container paths:
`/input/...` for data and `/output/<new-analysis-name>` for the output.

```bash
mkdir -p results
docker run --rm --network none \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$PWD/data,target=/input,readonly" \
  --mount "type=bind,source=$PWD/run-config.txt,target=/config.txt,readonly" \
  --mount "type=bind,source=$PWD/results,target=/output" \
  receptor-signaling:0.1.0 /config.txt
```

Create the host `data` directory and put your RDS there before running. The image
runs as UID/GID 10001 by default; overriding with your host IDs ensures output
ownership matches your account on Linux. The mounted output parent must be writable.
The `output` subdirectory specified in the configuration must not exist yet.
Configurations use plain-text `key = value` lines; they are not executed as R code.
For gene sets, set `gene_sets_gmt = /input/pathways.gmt` and place the GMT file
in the mounted input directory. Relative paths resolve from the config file directory.

## Windows PowerShell

Use Docker Desktop in Linux-container mode. Create `data`, `results` and
`run-config.txt` in the current directory, then:

```powershell
docker run --rm --network none `
  --mount "type=bind,source=$($PWD.Path)/data,target=/input,readonly" `
  --mount "type=bind,source=$($PWD.Path)/run-config.txt,target=/config.txt,readonly" `
  --mount "type=bind,source=$($PWD.Path)/results,target=/output" `
  receptor-signaling:0.1.0 /config.txt
```

## Transfer to another computer

No registry publication is required. Build for the receiving computer's CPU:

```bash
docker buildx build --platform linux/amd64 --load -t receptor-signaling:0.1.0 .
docker save -o receptor-signaling-amd64.tar receptor-signaling:0.1.0
# Copy the tar to another amd64 machine, then:
docker load -i receptor-signaling-amd64.tar
```

Use `linux/arm64` for ARM64 machines. Loading an image does not change its CPU
architecture. Docker/Podman emulation support is environment-dependent.
Do not copy the host `.Rlib` directory into an image; macOS libraries cannot run
inside Linux containers.

## HPC with Apptainer

On a Linux machine with the same CPU architecture as the image:

```bash
apptainer build receptor-signaling.sif docker-archive://receptor-signaling-amd64.tar
apptainer run --cleanenv \
  --bind "$PWD/data:/input:ro" \
  --bind "$PWD/run-config.txt:/config.txt:ro" \
  --bind "$PWD/results:/output" \
  receptor-signaling.sif /config.txt
```

Apptainer runs as your host user. Use the cluster's job scheduler and allocate
enough memory for the selected count matrix and per-gene model fits. The default
thread settings use one CPU thread; genes are still fitted sequentially.
Podman can also build/run this Dockerfile; use `--userns=keep-id` for output
ownership with rootless Podman, and configure SELinux bind labels where required.
These alternative-runtime commands should be validated on your cluster.

## Development

`container/entrypoint.sh` provides help, an example configuration, package versions
and the self-test. `.dockerignore` allows only runtime source and documentation
into the build context, excluding local libraries, Git history, raw data and results.
The analysis results and inference assumptions are described in `README.md`.

References: [Rocker platform support](https://rocker-project.org/images/versioned/r-ver.html),
[Docker multi-platform builds](https://docs.docker.com/build/building/multi-platform/).
