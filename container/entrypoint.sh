#!/bin/sh
set -eu
root=/opt/receptor-signaling
case "${1:---help}" in
  --help|-h)
    cat <<'EOF'
Usage: receptor-signaling /path/to/config.txt
       receptor-signaling --self-test
       receptor-signaling --versions
       receptor-signaling --example-config

Mount raw RDS data read-only, a configuration file, and a writable output parent.
Paths inside config.txt refer to paths INSIDE the container.
The configured output directory must not exist yet.
See CONTAINER.md for Docker, Podman and Apptainer examples.
EOF
    ;;
  --self-test)
    cd "$root"
    exec Rscript --vanilla tests/test_nb.R
    ;;
  --versions)
    cat "$root/packages.csv"
    ;;
  --example-config)
    cat "$root/container/config.txt"
    ;;
  *)
    exec Rscript --vanilla "$root/scripts/run_nb.R" "$@"
    ;;
esac
