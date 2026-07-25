#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=${1:-/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725}
TARBALL="$PROJECT_ROOT/software/downloads/ADFRsuite_x86_64Linux_1.0.tar.gz"
EXPECTED_SIZE=103496717
EXPECTED_SHA256=66539864ffb7ba87728d52ee50dde78586f551c7ca96aa455dae761875faf7e4
SOURCE_ROOT="$PROJECT_ROOT/software/src"
SOURCE_DIR="$SOURCE_ROOT/ADFRsuite_x86_64Linux_1.0"
INSTALL_DIR="$PROJECT_ROOT/software/ADFRsuite-1.0"
LOG="$PROJECT_ROOT/logs/adfrsuite_install_20260726.log"

mkdir -p "$SOURCE_ROOT" "$PROJECT_ROOT/logs"
exec > >(tee -a "$LOG") 2>&1
printf '%s START tarball=%s install_dir=%s\n' "$(date -Is)" "$TARBALL" "$INSTALL_DIR"

actual_size=$(stat -c %s "$TARBALL")
if [ "$actual_size" != "$EXPECTED_SIZE" ]; then
  printf '%s FAIL_SIZE expected=%s actual=%s\n' "$(date -Is)" "$EXPECTED_SIZE" "$actual_size"
  exit 2
fi
actual_sha=$(sha256sum "$TARBALL" | awk '{print $1}')
if [ "$actual_sha" != "$EXPECTED_SHA256" ]; then
  printf '%s FAIL_SHA256 expected=%s actual=%s\n' "$(date -Is)" "$EXPECTED_SHA256" "$actual_sha"
  exit 3
fi
tar -tzf "$TARBALL" >/dev/null

if [ -e "$SOURCE_DIR" ] || [ -e "$INSTALL_DIR" ]; then
  printf '%s FAIL_EXISTING_PATH source=%s install=%s\n' "$(date -Is)" "$SOURCE_DIR" "$INSTALL_DIR"
  exit 4
fi

tar -xzf "$TARBALL" -C "$SOURCE_ROOT"
(
  cd "$SOURCE_DIR"
  ./install.sh -d "$INSTALL_DIR" -c 0
)

required=(autosite autogrid4 prepare_receptor)
for command_name in "${required[@]}"; do
  if [ ! -x "$INSTALL_DIR/bin/$command_name" ]; then
    printf '%s FAIL_MISSING_EXECUTABLE path=%s\n' "$(date -Is)" "$INSTALL_DIR/bin/$command_name"
    exit 5
  fi
done

"$INSTALL_DIR/bin/autosite" --help >/dev/null 2>&1 || {
  rc=$?
  printf '%s WARN_AUTOSITE_HELP_EXIT rc=%s\n' "$(date -Is)" "$rc"
}
printf '%s PASS sha256=%s\n' "$(date -Is)" "$actual_sha"
