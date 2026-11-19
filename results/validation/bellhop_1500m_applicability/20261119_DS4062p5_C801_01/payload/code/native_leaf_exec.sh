#!/bin/bash
set -euo pipefail
printf '%s\n' "$$" > pid.txt
cat "/proc/$$/stat" > process_identity.txt
exec "$1" "$2"