#!/bin/sh
# Adds pyAF to the applications menu of the current user, from the folder of
# this script (no administrator rights needed). If the folder is moved, run it
# again.
#
#     ./install.sh               create the launcher
#     ./install.sh --uninstall   remove the launcher (not the application)

set -e

applications="${XDG_DATA_HOME:-$HOME/.local/share}/applications"
launcher="$applications/pyaf.desktop"

refresh() {
    if command -v update-desktop-database >/dev/null 2>&1; then
        update-desktop-database "$applications" >/dev/null 2>&1 || true
    fi
}

case "${1:-}" in
    --uninstall)
        rm -f "$launcher"
        refresh
        echo "pyAF removed from the applications menu."
        exit 0
        ;;
    "")
        ;;
    *)
        echo "Usage: $0 [--uninstall]" >&2
        exit 2
        ;;
esac

folder=$(cd "$(dirname "$0")" && pwd)
executable="$folder/pyAF"
icon="$folder/_internal/images/icon.png"

if [ ! -x "$executable" ]; then
    echo "pyAF not found in $folder: run this script from the pyAF folder." >&2
    exit 1
fi

# Exec key of the desktop entry: quoted, with its reserved characters escaped
exec_path=$(printf '%s' "$executable" |
    sed -e 's/[\\"`$]/\\\\&/g' -e 's/%/%%/g')

mkdir -p "$applications"
cat > "$launcher" <<DESKTOP
[Desktop Entry]
Type=Application
Name=pyAF
Comment=AFM force curve analysis
Exec="$exec_path"
Icon=$icon
Terminal=false
Categories=Science;
StartupWMClass=pyAF
DESKTOP
refresh

echo "pyAF added to the applications menu ($launcher)."
