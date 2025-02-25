import os
import pkgutil

# Path to your project
PROJECT_PATH = "/Users/sebastienjanel/Pasteur/Developpements/pyAF/pyAF"

# Set to store imported modules
imported_modules = set()

# Recursively scan Python files and extract imports
for root, _, files in os.walk(PROJECT_PATH):
    for file in files:
        if file.endswith(".py"):
            with open(os.path.join(root, file), "r", errors="ignore") as f:
                for line in f:
                    if line.startswith("import ") or line.startswith("from "):
                        parts = line.split()
                        if "import" in parts:
                            imported_modules.add(parts[1].split(".")[0])

# Find which installed packages match the imports
installed_packages = {pkg.name for pkg in pkgutil.iter_modules()}
used_packages = imported_modules & installed_packages

print("Packages needed for your project:", sorted(used_packages))
